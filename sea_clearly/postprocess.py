from glob import glob
import os,sys
import xarray as xr
import cartopy.crs as ccrs
from cartopy import feature
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
import pandas as pd
from IPython.display import HTML
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

from output_refinery import output_refinery_nearby
from create_masks import get_coastal_nodes
import settings

class postprocess:
    
    def __init__(self,settings,str_fwd,str_bwd):
        self.files_fwd = sorted(glob(os.path.join(settings.DIR_OUTPUT,str_fwd)))
        self.files_bwd = sorted(glob(os.path.join(settings.DIR_OUTPUT,str_bwd)))
        self.land_mask = xr.load_dataset(os.path.join(settings.DIR_INPUT,settings.NAME_LANDMASK))
        
    def refine(self,lon_refine,lat_refine,dlon=0.1,dlat=0.1,mode='fwd'):
        if mode == 'fwd':
            files = self.files_fwd
        elif mode == 'bwd':
            files = self.files_bwd
            
        for i1,file in enumerate(files):
            data = xr.open_zarr(file)
            if i1==0:
                data_refined = output_refinery_nearby(data,lon_refine,lat_refine,dlon=dlon/2,dlat=dlat/2)
            else:
                data_concat = output_refinery_nearby(data,lon_refine,lat_refine,dlon=dlon/2,dlat=dlat/2)
                data_refined = xr.concat([data_refined,data_concat],dim='traj')
        self.data_refined = data_refined
        self.n_release = len(files)
        
    def set_analysis_grids(self,dlon_coarse,dlat_coarse):
        X_land,Y_land = np.meshgrid(self.land_mask['lon'],self.land_mask['lat'])

        x_land_coarse = np.arange(self.land_mask['lon'].min(),self.land_mask['lon'].max()+dlon_coarse,dlon_coarse)
        y_land_coarse = np.arange(self.land_mask['lat'].min(),self.land_mask['lat'].max()+dlat_coarse,dlat_coarse)
        bins_x_land_coarse = np.insert( .5*(x_land_coarse[1:] + x_land_coarse[:-1]),[0,len(x_land_coarse)-1],[x_land_coarse[0]-.5*dlon_coarse,x_land_coarse[-1]+.5*dlon_coarse])
        bins_y_land_coarse = np.insert( .5*(y_land_coarse[1:] + y_land_coarse[:-1]),[0,len(y_land_coarse)-1],[y_land_coarse[0]-.5*dlat_coarse,y_land_coarse[-1]+.5*dlat_coarse])

        X_land_coarse,Y_land_coarse = np.meshgrid(x_land_coarse,y_land_coarse)
        self.analysis_landmask = griddata((X_land.ravel(),Y_land.ravel()),self.land_mask['mask_land'].values.ravel(),(X_land_coarse,Y_land_coarse),method='nearest')
        self.analysis_coastmask = get_coastal_nodes(self.analysis_landmask)
        self.analysis_X = X_land_coarse
        self.analysis_Y = Y_land_coarse
        self.analysis_bins_x = bins_x_land_coarse
        self.analysis_bins_y = bins_y_land_coarse
        self.X_pmesh,self.Y_pmesh = np.meshgrid(bins_x_land_coarse,bins_y_land_coarse)
        self.analysis_dlon = dlon_coarse
        self.analysis_dlat = dlat_coarse
        
    def calculate_likelihood(self,lon_refine,lat_refine,dlon=0.1,dlat=0.1,mode='fwd'):
        self.refine(lon_refine,lat_refine,dlon=dlon,dlat=dlat,mode=mode)
        
        hist = np.histogram2d(self.data_refined['lon'].values.ravel(),self.data_refined['lat'].values.ravel(),
                              bins=[self.analysis_bins_x,self.analysis_bins_y])[0].T
        
        ratio_analysis_release = (self.analysis_dlon*self.analysis_dlat) / (dlon*dlat)
        self.likelihood = (hist/hist.sum()) * (len(results.data_refined['obs']) / self.n_release) * ratio_analysis_release
        # self.likelihood = hist * ratio_analysis_release
        self.analysis_area_m2 = self.analysis_dlon*1.11e5 * (self.analysis_dlat*1.11e5*np.cos(lat_refine*(np.pi/180)))

        
    def calculate_priors(self,in_tot_coast=1525,in_tot_rivers=800): #in tonnes per year (see Kaandorp 2020)
        tyr_to_kgd = 1000/365
        
        input_coasts = xr.load_dataset('../data/coastal_MPW_input.nc')
        input_rivers = xr.load_dataset('../data/riverine_input.nc')

        X_coast,Y_coast = np.meshgrid(input_coasts['lon'],input_coasts['lat'])
        X_riv,Y_riv = np.meshgrid(input_rivers['lon'],input_rivers['lat'])

        mask_coast_coarse = self.analysis_coastmask == 1
        mask_input_coast = input_coasts['MPW'] > 0
        mask_input_rivers = input_rivers['input_rivers'] > 0
        
        vals_coast_coarse = griddata((X_coast[mask_input_coast],Y_coast[mask_input_coast]),
                                       input_coasts['MPW'].values[mask_input_coast],
                                       (self.analysis_X[mask_coast_coarse],self.analysis_Y[mask_coast_coarse]),method='nearest')
        
        vals_rivers_coarse = griddata((X_coast[mask_input_rivers],Y_coast[mask_input_rivers]),
                                       input_rivers['input_rivers'].values[mask_input_rivers],
                                       (self.analysis_X[mask_coast_coarse],self.analysis_Y[mask_coast_coarse]),method='nearest')
                
        self.input_coast = np.zeros(self.analysis_landmask.shape)
        self.input_coast[mask_coast_coarse] = vals_coast_coarse
        self.input_coast *= (in_tot_coast/self.input_coast.sum())*tyr_to_kgd
        
        self.input_rivers = np.zeros(self.analysis_landmask.shape)
        self.input_rivers[mask_coast_coarse] = vals_rivers_coarse
        self.input_rivers *= (in_tot_rivers/self.input_rivers.sum())*tyr_to_kgd  
        
        self.input_total = self.input_coast + self.input_rivers
        
    def calculate_posterior(self):
        self.posterior = self.likelihood * self.input_total
        self.concentration_g_m2 = (self.posterior.sum()*1000) / self.analysis_area_m2
          
    #method that can be used to sort results by time, which can be nice for animations but not used at the moment
    # at the moment, all release times are release at once
    def sort_results_by_time(self): 
        sim_fwd = results.data_refined['time'][0,1].values > results.data_refined['time'][0,0].values
        t_min = results.data_refined['time'].values.min()
        t_max = results.data_refined['time'].values.max()
        times_plot = pd.date_range(t_min,t_max,freq='D')
        if not sim_fwd:
            times_plot = times_plot[::-1]

        data_processed = {}
        for time_ in times_plot:
            data_processed[time_] = {}

            mask = (results.data_refined['time'] == time_).values
            data_processed[time_]['lons'] = results.data_refined['lon'].values[mask]
            data_processed[time_]['lats'] = results.data_refined['lat'].values[mask]

    def set_figure(self,xlim=(-6,36.5),ylim=(30,46),background_map='mesh'):
        self.fig = plt.figure(figsize=(14,4.5))
        gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[30,1])
        
        self.ax = self.fig.add_subplot(gs[0, 0],projection=ccrs.PlateCarree()) #plot with the map
        self.cax = self.fig.add_subplot(gs[0,1]) #plot with the colorbar
        self.cax.axis('off') #hide colorbar for now
    
        if background_map == 'mesh':
            self.ax.pcolormesh(self.analysis_X, self.analysis_Y, self.analysis_landmask[:-1,:-1],cmap='Greys',alpha=.5,zorder=0)
        else:
            self.ax.add_feature(feature.NaturalEarthFeature('physical','land','50m'),facecolor='grey',zorder=0)

        self.ax.set_xlim(xlim[0],xlim[1])
        self.ax.set_ylim(ylim[0],ylim[1])
    
    def update_fig_title(self,str_):
        if self.title:
            self.title.set_text(str_)
        else:
            self.title = self.ax.set_title(str_)
        self.fig.canvas.draw()
        
    def animate_particles(self,framedt=100):

        self.data_refined = self.data_refined.load()
        self.update_fig_title('Live integration of particle trajectories...')
        plt.pause(.5)
        self.points, = self.ax.plot([], [], 'o',markersize=2,zorder=1)     
        def animate(frame_num):
            self.points.set_data(( self.data_refined['lon'][:,frame_num],self.data_refined['lat'][:,frame_num]))
            
            if frame_num > n_start_opaque:
                n_opaque_frames = n_frames - n_start_opaque
                i_from_start = frame_num - n_start_opaque
                self.points.set_alpha(1-(i_from_start/n_opaque_frames))
            
            if frame_num == n_frames-1:
                self.update_fig_title('Animation done, click for histogram (likelihood)')

            return self.points
        
        n_frames = len(self.data_refined['obs'])-1
        n_start_opaque = n_frames - 10
        anim = FuncAnimation(self.fig, animate, frames=len(self.data_refined['obs'])-1, interval=framedt, repeat=False)
        return anim
    
    def analysis_and_plot(self,dlon_analysis,dlat_analysis,framedt=100,background_map='mesh'):
        # format log-scale colorbars
        def fmt(x, pos):
            float_ = 10**float(x)
            return r'%4.4f'%float_
     
    
        self.set_analysis_grids(dlon_analysis,dlat_analysis)
        self.set_figure(background_map=background_map)
        self.title = None
        
        reset_plot = False
        pos = []
        clicks = []

        def onclick(event):
            pos = [event.xdata,event.ydata]
            clicks.append(1)

            if len(clicks) == 1:
                self.pos_analysis = pos
                self.ax.plot(self.pos_analysis[0],self.pos_analysis[1],'kx',markersize=14,zorder=3)
                self.update_fig_title('Calculating particle statistics...')
                self.calculate_likelihood(pos[0],pos[1],dlon=dlon_analysis,dlat=dlat_analysis,mode='bwd')

                self.calculate_priors()
                self.calculate_posterior()
                
                self.update_fig_title('Big data server request (can take some time)...')
                global animation
                animation = self.animate_particles(framedt=framedt)
                animation.event_source.start()
                
            elif len(clicks) == 2:
                likelihood_plot = self.likelihood.copy()
                likelihood_plot[likelihood_plot == 0] = np.nan
                likelihood_plot = likelihood_plot * (100/np.nansum(likelihood_plot)) # convert to percentages
                int_min = int(np.floor(np.nanmin(np.log10(likelihood_plot))))
                int_max = int(np.ceil(np.nanmax(np.log10(likelihood_plot))))
                ticks = np.arange(int_min,int_max+1,1)
                self.plt_lik = self.ax.pcolormesh(self.analysis_X, self.analysis_Y, np.log10(likelihood_plot[:-1,:-1]),cmap=plt.cm.viridis,alpha=.7,zorder=2)
                self.cax.axis('on')
                cbar=self.fig.colorbar(self.plt_lik,cax=self.cax,ticks=ticks,format=ticker.FuncFormatter(fmt))
                cbar.set_label('Backtracking histogram [%]', rotation=90)
                self.update_fig_title('Click for posterior pdf (litter sources)...')
            
            elif len(clicks) == 3:
                posterior_plot = self.posterior.copy()
                posterior_plot[posterior_plot==0] = np.nan
                posterior_plot = posterior_plot * (100/np.nansum(posterior_plot)) # convert to percentages
                int_min = int(np.floor(np.nanmin(np.log10(posterior_plot))))
                int_max = int(np.ceil(np.nanmax(np.log10(posterior_plot))))
                ticks = np.arange(int_min,int_max+1,1)
                self.points.set_alpha(0.0)
                self.plt_lik.remove()
                self.plt_post = self.ax.pcolormesh(self.analysis_X, self.analysis_Y, np.log10(posterior_plot[:-1,:-1]),cmap=plt.cm.viridis,alpha=1.,zorder=2)
                cbar=self.fig.colorbar(self.plt_post,cax=self.cax,ticks=ticks,format=ticker.FuncFormatter(fmt))
                cbar.set_label('Pollution origin [%]', rotation=90)
                self.update_fig_title('Estimated concentration: %.2e g/m2' % self.concentration_g_m2)
                
            elif len(clicks) == 4:
                self.fig.canvas.mpl_disconnect(cid)
                reset_plot = True
                # stop_plotting = True
             
        cid = self.fig.canvas.mpl_connect('button_press_event', onclick)
