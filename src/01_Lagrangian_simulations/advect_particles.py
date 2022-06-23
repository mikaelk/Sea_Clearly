#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 12:01:53 2022

@author: kaandorp
"""
from parcels import FieldSet, ParticleSet, AdvectionRK4_3D, ErrorCode, ParticleFile, Field, \
    JITParticle, AdvectionRK4, DiffusionUniformKh, ParcelsRandom, VectorField
from parcels.tools.converters import Geographic, GeographicPolar 
import matplotlib.pyplot as plt
import sys
import numpy as np
import os
from glob import glob
import xarray as xr
import pandas as pd
from datetime import datetime, timedelta
from argparse import ArgumentParser
from scipy.interpolate import griddata

sys.path.insert(0,'..')
import settings
from kernels import unbeaching, delete_particle

from create_masks import make_landmask,get_coastal_nodes_diagonal,get_shore_nodes_diagonal,create_displacement_field
from utils.write_tools import to_netcdf


class Lagrangian_simulation:
    
    def __init__(self,settings):
        self.ufiles = sorted(glob( os.path.join(settings.DIR_INPUT,settings.DIR_UV,settings.PATTERN_U)))
        self.vfiles = sorted(glob( os.path.join(settings.DIR_INPUT,settings.DIR_UV,settings.PATTERN_V)))
        filenames = {'U': self.ufiles,
                     'V': self.vfiles}
        self.fieldset = FieldSet.from_netcdf(filenames, settings.VARS, settings.DIMS)
        data_sample_u = xr.load_dataset(self.ufiles[0])
        self.lons = data_sample_u['lon'].values
        self.lats = data_sample_u['lat'].values
        self.settings = settings 
        
    def set_landmask(self,to_fieldset=False):        
        outfile = os.path.join(settings.DIR_INPUT,settings.NAME_LANDMASK)
        if os.path.exists(outfile):
            ds = xr.open_dataset(outfile)
            self.landmask = np.array(ds['mask_land'],dtype=bool) 
        else:
            self.landmask = make_landmask(self.ufiles[0])
            to_netcdf(outfile,[self.landmask],['mask_land'],self.lons,self.lats,explanation='land mask')
    
        if to_fieldset:
            self.fieldset.add_field( Field('landMask',self.landmask,lon=self.lons,lat=self.lats,mesh='spherical') )
    
    def set_land_displacement(self,mag_u,to_fieldset=True):
        outfile = os.path.join(settings.DIR_INPUT,settings.NAME_LAND_U)
        if os.path.exists(outfile):
            ds = xr.open_dataset(outfile)
            self.u_land = np.array(ds['land_current_u'],dtype=float)    
            self.v_land = np.array(ds['land_current_v'],dtype=float)  
        else:
            self.set_landmask()
            self.u_land, self.v_land = create_displacement_field(self.landmask,mag_u)        
            to_netcdf(outfile,[self.u_land, self.v_land],['land_current_u','land_current_v'],self.lons,self.lats,explanation='land current, pusing particles on land back to the sea, magnitude of 1')
            
        if to_fieldset:
            U_unbeach = Field('U_unbeach',self.u_land,lon=self.lons,lat=self.lats,fieldtype='U',mesh='spherical')
            V_unbeach = Field('V_unbeach',self.v_land,lon=self.lons,lat=self.lats,fieldtype='V',mesh='spherical')
            # self.fieldset.add_field(U_unbeach)
            # self.fieldset.add_field(V_unbeach)
            vectorField_unbeach = VectorField('UV_unbeach',U_unbeach,V_unbeach)
            self.fieldset.add_vector_field(vectorField_unbeach)
        
    def set_turbulent_diffusion(self,K):
        K_m = K*np.ones([len(self.lats),len(self.lons)])
        K_z = K*np.ones([len(self.lats),len(self.lons)])

        Kh_meridional = Field('Kh_meridional', K_m,lon=self.lons,lat=self.lats,mesh='spherical')
        Kh_zonal = Field('Kh_zonal', K_z,lon=self.lons,lat=self.lats,mesh='spherical')

        self.fieldset.add_field(Kh_meridional)
        self.fieldset.add_field(Kh_zonal)
        
    def set_release(self,date_release):
                    
        if settings.RELEASE_MODE == 'uniform':
            dict_release = settings.DICT_RELEASE
            
            if dict_release['lon_min'] is None:
                dict_release['lon_min'] = min(self.lons)
            if dict_release['lon_max'] is None:
                dict_release['lon_max'] = max(self.lons)            
            if dict_release['lat_min'] is None:
                dict_release['lat_min'] = min(self.lats)            
            if dict_release['lat_max'] is None:
                dict_release['lat_max'] = max(self.lats)    
                
            lon_release = np.arange(dict_release['lon_min'],dict_release['lon_max'],dict_release['dlon'])
            lat_release = np.arange(dict_release['lat_min'],dict_release['lat_max'],dict_release['dlat'])
            
            
            if dict_release['remove_land']:
                self.set_landmask()
                X_release,Y_release = np.meshgrid(lon_release,lat_release)
                mesh_X,mesh_Y = np.meshgrid(self.lons,self.lats)
                land_mask_release = griddata((mesh_X.ravel(),mesh_Y.ravel()), self.landmask.ravel(), (X_release,Y_release), method='nearest').astype(bool)

                n_particles = len(X_release[~land_mask_release])
                lon_release = X_release[~land_mask_release].reshape([n_particles,1])
                lat_release = Y_release[~land_mask_release].reshape([n_particles,1])
                
            self.pset = ParticleSet.from_list(self.fieldset, JITParticle, 
                                         lon=lon_release,
                                         lat=lat_release,
                                         time=date_release.to_datetime64())
        else:
            raise NotImplementedError('Only uniform release implemented for now')
        
    def set_kernels(self,list_kernels):
        self.kernels = self.pset.Kernel(list_kernels[0])
        
        if len(list_kernels) > 1:
            for kernel_ in list_kernels[1:]:
                self.kernels += self.pset.Kernel(kernel_)
        
    def execute(self,date_start,date_end,dt_write,output_filename):
        pfile = ParticleFile(os.path.join(settings.DIR_OUTPUT,output_filename), self.pset, outputdt=timedelta(days=dt_write))
        
        if date_end > date_start:
            print('Running forwards simulation...')
            dt = timedelta(minutes=20)
        else:
            print('Running backwards simulation...')
            dt = timedelta(minutes=-20)
            
        runtime = abs(date_end - date_start)
        
        self.pset.execute(self.kernels,runtime=runtime,dt=dt,output_file=pfile,
                          verbose_progress=True,recovery={ErrorCode.ErrorOutOfBounds: delete_particle, ErrorCode.ErrorInterpolation: delete_particle})
        
        
if __name__ == "__main__":
    p = ArgumentParser(description="""Parcels runs to construct global transition matrices""")
    p.add_argument('-K_horizontal', '--K_horizontal', default=0, type=float, help='amount of horizontal diffusive mixing [m2/s], 0 for none')
    p.add_argument('-date_start', '--date_start', default='2014-01-01-12', type=str, help='Advection starting date')    
    p.add_argument('-date_end', '--date_end', default='2015-01-01-12', type=str, help='Advection end date')    
    p.add_argument('-dt_write', '--dt_write', default=1, type=float, help='Output dt (keep small for smooth particle simulation plotting)')    
    p.add_argument('-u_mag_land', '--u_mag_land', default=1, type=float, help='land current magnitude m/s')    
    
    args = p.parse_args()

    date_start = pd.Timestamp(args.date_start)
    date_end = pd.Timestamp(args.date_end)

    dt_write = args.dt_write   
    K_horizontal = args.K_horizontal
    u_mag_land = args.u_mag_land


    simulation = Lagrangian_simulation(settings)
    
    simulation.set_release(date_start)
    
    list_kernels = [AdvectionRK4]
    
    if u_mag_land > 0:
        simulation.set_land_displacement(u_mag_land,to_fieldset=True)
        list_kernels.append(unbeaching)
    
    if K_horizontal > 0:
        simulation.set_turbulent_diffusion(K_horizontal)
        list_kernels.append(DiffusionUniformKh)
                        
    simulation.set_kernels(list_kernels)

    output_filename = 'Test.nc'
    simulation.execute(date_start,date_end,dt_write,output_filename)