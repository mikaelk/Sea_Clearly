import numpy as np
import geopandas as gpd
# import settings
import fiona
import pandas as pd
from shapely.geometry import Polygon
from shapely.geometry import Point
import xarray as xr


def to_netcdf(output_filename,data,data_name,lons,lats,explanation=''):
    """
    

    Parameters
    ----------
    output_filename : str
        output file path
    data : list
        list containing the data to write to netcdf
    data_name : list
        list containing the variable names to be written
    lons : np.array
        longitudes
    lats : np.array
        latitudes
    explanation : str
        explanation of the netcdf contents

    Returns
    -------
    None.

    """
    
    dict_data = {}
    for name_, data_ in zip(data_name, data):
        dict_data[name_] = (( "lat", "lon"), data_ )
    dict_data['explanation'] = explanation
        
    ds = xr.Dataset(
        dict_data,
        coords={
            "lon": lons,
            "lat": lats,
        },
    )   
    ds.to_netcdf(output_filename)
    

filedir = '/storage/shared/oceanparcels/output_data/data_Darshika/SeaClearlyStuff/'
gdb_file = filedir + "EMODnet_HA_Environment_Natura2000_end2020_20210909/EMODnet_HA_Environment_Natura2000_end2020_20210909.gdb" #EMODnet_HA_Environment_Natura2000_end2020_20210909/


top_lat, bottom_lat = 46, 30
left_lon, right_lon = -7, 37


full_data = gpd.read_file(gdb_file, bbox=(left_lon, bottom_lat, right_lon, top_lat))
data = full_data[full_data.COAST_MAR == 1]
shapes = data.geometry

landmask = xr.load_dataset('../data/CMEMS_MED_landmask.nc')

MPA_field = np.zeros(landmask['mask_land'].shape)

for i_lon, lon_ in enumerate(landmask['lon']):
    print(i_lon)
    for i_lan, lat_ in enumerate(landmask['lat']):
        for i_shape, shape_ in enumerate(shapes):
            for geom_ in shape_.geoms:
                contains = Polygon(geom_).contains(Point(lon_,lat_))
                if contains:
                    MPA_field[i_lan, i_lon] = i_shape
                    
                    
to_netcdf('../data/MPA_field.nc', [MPA_field], ['MPA_field'],
          landmask['lon'].values,
         landmask['lat'].values,
         explanation='Field containig the binned MPAs.')