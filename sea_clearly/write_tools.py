#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:05:22 2022

@author: kaandorp
"""
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