#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 12:01:06 2022

@author: kaandorp
"""
import socket

if 'Kaandorp' in socket.gethostname():
    DIR_INPUT = '/Users/kaandorp/Data/SeaClearly/Input/'
    DIR_OUTPUT = '/Users/kaandorp/Data/SeaClearly/Output/'
elif 'lorenz' or 'node0' in socket.gethostname():
    DIR_INPUT = '/storage/shared/oceanparcels/output_data/data_Mikael/SeaClearly/Input/'
    DIR_OUTPUT = '/storage/shared/oceanparcels/output_data/data_Mikael/SeaClearly/Output/'
elif 'jupyter' or 'node0' in socket.gethostname():
    DIR_INPUT = '/home/jovyan/Data/Input/'
    DIR_OUTPUT = '/home/jovyan/Data/Output/'
else:
    raise RuntimeError('Check socket hostname, and add folder paths to src/settings.py file')


DIR_UV = 'CMEMS_MED'

PATTERN_U = '2018*MEDSEA*'
PATTERN_V = '2018*MEDSEA*'
VARS = {'U': 'uo',
        'V': 'vo'}
DIMS = {'lat': 'lat',
        'lon': 'lon',
        'time': 'time'}

# land masks and displacement velocity fields are saved in DIR_INPUT for quick processing.
NAME_LANDMASK = 'landmask.nc'
NAME_LAND_U = 'land_displacement.nc'

# Release of particles
# you don't have to change any of this if you're using pre-run data (as should proabably be the case)
RELEASE_MODE = 'uniform'
DICT_RELEASE = {'lon_min': None, #set uniform release bounds here, or set to None to use native grid bounds
                'lon_max': None,
                'lat_min': None,
                'lat_max': None,
                'n_gridcell': 2, #set amount of particles released per gridcell both in lon and lat dir (i.e. total = n_gridcell*n_gridcell). Prioritized over dlon/dlat below
                'dlon': 0.1, #spacing of uniform release grid. n_gridcell above is prioritized
                'dlat': 0.1,
                'remove_land': True} #use landmask to remove particles on land