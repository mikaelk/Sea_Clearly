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
else:
    raise RuntimeError('Check socket hostname, and add folder paths to src/settings.py file')


DIR_UV = 'CMEMS_MED/'

PATTERN_U = '[2011-2015]*RFVL*'
PATTERN_V = '[2011-2015]*RFVL*'
VARS = {'U': 'uo',
        'V': 'vo'}
DIMS = {'lat': 'lat',
        'lon': 'lon',
        'time': 'time'}

# land masks and displacement velocity fields are saved in DIR_INPUT for quick processing.
NAME_LANDMASK = 'landmask.nc'
NAME_LAND_U = 'land_displacement.nc'

# Release of particles
RELEASE_MODE = 'uniform'
# release bounds can be set to None to release from min/max lon and lat of the numerical data
DICT_RELEASE = {'lon_min': None,
                'lon_max': None,
                'lat_min': None,
                'lat_max': None,
                'dlon': 0.1,
                'dlat': 0.1,
                'remove_land': True}