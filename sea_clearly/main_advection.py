#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 12:05:23 2022

@author: kaandorp
"""

from argparse import ArgumentParser
import os,sys
from datetime import datetime, timedelta
import pandas as pd
import subprocess
from multiprocessing import Pool

module_path = os.path.abspath(os.path.join('.'))
if module_path not in sys.path:
    sys.path.append(module_path)
    
if __name__ == "__main__":
    p = ArgumentParser(description="""Run additional parcels simulations""")
    p.add_argument('-K_horizontal', '--K_horizontal', default=0, type=float, help='amount of horizontal diffusive mixing [m2/s], 0 for none')
    p.add_argument('-date_start_release', '--date_start_release', default='2015-01-01-12', type=str, help='Release starting date')    
    p.add_argument('-date_end_release', '--date_end_release', default='2015-02-01-12', type=str, help='Release end date')    
    p.add_argument('-n_days', '--n_days', default=365, type=int, help='Amount of advection days')    
    p.add_argument('-dt_write', '--dt_write', default=1, type=float, help='Output dt (keep small for smooth particle simulation plotting)')   
    p.add_argument('-u_mag_land', '--u_mag_land', default=1, type=float, help='land current magnitude m/s')    
    
    args = p.parse_args()

    date_start_release = pd.Timestamp(args.date_start_release)
    date_end_release = pd.Timestamp(args.date_end_release)
    
    if date_end_release < date_start_release: #for the pd.date_range to work, starting date should be smaller than end date
        assert args.n_days<0, 'Do you want a backwards simulation? Set n_days to a negative value'
        date_start_release, date_end_release = date_end_release, date_start_release
    
    dates_release = pd.date_range(date_start_release,date_end_release,freq='D') #daily release
    
    n_parallel = 10
    strings_execute = []
    for i1,date_start in enumerate(dates_release):
          
        print('Command:')
        str_cmd = ('python sea_clearly/advect_particles.py -date_start %s -n_days %s -K_horizontal %s -dt_write %s -u_mag_land %s' % 
                   (str(date_start).replace(' ','-'),args.n_days,
                    args.K_horizontal,args.dt_write,args.u_mag_land))
        print(str_cmd)
        strings_execute.append(str_cmd)
  
    def run_command(command):
        subprocess.call(command, shell=True)
    
    pool = Pool(processes=n_parallel)
    pool.map(run_command, strings_execute)

#         if i1 % n_parallel == (n_parallel-1):
#             procs = [ subprocess.Popen(i, shell=True, stdout=subprocess.PIPE) for i in strings_execute ]
#             for p in procs:
#                 for line in p.stdout: #the output of the first process is used as output for monitoring
#                     print(line)
#                 p.wait()
        
#             strings_execute = []
