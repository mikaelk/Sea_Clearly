#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 15:06:46 2022

@author: kaandorp
"""
import math

def unbeaching(particle, fieldset, time):
    (vel_u, vel_v) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
    
    if math.fabs(vel_u) < 1e-14 and math.fabs(vel_v) < 1e-14:
        U_ub = fieldset.unbeach_U[time, particle.depth, particle.lat, particle.lon]
        V_ub = fieldset.unbeach_V[time, particle.depth, particle.lat, particle.lon]

        dlon = U_ub * particle.dt
        dlat = V_ub * particle.dt  

        particle.lon += dlon
        particle.lat += dlat
        
def delete_particle(particle, fieldset, time):
    """Kernel for deleting particles if they are out of bounds."""
    if fieldset.verbose_delete == 1:
        print('particle is deleted out of bounds at lon = ' + str(particle.lon) + ', lat =' + str(
            particle.lat) + ', depth =' + str(particle.depth))

    particle.delete()