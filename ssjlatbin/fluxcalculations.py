# (C) 2021 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam Kilcommons - University of Colorado, Boulder - Colorado Center for Astrodynamics Research
# Mar 2021
import numpy as np
import pandas as pd
import datetime
from dateutil.relativedelta import *
from functools import partial

from logbook import Logger
log = Logger('fluxcalculations')

#Center energy of each channel in eV
CHANNEL_ENERGIES = [ 30000.,  20400.,  13900.,   9450.,   6460.,   4400.,   3000.,
                        2040.,   1392.,    949.,    646.,    440.,    300.,    204.,
                        139.,     95.,     65.,     44.,     30.]

def _hardy_channel_widths():
    """Calculate the 'width' (in energy) of each
    DMSP SSJ channel in order to integrate energy from each channel
    
    Hardy, D. A., Holeman, E. G., Burke, W. J., Gentile, L. C., 
    and Bounar, K. H. (2008), 
    Probability distributions of electron precipitation at 
    high magnetic latitudes, J. Geophys. Res., 113, A06305, 
    doi:10.1029/2007JA012746. 

    """
    channel_energies = CHANNEL_ENERGIES
    channel_width_energy = np.zeros((19,1))
    channel_width_energy[0] = channel_energies[0]-channel_energies[1]
    for k in range(1,18):
        channel_width_energy[k] = (channel_energies[k-1]-channel_energies[k+1])/2.

    channel_width_energy[18] = channel_energies[17]-channel_energies[18]
    return channel_width_energy

def _energy_or_number_check(energy_or_number):
    if energy_or_number not in ['energy','number']:
        raise ValueError('Invalid type of flux {}'.format(energy_or_number))

def _integrated_flux_std(diff_flux,diff_flux_rel_uncert,channels,energy_or_number):
    """Calculate the appoximate standard deviation of a full or partial
    integration (over channel numbers 'channels') 
    of the differential energy or number flux"""
    channel_energies = CHANNEL_ENERGIES
    channel_width_energy = _hardy_channel_widths()

    #Sigma is given as relative error, make absolute
    sigma = diff_flux_rel_uncert*diff_flux  #Standard deviation of each channel

    integral_flux_std = np.zeros((sigma.shape[0],1)).flatten()
    #Make it so a NaN doesn't NaN the entire computation
    for c in channels:
        nn = np.logical_not(np.isfinite(sigma[:,c]))
        sigma[nn,c] = 0. #Zero it out so it doesn't break the sum
        ccs = sigma[:,c]/channel_energies[c]*channel_width_energy[c]*np.pi# #/cm^2/s/sr/eV -> #/cm^2/s
        ccs *= 1e4 # #/cm^2/s -> #/m^2/s
        if energy_or_number == 'energy':
            ccs *= channel_energies[c]*1.6e-19 # #/m^2/s -> J/m^2/s == W/m^2
            ccs *= 1e3 #W/m^/s -> mW/m^2/s
        ccs = ccs**2. #Uncertainties add in quadrature
        integral_flux_std = integral_flux_std + ccs.flatten()
    return np.sqrt(integral_flux_std)

def _integrate_flux(diff_flux,channels,energy_or_number):
    """Calculate the total energy or number flux for a set of channels
    from the differential (channel-by-channel) energy flux
    """
    _energy_or_number_check(energy_or_number)
        
    channel_energies = CHANNEL_ENERGIES
    channel_width_energy = _hardy_channel_widths()

    #Integrate the flux across the chosen channels
    integral_flux = np.zeros((diff_flux.shape[0],1)).flatten()

    for c in channels:
        #Make it so a NaN doesn't NaN the entire computation
        nn = np.logical_not(np.isfinite(diff_flux[:,c]))
        diff_flux[nn,c] = 0. #Zero it out so it doesn't break the sum
        cc = diff_flux[:,c]/channel_energies[c]*channel_width_energy[c]*np.pi# #/cm^2/s/sr/eV -> #/cm^2/s
        cc *= 1e4 # #/cm^2/s -> #/m^2/s
        if energy_or_number == 'energy':
            cc *= channel_energies[c]*1.6e-19 # #/m^2/s -> J/m^2/s == W/m^2
            cc *= 1e3 #W/m^/s -> mW/m^2/s
        #'Mean contribution from ssj channel %d was %.3f mW/m^2' % (c,np.nanmean(cc)) )
        integral_flux = integral_flux + cc.flatten()
            
    return integral_flux

def integrate_flux(diff_flux,channels,energy_or_number,diff_flux_rel_uncert=None,uncertainty_tolerance=None):
    integral_flux = _integrate_flux(diff_flux,channels,energy_or_number)
    if diff_flux_rel_uncert is not None and uncertainty_tolerance is not None:
        integral_flux_std = _integrated_flux_std(diff_flux,
                                                diff_flux_rel_uncert,
                                                channels,
                                                energy_or_number)
        #uncertainty tolerance is a percent
        too_uncert = integral_flux_std/integral_flux*100>=uncertainty_tolerance
        integral_flux[too_uncert]=0.
        log.debug('Removed %d/%d flux values due to high uncertainty' % (
                                                np.count_nonzero(too_uncert),
                                                len(too_uncert)))
    else:
        log.debug('Uncertainty not checked')
    return integral_flux

def average_particle_energy(integrated_energy_flux,integrated_number_flux):
    """Calculate the average particle energy (eV) for output from above function 
    (only valid if same channel set used for energy and number flux)"""
    #Calculate average energy in electron-volts
    eavg = (integrated_energy_flux/1.6e-19/1000)/integrated_number_flux #mW/m^2->eV/m^2
    return eavg