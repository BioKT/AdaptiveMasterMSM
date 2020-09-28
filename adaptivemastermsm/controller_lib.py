"""
This file is part of the AdaptiveMasterMSM package.

"""

#def driver_equilibration(run):
#
#    equil = set_equilibration()
#    equil['total_steps'] = \
#        int(equil['total_time'] / equil['timestep'])
#    for p in equil.keys():
#        run[p] = equil[p]
#    return run
#
#def driver_production(run):
#
#    prod = set_production()
#    prod['total_steps'] = \
#        int(prod['total_time'] / prod['timestep'])
#    for p in prod.keys():
#        run[p] = prod[p]
#    return run
#
#def set_production():
#    """
#    Parameters for the MD production run
#
#    """
#    production = {
#        'timestep'       : 0.002,         # ps
#        'total_time'     : 1,             # ps / 10ns
#        'log_freq'       : 5000,          # timesteps / 20ps
#        'xtc_freq'       : 5000,          # timesteps
#        'temperature'    : 300,           # Kelvin
#        'thermostat'     : 'V-rescale',   # type
#        'box_size'       : 1.0,           # Angstroms
#        'barostat'       : 'No',          # No barostat
#        'pressure'       : 1,             # pressure
#        'Cl'             : 0,             # number Cl's
#        'Na'             : 0              # number Na's
#    }
#    return production
#
#def set_equilibration():
#    """
#    Parameters for the MD equilibration run
#
#    """
#    equilibration = {
#        'timestep'       : 0.002,         # ps
#        'total_time'     : 1,             # ps
#        'barostat'       : 'berendsen',   # type
#        'pressure'       : 1,             # bar
#        'log_freq'       : 10000,         # timesteps / 20ps
#        'xtc_freq'       : 100,           # timesteps
#        'temperature'    : 300,           # Kelvin
#        'thermostat'     : 'Nose-Hoover', # type
#        'box_size'       : 1.0,           # Angstroms
#        'Cl'             : 0,             # number Cl's
#        'Na'             : 0              # number Na's
#    }
#    return equilibration

