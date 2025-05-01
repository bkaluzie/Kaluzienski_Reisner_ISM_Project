###################################################################################################
#
# Authors: Brooke Kaluzienski and Charles Reisner
#
# Created for use with the ULULA code.
#
###################################################################################################

import numpy as np

import ulula.core.setup_base as setup_base

###################################################################################################

class SetupBondiParker(setup_base.Setup):
    
    """
    Spherical accretion/outflow setup
    
    In this setup, a star is placed at the center of a spherically-symmetric gas cloud.
    The star is at rest relative to distant regions of the cloud. Far from the star, 
    the cloud has a uniform density rho_inf and a uniform pressure P_inf. 
    This setup demonstrates
    
    * Ideal equation of state
    * Sub- and super-sonic flows
    * Polar coordinate solver.
    
    Parameters
    ----------
    unit_l: float
        Code unit for length in units of centimeters.
    unit_t: float
        Code unit for time in units of seconds.
    unit_m: float
        Code unit for mass in units of gram.
    
    """
    
    def __init__(self, unit_l = 1.0, unit_t = 1.0, unit_m = 1.0):
        
        setup_base.Setup.__init__(self, unit_l = unit_l, unit_t = unit_t, unit_m = unit_m)
        
        # stuff here
        
        return
    
    # ---------------------------------------------------------------------------------------------

    def shortName(self):
        
        return 'bondi_parker'

    # ---------------------------------------------------------------------------------------------
    
    def setInitialData(self, sim, nx):
        
        # stuff here
        
        return
    
    # ---------------------------------------------------------------------------------------------
    
    