# -*- coding: utf-8 -*-
"""
Authors: Brooke Kaluzienski and Charles Reisner

Created for use with the ULULA code.
"""
import numpy as np

import ulula.core.setup_base as setup_base

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
    param: type
        description
    
    """