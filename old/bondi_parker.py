###################################################################################################
#
# Authors: Brooke Kaluzienski and Charles Reisner
#
# Created for use with the ULULA code.
#
###################################################################################################

import numpy as np

import ulula.core.setup_base as setup_base
import ulula.core.run as ulula_run

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
        self.rho_inf = 1.0
        self.P_inf = 1.0
        self.g = 1.0

        return
    
    # ---------------------------------------------------------------------------------------------
    def shortName(self):
        
        return 'bondi_parker'

    # ---------------------------------------------------------------------------------------------
    
    def setInitialData(self, sim, nx):

        sim.setGravityMode(gravity_mode = 'fixed_acc', g = self.g)
        sim.setDomain(nx, nx, xmin = 0.0, xmax = 1.0, ymin = 0.0, bc_type = 'outflow')
            
        DN = sim.q_prim['DN']
        PR = sim.q_prim['PR']
        VX = sim.q_prim['VX']
        VY = sim.q_prim['VY']
        VR = np.sqrt(VX**2 + VY**2)

        x, y = sim.xyGrid()
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)

        sim.V[DN] = 1.0
        sim.V[PR] = 1.0
        


        def u(r):
            Mdot = 1.0
            gamma = 5/3
            cs_inf = np.sqrt(gamma * self.P_inf / self.rho_inf)
            cs_r = sim.soundSpeed(sim.V)

            return (Mdot / (4.0 * np.pi * r**2)) * (cs_r / cs_inf)**(-2 / (gamma - 1.0))
        

        sim.V[VX] = u(np.sqrt(x**2 + y**2)) * np.cos(phi)
        sim.V[VY] = u(np.sqrt(x**2 + y**2)) * np.sin(phi)
    
        return
    
    # ---------------------------------------------------------------------------------------------
    
    def boundaryConditions(self, sim):
        
        # stuff here
        
        return

def runBondi():
    setup = SetupBondiParker()
    ulula_run.run(setup, nx = 128, nsteps = 1000, dt = 0.01, plot_interval = 10, plot_vars = ['DN', 'PR', 'VX', 'VY'], plot_geometry = 'polar')
    return
