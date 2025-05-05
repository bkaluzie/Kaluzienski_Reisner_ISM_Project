###################################################################################################
#
# Authors: Brooke Kaluzienski and Charles Reisner
#
# Created for use with the ULULA code.
#
###################################################################################################

import numpy as np
import matplotlib.pyplot as plt

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
        
        self.rho_inf = 1.0
        self.P_inf = 1.0
        self.g = 1.0
        self.gamma = 3.0 / 5.0
        
        self.cs_inf = np.sqrt(self.gamma * self.P_inf / self.rho_inf)
        
        self.Mdot = ((np.pi * self.rho_inf) / self.cs_inf**3) * (2 / (5 - 3 * self.gamma))**((5 - 3*self.gamma) / (2*self.gamma - 2))

        return
    
    # ---------------------------------------------------------------------------------------------
    def shortName(self):
        
        return 'bondi_parker_1d'

    # ---------------------------------------------------------------------------------------------
    
    def setInitialData(self, sim, nx):

        sim.setEquationOfState(eos_mode = 'ideal', gamma = self.gamma)
        sim.setGravityMode(gravity_mode = 'fixed_acc', g = self.g)
        sim.setDomain(nx, 1, xmin = 0, xmax = 2.0, bc_type = 'wall')
            
        DN = sim.q_prim['DN']
        PR = sim.q_prim['PR']
        VX = sim.q_prim['VX']

        x, _ = sim.xyGrid()

        sim.V[DN] = self.rho_inf
        sim.V[PR] = self.P_inf
        sim.V[VX] = 0
        
        sim.setGravityPotentials()

        return
    
    def u(self, sim, r):
        cs_r = sim.soundSpeed(sim.V)
        
        return (self.Mdot / (4.0 * np.pi * r**2)) * (cs_r / self.cs_inf)**(-2 / (self.gamma - 1.0))
    
    def boundaryConditions(self, sim):
        xhi = sim.xhi
        
        sim.V[sim.DN, xhi:, 0] = self.rho_inf
        sim.V[sim.PR, xhi:, 0] = self.P_inf
        sim.V[sim.VX, xhi:, 0] = 0.0
        
        sim.primitiveToConserved(sim.V[:, xhi:, :], sim.U[:, xhi:, :])
    
    '''
    def trueSolution(self, sim, x, q_plot, plot_geometry):

        xlo = sim.xlo
        xhi = sim.xhi
        
        sol_list = []
        for i in range(len(q_plot)):
            q = q_plot[i]
            sol = np.zeros((len(x)), float)
            if q == 'VX':
                sol = self.u(sim, x)[xlo:xhi+1]
            else:
                sol = None
            sol_list.append(sol)
        
        return sol_list
    '''
        
# ---------------------------------------------------------------------------------------------

def sonicSolutions():
    rho_inf = 1.0
    P_inf = 1.0
    gamma = 3.0 / 5.0
    cs_inf = np.sqrt(gamma * P_inf / rho_inf)
    
    r_B = (5 - 3 * gamma) / (4 * cs_inf**2)
    rho_B = rho_inf * (2 / (5 - 3*gamma))**(1 / (gamma - 1))
    cs_B = cs_inf * np.sqrt(2 / (5 - 3*gamma))
    
    return r_B, rho_B, cs_B
    

# ---------------------------------------------------------------------------------------------

def plotCallback(sim, fig, panels, plot_type):
    
    r_B, rho_B, cs_B= sonicSolutions()
    
    for panel in panels:
        panel.axvline(r_B)
        
    panels[0].axhline(rho_B)
    panels[2].axhline(cs_B)
    
    return

# ---------------------------------------------------------------------------------------------

def runBondi():
    setup = SetupBondiParker()
    ulula_run.run(setup, tmax = 5, nx = 200, q_plot = ['DN', 'PR', 'VX'], 
                  plot_time = 1, save_plots = False,
                  plot_callback_func = plotCallback)
    return