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
        
        self.gamma = 7.0 / 5.0
        self.g = 1.0
    
        self.r_B = 0.5 # define the sonic radius
        
        self.cs_inf = np.sqrt((5.0 - 3.0 * self.gamma) / (4.0 * self.r_B))
        self.rho_inf = 1.0
        self.P_inf = self.cs_inf * self.rho_inf / self.gamma
        
        self.rho_B = self.rho_inf * (2 / (5 - 3*self.gamma))**(1 / (self.gamma - 1))
        self.P_B = self.P_inf  * (self.rho_B / self.rho_inf)**self.gamma
        self.cs_B = self.cs_inf * np.sqrt(2 / (5 - 3*self.gamma))
        
        self.Mdot = ((np.pi * self.rho_inf) / self.cs_inf**3) * (2 / (5 - 3 * self.gamma))**((5 - 3*self.gamma) / (2*self.gamma - 2))
        
        
        '''
        self.rho_inf = 1.0
        self.P_inf = 1.0
        self.g = 1.0
        self.gamma = 7.0 / 5.0
        
        self.cs_inf = np.sqrt(self.gamma * self.P_inf / self.rho_inf)
        
        # define bondi radius quantities
        self.r_B = (5 - 3 * self.gamma) / (4 * self.cs_inf**2)
        self.rho_B = self.rho_inf * (2 / (5 - 3*self.gamma))**(1 / (self.gamma - 1))
        self.cs_B = self.cs_inf * np.sqrt(2 / (5 - 3*self.gamma))
        
        self.Mdot = ((np.pi * self.rho_inf) / self.cs_inf**3) * (2 / (5 - 3 * self.gamma))**((5 - 3*self.gamma) / (2*self.gamma - 2))
        '''
        
        return
    
    # ---------------------------------------------------------------------------------------------
    def shortName(self):
        
        return 'bondi_parker_1d'

    # ---------------------------------------------------------------------------------------------
    
    def setInitialData(self, sim, nx):

        sim.setEquationOfState(eos_mode = 'ideal', gamma = self.gamma)
        sim.setGravityMode(gravity_mode = 'fixed_acc', g = self.g)
        sim.setDomain(nx, 1, xmin = 0, xmax = 1.0, bc_type = 'wall')
            
        DN = sim.q_prim['DN']
        PR = sim.q_prim['PR']
        VX = sim.q_prim['VX']

        x, _ = sim.xyGrid()
        
        def linear(x, a_x, a_y, b_x, b_y):
            m = (b_y - a_y) / (b_x - a_x)
            b = a_y - m*a_x
            return m*x + b

        # idk i just made all these linear -- solution u(r) does NOT work and isn't correct ?
        # tbf this also doesn't work.
        sim.V[DN] = linear(x, self.r_B, self.rho_B, x[sim.xhi], self.rho_inf)
        sim.V[PR] = linear(x, self.r_B, self.P_B, x[sim.xhi], self.P_inf)
        sim.V[VX] = linear(x, self.r_B, self.cs_B, x[sim.xhi], 0.0)
        
        sim.setGravityPotentials()

        return
    
    # ---------------------------------------------------------------------------------------------
    
    def u(self, sim, r):
        cs_r = sim.soundSpeed(sim.V)
        
        return (self.Mdot / (4.0 * np.pi * r**2)) * (cs_r / self.cs_inf)**(-2 / (self.gamma - 1.0))
    
    # ---------------------------------------------------------------------------------------------
    
    def boundaryConditions(self, sim):
        xhi = sim.xhi
        xlo = sim.xlo
        dx = sim.dx
        
        ind_B = int(xlo + self.r_B / dx)
        
        sim.V[sim.DN, xhi:, 0] = self.rho_inf
        sim.V[sim.PR, xhi:, 0] = self.P_inf
        sim.V[sim.VX, xhi:, 0] = 0.0 # velocity at infinity is 0
        
        # set sonic radius conditions
        sim.V[sim.DN, ind_B, 0] = self.rho_B
        sim.V[sim.PR, ind_B, 0] = self.P_B
        sim.V[sim.VX, ind_B, 0] = self.cs_B 
        
        sim.primitiveToConserved(sim.V[:, xhi:, :], sim.U[:, xhi:, :])
        
        return
    
    # ---------------------------------------------------------------------------------------------
    
    def plotLimits(self, q_plot, plot_geometry):
        
        vmin = []
        vmax = []

        for q in q_plot:
            if q == 'DN':
                vmin.append(0.0)
                vmax.append(20.0)
            elif q in ['VX']:
                vmin.append(0.0)
                vmax.append(10.0)
            elif q == 'PR':
                vmin.append(0.0)
                vmax.append(20.0)
            else:
                vmin.append(None)
                vmax.append(None)
        
        return vmin, vmax, None
        
        return
    
    # ---------------------------------------------------------------------------------------------
    
    ''' this is behaving VERY weirdly---several lines per plot???
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
    r_B = 0.5 # define the sonic radius
    gamma = 7.0 / 5.0
    
    cs_inf = np.sqrt((5.0 - 3.0 * gamma) / (4.0 * r_B))
    rho_inf = 1.0
    P_inf = cs_inf * rho_inf / gamma
    
    rho_B = rho_inf * (2 / (5 - 3*gamma))**(1 / (gamma - 1))
    P_B = P_inf  * (rho_B / rho_inf)**gamma
    cs_B = cs_inf * np.sqrt(2 / (5 - 3*gamma))
    
    return r_B, rho_B, P_B, cs_B
    

# ---------------------------------------------------------------------------------------------

def plotCallback(sim, fig, panels, plot_type):
    
    r_B, rho_B, P_B, cs_B = sonicSolutions()
    
    for panel in panels:
        panel.axvline(r_B)
        panel.set_xlabel('$r$ (CU)')
        
    panels[0].axhline(rho_B)
    panels[1].axhline(P_B)
    panels[2].axhline(cs_B)
    panels[2].set_ylabel('$v_r$ (CU)')
    
    return

# ---------------------------------------------------------------------------------------------

def runBondi():
    setup = SetupBondiParker()
    ulula_run.run(setup, tmax = 5, nx = 200, q_plot = ['DN', 'PR', 'VX'], 
                  plot_time = 0.1, save_plots = False,
                  plot_callback_func = plotCallback)
    return
