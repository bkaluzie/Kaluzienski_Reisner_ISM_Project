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
        
        self.gamma = 5.0 / 3.0
        self.g = 1.0
        self.rho_inf = 1.0
        self.P_inf = 1.0
        self.cs_inf = np.sqrt(self.gamma * self.P_inf / self.rho_inf)
        
        return
    
    # ---------------------------------------------------------------------------------------------
    def shortName(self):
        
        return 'bondi_parker'

    # ---------------------------------------------------------------------------------------------
    
   # need an update function to change the central mass (maybe?) 
   # start in uniform density, add material + velocity at the center to reinforce conditions
   # see how it moves out and how solution approaches parker
   # accretion -- set velocity at boundaries, let material flow in, update mass of object
   # can't add multiple velocities in single cell --- use gaussian (look at gresho vortex)
    
    def setInitialConditions(self, sim, nx):

        sim.setEquationOfState(eos_mode = 'ideal', gamma = self.gamma, pressure_floor = 0.01)
        sim.setGravityMode(gravity_mode = 'fixed_pot', g = self.g)
        sim.setDomain(nx, nx, xmin = -1.0, xmax = 1.0, ymin = -1.0, bc_type = 'outflow')
            
        DN = sim.q_prim['DN']
        PR = sim.q_prim['PR']
        VX = sim.q_prim['VX']
        VY = sim.q_prim['VY']
        GP = sim.q_prim['GP']

        x, y = sim.xyGrid()
        r = np.sqrt(x**2 + y**2)

        sim.V[DN] = self.rho_inf
        sim.V[PR] = self.P_inf
        sim.V[VX] = 0.0
        sim.V[VY] = 0.0
        
        sim.V[GP] = -0.1 / (0.01 + r) # potential from TDE setup

        return
    
    # ---------------------------------------------------------------------------------------------
    '''
    def u(self, sim, r):
        cs_r = sim.soundSpeed(sim.V)
        
        return (self.Mdot / (4.0 * np.pi * r**2)) * (cs_r / self.cs_inf)**(-2 / (self.gamma - 1.0))
    '''
    # ---------------------------------------------------------------------------------------------
    
    def updateFunction(self, sim):
        
        x, y = sim.xyGrid()
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        
        V_rad = 1.0 # inject some radial velocity
        V_x = V_rad * np.cos(phi)
        V_y = V_rad * np.sin(phi)
        
        sigma = 5.0 / sim.xhi
        
        zero_ind = int((sim.xhi + sim.xlo) / 2)
        
        rho_0 = sim.V[sim.DN, zero_ind + 6, zero_ind + 6]
        
        sim.V[sim.DN, zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] = rho_0 + 1.0 * np.exp(- 0.5 * r**2 / sigma**2)[zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] # gaussian density profile on top of whatever's there
        sim.V[sim.VX, zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] = V_x[zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5]
        sim.V[sim.VY, zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] = V_y[zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5]
        
        sim.primitiveToConserved(sim.V, sim.U)
        
        return

    # ---------------------------------------------------------------------------------------------
    '''
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
    '''
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

'''
def sonicSolutions():
    gamma = 5.0 / 3.0
    
    rho_inf = 1.0
    P_inf = 1.0
    cs_inf = np.sqrt(gamma * P_inf / rho_inf)
    
    q = 1/4 # for gamma = 5/3

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
'''
# ---------------------------------------------------------------------------------------------

def runBondi():
    setup = SetupBondiParker()
    ulula_run.run(setup, tmax = 15, nx = 200, q_plot = ['DN', 'PR', 'VR'], 
                  movie = True, movie_file_ext = 'gif')
    return
