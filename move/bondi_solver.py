###################################################################################################
#
# Authors: Brooke Kaluzienski and Charles Reisner
#
# Created for use with the ULULA code.
#
# written by Charles in an attempt to get a spherical inflow working until his dumbass can't figure it out
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
        self.rho_inf = 0.5
        self.P_inf = 0.5
        self.cs_inf = np.sqrt(self.gamma * self.P_inf / self.rho_inf)
        
        return
    
    # ---------------------------------------------------------------------------------------------
    def shortName(self):
        
        return 'bondi'

    # ---------------------------------------------------------------------------------------------
    
   # need an update function to change the central mass (maybe?) 
   # start in uniform density, add material + velocity at the center to reinforce conditions
   # see how it moves out and how solution approaches parker
   # accretion -- set velocity at boundaries, let material flow in, update mass of object
   # can't add multiple velocities in single cell --- use gaussian (look at gresho vortex)
    
    def setInitialConditions(self, sim, nx):

        sim.setEquationOfState(eos_mode = 'ideal', gamma = self.gamma, pressure_floor = 0.01)
        sim.setGravityMode(gravity_mode = 'fixed_pot', g = self.g)
        sim.setDomain(nx, nx, xmin = -3.0, xmax = 3.0, ymin = -3.0, bc_type = 'outflow')
            
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
        
        sim.V[GP] = -self.g / (0.03 + r) # potential from TDE setup

        return
    
    # ---------------------------------------------------------------------------------------------
    
    def u(self, sim, r):
        cs_r = sim.soundSpeed(sim.V)
        
        return -(self.Mdot / (4.0 * np.pi * r**2)) * (cs_r / self.cs_inf)**(-2 / (self.gamma - 1.0))
    
    # ---------------------------------------------------------------------------------------------
    
#    def updateFunction(self, sim):
#        
#        x, y = sim.xyGrid()
#        r = np.sqrt(x**2 + y**2)
#        phi = np.arctan2(y, x)
#        
#        V_rad = -1.0 # inject some radial velocity
#        V_x = V_rad * np.cos(phi)
#        V_y = V_rad * np.sin(phi)
#
#        sigma = 5.0 / sim.xhi
#        
#        zero_ind = int((sim.xhi + sim.xlo) / 2)
#        r_circle = int(5)   # input as integer
#        #circle_mask = (r[zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] <= r_circle)
#        ym, xm = np.ogrid[-r_circle:r_circle, -r_circle:r_circle]
#        circ_mask = (xm**2 + ym**2 < (1*r_circle)**2)
#
#      
#        #rho_0 = sim.V[sim.DN, zero_ind + r_circle + 1, zero_ind + r_circle + 1]
#        rho_0 = sim.V[sim.DN, zero_ind+6, zero_ind+6]
#        
#        #sim.V[sim.DN, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = rho_0 + 1.0 * np.exp(- 0.5 * r**2 / sigma**2)[zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] # gaussian density profile on top of whatever's there
#        #sim.V[sim.VX, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = V_x[zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] 
#        #sim.V[sim.VY, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = V_y[zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] 
#
#
#        sim.V[sim.DN, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = 1
#        sim.V[sim.VX, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = 0
#        sim.V[sim.VY, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = 0
#
#
#
#        sim.primitiveToConserved(sim.V, sim.U)
#        
#        return

    def setBoundaryConditions(self, sim):
        xhi = sim.xhi
        xlo = sim.xlo
        dx = sim.dx

        x, y = sim.xyGrid()
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        
        #ind_B = int(xlo + self.r_B / dx)
        
        sim.V[sim.DN, xhi:, 0] = self.rho_inf
        sim.V[sim.PR, xhi:, 0] = self.P_inf
        sim.V[sim.VX, xhi:, 0] = 0.0
        sim.V[sim.VY, xhi:, 0] = 0.0
        
        zero_ind = int((sim.xhi + sim.xlo) / 2)
        r_circle = int(3)   # input as integer
        #circle_mask = (r[zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] <= r_circle)
        ym, xm = np.ogrid[-r_circle:r_circle, -r_circle:r_circle]
        circ_mask = (xm**2 + ym**2 < (1*r_circle)**2)

      
        #rho_0 = sim.V[sim.DN, zero_ind + r_circle + 1, zero_ind + r_circle + 1]
        rho_0 = sim.V[sim.DN, zero_ind+6, zero_ind+6]
        
        #sim.V[sim.DN, zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] = rho_0 + 1.0 * np.exp(- 0.5 * r**2 / sigma**2)[zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] # gaussian density profile on top of whatever's there
        #sim.V[sim.VX, zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] = V_x[zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] 
        #sim.V[sim.VY, zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] = V_y[zero_ind-5:zero_ind+5, zero_ind-5:zero_ind+5] 
        
        sim.V[sim.DN, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = 1.0
        sim.V[sim.VX, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = 0
        sim.V[sim.VY, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = 0


        ## set sonic radius conditions
        #sim.V[sim.DN, ind_B, 0] = self.rho_B
        #sim.V[sim.PR, ind_B, 0] = self.P_B
        #sim.V[sim.VX, ind_B, 0] = self.cs_B 
#
        #sim.V[sim.VX, :xlo, 0] = min(sim.V[sim.VX, xlo, 0], -1e-10)
        
        sim.primitiveToConserved(sim.V[:, xhi:, :], sim.U[:, xhi:, :])
        
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


def runBondi2d():
    setup = SetupBondiParker()
    ulula_run.run(setup, tmax = 0.5, nx = 275, max_steps = 6000, q_plot = ['DN', 'PR', 'VR'], 
                  movie = True, movie_file_ext = 'gif', plot_2d=True, save_plots=True, plot_ghost_cells=True)
    
    return

def runBondi1d():
    setup = SetupBondiParker()
    ulula_run.run(setup, tmax = 0.5, nx = 275, max_steps=6000, q_plot = ['DN', 'PR', 'VR'],
                  save_plots=False, plot_1d=True, plot_2d=False, plot_step=50)
    return


#def runBondi():
#    setup = SetupBondiParker()
#    ulula_run.run(setup, tmax = 10, nx = 200, q_plot = ['DN', 'PR', 'VR'], 
#                  save_plots=False, plot_1d=True, plot_2d=False, plot_time=0.5)
#    return
