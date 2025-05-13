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

class SetupParker(setup_base.Setup):
    
    """
    Spherical outflow setup
    
    In this setup, a star is placed at the center of a spherically-symmetric gas cloud.
    The star is at rest relative to distant regions of the cloud. Far from the star, 
    the cloud has a uniform density rho_inf and a uniform pressure P_inf. 
    This setup demonstrates
    
    * Ideal equation of state
    * Sub- and super-sonic flows
    * Polar coordinate solver.
    
    Parameters
    ----------
    name: string
        'bondi' runs the setup for spherical accretion; 'parker' runs the setup for spherical outflow.
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
        
        return 'parker'

    # ---------------------------------------------------------------------------------------------
    
    def setInitialConditions(self, sim, nx):

        sim.setEquationOfState(eos_mode = 'ideal', gamma = self.gamma, pressure_floor = 1E-10)
        sim.setGravityMode(gravity_mode = 'fixed_pot')
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
    
    def updateFunction(self, sim):
        x, y = sim.xyGrid()
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        
        sigma = 5.0 / sim.xhi
        
        V_rad = 2.0 #* np.exp(- 0.5 * r**2 / sigma**2) # inject some radial velocity
        
        V_x = V_rad * np.cos(phi)
        V_y = V_rad * np.sin(phi)
        
        #V_x = V_rad * r / (x + y**2/x)
        #V_y = V_x * y/x
        
        
        zero_ind = int((sim.xhi + sim.xlo) / 2)
        r_circle = int(5)   # input as integer
        ym, xm = np.ogrid[-r_circle:r_circle, -r_circle:r_circle]
        circ_mask = (xm**2 + ym**2 < (1*r_circle)**2)
      
        rho_0 = sim.V[sim.DN, zero_ind + r_circle, zero_ind]
        
        sim.V[sim.DN, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = rho_0 + 2.0 * np.exp(- 0.5 * r**2 / sigma**2)[zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] # gaussian density profile on top of whatever's there
        sim.V[sim.VX, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = V_x[zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] 
        sim.V[sim.VY, zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] = V_y[zero_ind-r_circle:zero_ind+r_circle, zero_ind-r_circle:zero_ind+r_circle][circ_mask] 

        sim.primitiveToConserved(sim.V, sim.U)
       
        return
    
    # ---------------------------------------------------------------------------------------------
    
    def setBoundaryConditions(self, sim):
        
        # implement pressure and density caps to prevent breaking
        pressure_cap = 50.0
        density_cap = 50.0
        
        if np.any(sim.V[sim.PR] > pressure_cap):
            sim.V[sim.PR][sim.V[sim.PR] > pressure_cap] = pressure_cap
            
        if np.any(sim.V[sim.DN] > density_cap):
            sim.V[sim.DN][sim.V[sim.DN] > density_cap] = density_cap
        
        sim.primitiveToConserved(sim.V, sim.U)
        
        return
    
    # ---------------------------------------------------------------------------------------------
    
    def plotLimits(self, q_plot, plot_geometry):
        
        vmin = []
        vmax = []

        for q in q_plot:
            if q in ['DN', 'PR']:
                vmin.append(0.0)
                vmax.append(15.0)
            elif q in ['VX', 'VY', 'VR', 'VA']:
                vmin.append(0.0)
                vmax.append(2.0)
            
            else:
                vmin.append(None)
                vmax.append(None)
        
        return vmin, vmax, None
    
    # ---------------------------------------------------------------------------------------------
    
    def plotColorMaps(self, q_plot):
        
        cmap = []
        
        for q in q_plot:
            cmap.append('viridis')
            
        return cmap
    
# ---------------------------------------------------------------------------------------------

def runParker2d():
    setup = SetupParker()
    ulula_run.run(setup, tmax = 1.0, nx = 275, max_steps = 6000, q_plot = ['DN', 'PR', 'VR'], 
                  #plot_2d=True, save_plots=False, plot_time = 0.2)
                  movie = True, movie_file_ext='gif', movie_length=5.0, movie_fps=25, movie_dpi=500)
    return

def runParker1d():
    setup = SetupParker()
    ulula_run.run(setup, tmax = 1.0, nx = 275, max_steps=6000, q_plot = ['DN', 'PR', 'VR'],
                  save_plots=False, plot_1d=True, plot_2d=False, plot_time = 0.2)
    return
