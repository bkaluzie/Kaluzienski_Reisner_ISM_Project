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

class SetupBondi(setup_base.Setup):
    
    """
    Spherical accretion setup
    
    In this setup, a gravitational potential is placed at the center of a spherically-symmetric gas cloud.
    The potential is at rest relative to distant regions of the cloud. Far from the potential, 
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
        
        return 'bondi'

    # ---------------------------------------------------------------------------------------------
    
   # need an update function to change the central mass (maybe?) 
   # start in uniform density, add material + velocity at the center to reinforce conditions
   # see how it moves out and how solution approaches parker
   # accretion -- set velocity at boundaries, let material flow in, update mass of object
   # can't add multiple velocities in single cell --- use gaussian (look at gresho vortex)
    
    def setInitialConditions(self, sim, nx):

        sim.setEquationOfState(eos_mode = 'ideal', gamma = self.gamma, pressure_floor = 0.001)
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
        
        sim.V[GP] = -self.g / (0.03 + r) # potential from TDE setup

        return
    
    # ---------------------------------------------------------------------------------------------
    
    def u(self, sim, r):
        cs_r = sim.soundSpeed(sim.V)
        
        return -(self.Mdot / (4.0 * np.pi * r**2)) * (cs_r / self.cs_inf)**(-2 / (self.gamma - 1.0))
    
    # ---------------------------------------------------------------------------------------------
    
    def setBoundaryConditions(self, sim):
        xhi = sim.xhi
        xlo = sim.xlo
        ylo = sim.ylo
        yhi = sim.yhi
        
        # set boundary conditions in all ghost cells
        sim.V[sim.DN, xhi+1:, yhi+1:] = self.rho_inf
        sim.V[sim.DN, xhi+1:, :ylo] = self.rho_inf
        sim.V[sim.DN, :xlo, yhi:] = self.rho_inf
        sim.V[sim.DN, :xlo, :ylo] = self.rho_inf
        
        sim.V[sim.PR, xhi+1:, yhi+1:] = self.P_inf
        sim.V[sim.PR, xhi+1:, :ylo] = self.P_inf
        sim.V[sim.PR, :xlo, yhi+1:] = self.P_inf
        sim.V[sim.PR, :xlo, :ylo] = self.P_inf
        
        
        sim.V[sim.VX, xhi+1:, yhi+1:] = 0
        sim.V[sim.VX, xhi+1:, :ylo] = 0
        sim.V[sim.VX, :xlo, yhi+1:] = 0
        sim.V[sim.VX, :xlo, :ylo] = 0
        
        sim.V[sim.VY, xhi+1:, yhi+1:] = 0
        sim.V[sim.VY, xhi+1:, :ylo] = 0
        sim.V[sim.VY, :xlo, yhi+1:] = 0
        sim.V[sim.VY, :xlo, :ylo] = 0
        
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
                vmax.append(10.0)
            elif q in ['VX', 'VY', 'VR', 'VA']:
                vmin.append(-5.0)
                vmax.append(0.0)
            
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

def runBondi2d():
    setup = SetupBondi()
    ulula_run.run(setup, tmax = 1.0, nx = 275, max_steps = 6000, q_plot = ['DN', 'PR', 'VR'], 
                  #plot_2d=True, save_plots=False, plot_time = 0.2)
                  movie = True, movie_file_ext='gif', movie_length=5.0, movie_fps=25, movie_dpi=500)
    return

def runBondi1d():
    setup = SetupBondi()
    ulula_run.run(setup, tmax = 1.0, nx = 275, max_steps=6000, q_plot = ['DN', 'PR', 'VR'],
                  save_plots=False, plot_1d=True, plot_2d=False, plot_time = 0.2)
    return