import numpy as np
from scipy.interpolate import interp1d
from .constants_list import *

class Model:
    #forward and backward shocks are "fwd" and "bwd"
    def __init__(self, filename: str, simtype: str, eps_B: float=0.1, eps_E: float=0.1, p_exp: float=3., gamma_min: float=1.):
        d = np.load(filename) #e.g. shock_data.npz
        self.times = d['times']
        self.vsh = d['vsh_of_t']
        self.rsh = d['rsh_of_t']

        self.int_rhofwd_sq_dr = d['int_rhofl_1_sq_dr'] #g^2/cm^5

        self.eps_B = eps_B
        self.eps_E = eps_E
        self.p = p_exp
        self.gamma_min = gamma_min

        self.rsh_t0 = self.rsh[0] # cm
        self.vsh_t0 = self.vsh[0] # cm/s
        self.tdyn_t0 = self.rsh[0] / self.vsh[0] # s

        if simtype == 'flare_flare':
            # this one will run twice for fwd and reverse shocks (uses deltav_1/2 and rhofl_1/2)
            self.vsh_fwd = self.deltav_1 =d['deltav1_of_t'] #cm/s
            self.rho_fwd = self.rhofl_1 = d['rho1_of_t'] #g/cm^3

            self.rho_bwd = self.rhofl_2 = d['rho2_of_t']
            self.vsh_bwd = self.deltav_2 = d['deltav2_of_t']

        elif simtype == 'flare_ism':
            # this one will run once for fwd shock only (uses rhoISM and vsh)
            self.vsh_fwd = self.vsh
            self.rho_fwd = self.rhoISM = d['rho_ism'] #g/cm^3

        else:
            raise ValueError("Invalid simulation type. Use 'flare_flare' or 'flare_ism'.")

        self.B_fwd = self.compute_B(self.rho_fwd,self.vsh_fwd)
        self.n0_fwd = self.compute_n0(self.rho_fwd,self.vsh_fwd) 

        if simtype == 'flare_flare':
            self.B_bwd = self.compute_B(self.rho_bwd,self.vsh_bwd)
            self.n0_bwd = self.compute_n0(self.rho_bwd,self.vsh_bwd) 
        #N_t0 = n0_t0 * rsh_t0**3, n0_t0 = n0[0]

    def compute_B(self,rho,v):
        if rho is None: return None
        return np.sqrt(8 * np.pi * self.eps_B * rho * v**2)
    def compute_n0(self,rho,v):
        if rho is None:  return None
        return self.eps_E*((self.p-2)*self.gamma_min**(self.p-2))*rho*v**2/(m_e*c**2)
    
    def generate_interp_funcs(self,simtype):
        self.rsh_func = interp1d(self.times, self.rsh, kind='linear')
        self.vsh_tot_func = interp1d(self.times, self.vsh, kind='linear')

        self.int_rhofwd_sq_dr_func = interp1d(self.times, self.int_rhofwd_sq_dr, kind='linear')

        self.vsh_fwd_func = interp1d(self.times, self.vsh_fwd, kind='linear')
        self.rho_fwd_func = interp1d(self.times, self.rho_fwd, kind='linear')
        self.B_fwd_func = interp1d(self.times, self.B_fwd, kind='linear')
        self.n0_fwd_func = interp1d(self.times, self.n0_fwd, kind='linear')

        if simtype == 'flare_flare':
            self.vsh_bwd_func = interp1d(self.times, self.vsh_bwd, kind='linear') 
            self.rho_bwd_func = interp1d(self.times, self.rho_bwd, kind='linear') 
            self.B_bwd_func = interp1d(self.times, self.B_bwd, kind='linear') 
            self.n0_bwd_func = interp1d(self.times, self.n0_bwd, kind='linear') 
        return

    def generate_ND_interp_funcs(self,simtype):
        self.rsh_ND_func = interp1d(self.times/self.tdyn_t0, self.rsh/self.rsh_t0, kind='linear')
        self.vsh_tot_ND_func = interp1d(self.times/self.tdyn_t0, self.vsh/self.vsh_t0, kind='linear')

        self.vsh_fwd_ND_func = interp1d(self.times/self.tdyn_t0, self.vsh_fwd/self.vsh_t0, kind='linear')
        self.n0_fwd_ND_func = interp1d(self.times/self.tdyn_t0, self.n0_fwd/self.n0_fwd[0], kind='linear')
        self.B_fwd_ND_func = interp1d(self.times/self.tdyn_t0, self.B_fwd, kind='linear')

        if simtype == 'flare_flare':
            self.vsh_bwd_ND_func = interp1d(self.times/self.tdyn_t0, self.vsh_bwd/self.vsh_t0, kind='linear') 
            self.n0_bwd_ND_func = interp1d(self.times/self.tdyn_t0, self.n0_bwd/self.n0_bwd[0], kind='linear') 
            self.B_bwd_ND_func = interp1d(self.times/self.tdyn_t0, self.B_bwd, kind='linear')
        return

        

