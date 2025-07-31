import numpy as np
from scipy.interpolate import interp1d
from constants_list import *

class Model:

    def __init__(self, filename: str, simtype: str, eps_B: float=0.1, eps_E: float=0.1, p_exp: float=3., gamma_min: float=1.):
        d = np.loadtxt(filename) 
        self.times = d[:,0] #s
        self.rsh = d[:,1] #cm
        self.vsh = d[:,2] #cm/s

        self.eps_B = eps_B
        self.eps_E = eps_E
        self.p = p_exp
        self.gamma_min = gamma_min

        self.rsh_t0 = self.rsh[0] # cm
        self.vsh_t0 = self.vsh[0] # cm/s
        self.tdyn_t0 = self.rsh[0] / self.vsh[0] # s

        if simtype == 'flare_flare':
            # this one will run twice for fwd and reverse shocks (uses deltav_1/2 and rhofl_1/2)
            self.deltav_2 = d[:,4] #cm/s
            self.rhofl_2 = d[:,6] #g/cm^3

            self.deltav_1 = d[:,3] #cm/s
            self.rhofl_1 = d[:,5] #g/cm^3

            self.int_rhofwd_sq_dr = d[:,7] #g^2/cm^5 

            self.rho_fwd = self.rhofl_1
            self.rho_bwd = self.rhofl_2

            self.v_fwd = self.deltav_1
            self.v_bwd = self.deltav_2
            #forward and backward shocks are "fwd" and "bwd"
        elif simtype == 'flare_ism':
            # this one will run once for fwd shock only (uses rhoISM and vsh)
            self.deltav_2 = d[:,3] #cm/s
            self.rhoISM = d[:,4] # g/cm^3

            self.rhofl_2 = d[:,5] # g/cm^3
            self.int_rhofwd_sq_dr = d[:,6] #g^2/cm^5 

            self.rho_fwd = self.rhoISM
            self.v_fwd = self.vsh
            
            self.rho_bwd = None
            self.v_bwd = None
        else:
            raise ValueError("Invalid simulation type. Use 'flare_flare' or 'flare_ism'.")

        self.B_fwd = self.compute_B(self.rho_fwd,self.v_fwd)
        self.n0_fwd = self.compute_n0(self.rho_fwd,self.v_fwd) 
        self.B_bwd = self.compute_B(self.rho_bwd,self.v_bwd)
        self.n0_bwd = self.compute_n0(self.rho_bwd,self.v_bwd) 
        #N_t0 = n0_t0 * rsh_t0**3, n0_t0 = n0[0]

    def compute_B(self,rho,v):
        if rho is None: return None
        return np.sqrt(8 * np.pi * self.eps_B * rho * v**2)
    def compute_n0(self,rho,v):
        if rho is None:  return None
        return self.eps_E*((self.p-2)*self.gamma_min**(self.p-2))*rho*v**2/(m_e*c**2)
    
    def generate_interp_funcs(self):
        self.rsh_func = interp1d(self.times, self.rsh, kind='linear')
        self.vsh_tot_func = interp1d(self.times, self.vsh, kind='linear')

        self.int_rhofwd_sq_dr_func = interp1d(self.times, self.int_rhofwd_sq_dr, kind='linear')

        self.vsh_fwd_func = interp1d(self.times, self.v_fwd, kind='linear')
        self.rho_fwd_func = interp1d(self.times, self.rho_fwd, kind='linear')
        self.B_fwd_func = interp1d(self.times, self.B_fwd, kind='linear')
        self.n0_fwd_func = interp1d(self.times, self.n0_fwd, kind='linear')

        self.vsh_bwd_func = interp1d(self.times, self.v_bwd, kind='linear') if self.v_bwd is not None else None
        self.rho_bwd_func = interp1d(self.times, self.rho_bwd, kind='linear') if self.rho_bwd is not None else None
        self.B_bwd_func = interp1d(self.times, self.B_bwd, kind='linear') if self.B_bwd is not None else None
        self.n0_bwd_func = interp1d(self.times, self.n0_bwd, kind='linear') if self.n0_bwd is not None else None
        return

    def generate_ND_interp_funcs(self):
        self.rsh_ND_func = interp1d(self.times, self.rsh/self.rsh_t0, kind='linear')
        self.vsh_tot_ND_func = interp1d(self.times, self.vsh/self.vsh_t0, kind='linear')

        self.vsh_fwd_ND_func = interp1d(self.times, self.v_fwd/self.vsh_t0, kind='linear')
        self.n0_fwd_ND_func = interp1d(self.times, self.n0_fwd/self.n0_fwd[0], kind='linear')

        self.vsh_bwd_ND_func = interp1d(self.times, self.v_bwd/self.vsh_t0, kind='linear') if self.v_bwd is not None else None
        self.n0_bwd_ND_func = interp1d(self.times, self.n0_bwd/self.n0_bwd[0], kind='linear') if self.n0_bwd is not None else None
        return

        

