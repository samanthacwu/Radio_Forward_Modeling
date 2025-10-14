import numpy as np
from scipy.interpolate import interp1d
from .constants_list import *

class Model:
    #forward and backward shocks are "fwd" and "bwd"
    def __init__(self, filename: str, simtype: str, eps_B: float=0.1, eps_E: float=0.1, p_exp: float=3., gamma_min: float=1., N_gamma: int=256, f_omega: float=1.,
                 X_H: float=0.7, X_He: float=0.3,
                calculate_SBO: bool=False, integrated_Bfield: bool=False):
        d = np.load(filename) #e.g. shock_data.npz
        self.directory_loc = filename.split('shock_data.npz')[0]
        self.times = d['times']
        self.vsh = d['vsh_of_t']
        self.rsh = d['rsh_of_t']

        self.int_rhofwd_sq_dr = d['int_rhofl_1_sq_dr'] #g^2/cm^5

        self.eps_B = eps_B
        self.eps_E = eps_E
        self.p = p_exp
        self.gamma_min = gamma_min
        self.N_g = N_gamma

        self.X_H = X_H
        self.X_He = X_He
        self.kappa = 0.2*(1+self.X_H) #cm^2/g, opacity

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
        elif simtype == 'SNejecta_CSM':
            self.vsh_fwd = self.vsh
            self.rho_fwd = self.rhoCSM = d['rho_csm'] #g/cm^3
            self.f_omega = f_omega
        else:
            raise ValueError("Invalid simulation type. Use 'flare_flare', 'flare_ism', or 'SNejecta_CSM'.")

        if calculate_SBO:
            i_start = self.shock_breakout_index(self.kappa)
            print('old t0 (yr)',self.times[0]/secinyear,'new t0 (yr)', self.times[i_start]/secinyear)
            self.rsh_t0 = self.rsh[i_start] # cm
            self.vsh_t0 = self.vsh[i_start] # cm/s
            self.tdyn_t0 = self.rsh[i_start] / self.vsh[i_start] # s

            self.times = self.times[i_start:]-self.times[i_start]
            self.vsh_fwd = self.vsh_fwd[i_start:]
            self.rsh = self.rsh[i_start:]
            self.rho_fwd = self.rho_fwd[i_start:]
            self.int_rhofwd_sq_dr = self.int_rhofwd_sq_dr[i_start:]

            print('new t0 (yr)', self.times[0]/secinyear)
        else:
            self.rsh_t0 = self.rsh[0] # cm
            self.vsh_t0 = self.vsh[0] # cm/s
            self.tdyn_t0 = self.rsh[0] / self.vsh[0] # s

        if integrated_Bfield:
            self.B_fwd =  self.Bfield_integrated()
        else:
            self.B_fwd = self.compute_B(self.rho_fwd,self.vsh_fwd)
        
        
        self.n0_fwd = self.compute_n0(self.rho_fwd,self.vsh_fwd) 

        if simtype == 'flare_flare':
            self.B_bwd = self.compute_B(self.rho_bwd,self.vsh_bwd)
            self.n0_bwd = self.compute_n0(self.rho_bwd,self.vsh_bwd) 
        #N_t0 = n0_t0 * rsh_t0**3, n0_t0 = n0[0]

    def compute_B(self,rho,v):
        # if want to use integrated B field, can compute that as well (see above)
        return np.sqrt(8 * np.pi * self.eps_B * rho * v**2)
    
    def compute_n0(self,rho,v):
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

    def shock_breakout_index(self,kappa=0.2):
        print('Calculating condition for onset of shock breakout')
        try:
            from scipy.integrate import simpson
        except:
            from scipy.integrate import simps as simpson

        #kappa value is  0.2 cm^2/g for H-poor fully ionized gas
        #calculate initial radius to start integrations for this model.
        rhosh_of_r = interp1d(self.rsh,self.rho_fwd,kind='linear')
        integrand = lambda r: kappa*rhosh_of_r(r)
        RHS = c/self.vsh[0]
        R0 = self.rsh[0]
        i_start = 0
        for i in np.arange(len(self.rsh)): #should be pretty close to the interior
            ans = simpson(integrand(self.rsh[i:]),self.rsh[i:])-RHS
            if ans < 0:
                R0 = self.rsh[i]
                i_start = i
                print('index',i, 'difference from c/v', ans)
                break
        return i_start

    def Bfield_integrated(self):
       
        import os
        from scipy.integrate import quad
        def Bsq_RHS(t): #inputs in CGS
            return self.f_omega*rhosh_func(t)*vsh_func(t)**3*rsh_func(t)**2 * rsh_func(t) #adiabatic expansion factor
        def coeff_of_t(time):
            vol = self.f_omega*(4./3.)*np.pi*rsh_func(time)**3
            return 8.*np.pi*self.eps_B*4*np.pi/vol/rsh_func(time)

        vsh_func = interp1d(self.times,self.vsh_fwd,kind='linear')
        rsh_func = interp1d(self.times,self.rsh,kind='linear')
        rhosh_func = interp1d(self.times,self.rho_fwd,kind='linear')    


        coeff_vals = np.zeros(len(self.times))
        dBsq_vals = np.zeros(len(self.times))
        B_vals = np.zeros(len(self.times))
        print('length of array',len(self.times))
        if not os.path.exists(self.directory_loc + 'Bfield_vs_t.npz'):
            print("integrating B field. Make sure shock_breakout condition is true.")
            for i, time in enumerate(self.times):
                if i==0:
                    dBsq_vals[i] = 0
                    continue
                if i % 500 == 0:
                    print('Bfield integration: ______i______',i)
                    print('time left',(self.times[-1]-time)/secinyear)
                    print('vsh',vsh_func(time)/1e5, 'km/s')
                    print('rsh',rsh_func(time)/AU_cm,'AU')
                    print('rho',rhosh_func(time), 'g/cm^3')
                    print('dBsq_vals',dBsq_vals[i-1])
                    # print('Bsq_scaled',Bsq_RHS(time)*coeff_of_t)
                # if i % 10 == 0: 
                Bsq_RHS_scaled = lambda t: Bsq_RHS(t)
                ans, err = quad(Bsq_RHS_scaled,self.times[i-1],time,limit=1000)
                    # if ans !=0:
                    #     if err/ans > 1e-1:
                    #         print('Warning: error in integration is large:',err/ans)
                dBsq_vals[i] = ans
                coeff_vals[i]=coeff_of_t(time)
                
            B_vals = np.sqrt(np.cumsum(dBsq_vals)*coeff_vals)
            np.savez(self.directory_loc + 'Bfield_vs_t.npz',times=self.times,B_vals=B_vals)
        else:
            print('Loading previously calculated B field values')
            Bfile = np.load(self.directory_loc + 'Bfield_vs_t.npz')
            B_vals = Bfile['B_vals']
        return B_vals