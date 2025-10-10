#  Synchrotron emission 
import numpy as np
from astropy.constants import e,h,m_p,m_e,c


def gamma_ph_func(gamma_ph_min=1.5e-18,gamma_ph_max=1.5,N_ph=256):
    d_ln_gamma_ph = (np.log(gamma_ph_max)-np.log(gamma_ph_min))/(N_ph-1)
    gamma_ph_vals = np.zeros(N_ph)
    gamma_ph_vals[0] = gamma_ph_min
    for i in np.arange(1,N_ph):
        gamma_ph_vals[i] = gamma_ph_vals[i-1] + (np.exp(d_ln_gamma_ph)-1)*gamma_ph_vals[i-1]

    return gamma_ph_vals

def nu_ph_func(gamma_ph_vals):
    #photon frequencies
   return gamma_ph_vals*m_e.cgs.value*c.cgs.value**2/h.cgs.value

def omega_c(B,gamma,sin_alpha=2./3.):#critical angular frequency
    return 3 * gamma**2 * e.gauss.value * B * sin_alpha/(2*m_e.cgs.value*c.cgs.value)
def sync_x_vals(B,gamma,nu_ph_vals):
    return nu_ph_vals/(omega_c(B,gamma)/2*np.pi)

def syn_func_fit(x):
    # /* analytical fitting of synchrotron function F(x) */
    # /* see http://arxiv.org/pdf/1301.6908.pdf */
    GAMMA13 = 2.67893 # Gamma(1/3) 
    F1 = np.pi * 2.0**(5.0/3.0) /np.sqrt(3.0)/GAMMA13 * x**(1.0/3.0)
    F2 = np.sqrt(np.pi/2.0)*np.exp(-x)*x**(1.0/2.0)

    a1_1 = -0.97947838884478688
    a1_2 = -0.83333239129525072
    a1_3 = 0.1554179602681624
    H_1 = a1_1 * x**1./1. + a1_2 * x**(1.0/2.0) + a1_3 * x**(1.0/3.0)
    delta_1 =np.exp(H_1)

    a2_1 = -0.0469247165562628882
    a2_2 = -0.70055018056462881
    a2_3 = 0.0103876297841949544
    H_2 = a2_1 * x**1./1. + a2_2 * x**(1.0/2.0) + a2_3 * x**(1.0/3.0)
    delta_2 = 1.0 - np.exp(H_2)

    return F1*delta_1+F2*delta_2

def emission_absorption_at_time(t,dNdgamma_vals,gamma_e_vals,delta_gamma_e,N_g,nu_ph_vals,
                                rsh_ND_func,rsh_t0,B_ND_func,f_omega=1):
    #t is in units of dynamical times
    # dNdgamma_vals is the electron spectrum at time t
    # gamma_e_vals is the electron Lorentz factor values, N_g is the number of gamma_e_vals
    # delta_gamma_e is the width of the electron Lorentz factor bins
    # rsh_ND_func,  B_ND_func are functions that return the shock radius and magnetic field at time t. 
    # rsh_t0 is the normalization constant for rsh_ND_func.

    sin_alpha=2./3. # angle averaged value
    L_nu_integral = 0
    alpha_SSA_integral = 0
    vol_shock_div_dr = f_omega*4*np.pi*(rsh_ND_func(t)*rsh_t0)**2
    for i in np.arange(0,N_g-1):
        gamma_e = gamma_e_vals[i]
        x_vals_tmp = sync_x_vals(B_ND_func(t),gamma_e,nu_ph_vals)
        Pnu = syn_func_fit(x_vals_tmp) * sin_alpha*np.sqrt(3)*e.gauss.value**3*B_ND_func(t)/(m_e.cgs.value*c.cgs.value**2)
        L_nu_integral += Pnu*dNdgamma_vals[i]*delta_gamma_e[i]
        d_dgamma_term = dNdgamma_vals[i+1]/gamma_e_vals[i+1]**2 - dNdgamma_vals[i]/gamma_e_vals[i]**2
        d_dgamma_term = d_dgamma_term/delta_gamma_e[i]
        alpha_SSA_integral += gamma_e**2 * Pnu * d_dgamma_term * delta_gamma_e[i]
    for i in [0,N_g-1]:
        gamma_e = gamma_e_vals[i]
        x_vals_tmp = sync_x_vals(B_ND_func(t),gamma_e,nu_ph_vals)
        Pnu = syn_func_fit(x_vals_tmp) * sin_alpha*np.sqrt(3)*e.gauss.value**3*B_ND_func(t)/(m_e.cgs.value*c.cgs.value**2)
        L_nu_integral += 0.5*Pnu*dNdgamma_vals[i]*delta_gamma_e[i]
        d_dgamma_term = dNdgamma_vals[i]/gamma_e_vals[i]**2
        d_dgamma_term = d_dgamma_term/delta_gamma_e[i]
        alpha_SSA_integral += 0.5*gamma_e**2 * Pnu * d_dgamma_term * delta_gamma_e[i]

    L_nu = L_nu_integral
    tau_SSA = -alpha_SSA_integral/(8*np.pi*nu_ph_vals**2*m_e.cgs.value)/vol_shock_div_dr #since tau = integral of alpha_SSA dr, using vol_shock_div_dr gives tau_SSA directly
    #assume helium CSM: n_e = rho/2/m_p, n_i = rho/4/m_p
    return L_nu, tau_SSA #each is an array over nu_ph

def Lnu_obs_spectrum(Lnu_syn,tau_ff,tau_ssa):
    return Lnu_syn*np.exp(-tau_ff)*(1-np.exp(-tau_ssa))/tau_ssa

def tau_ff(densitysq_integral,nu_ph,T_e_csm=1e4,X_h=0.7,Y_he=0.3): 
    # need to integrate alpha_ff for all densities outside r_shock(t). this is done in the densitysq_integral term.
    # currently assumes metal fraction is close to zero
    mu_e = (X_h+0.5*Y_he)**(-1)
    ni_Zsq_term = X_h*1**2/1 + Y_he*2**2/4 #this term is \sum(n_i Z^2) without the factor of rho/m_p. it equals 1 unless there is significant metals
    nu_scale = 10*1e9 #10 GHz
    tau_ff_arr = 8.4e-28 * (T_e_csm/1e4)**(-1.35) * (nu_ph/nu_scale)**(-2.1) * densitysq_integral*ni_Zsq_term/(mu_e*m_p.cgs.value**2) 
    return tau_ff_arr #array vs. nu_ph
