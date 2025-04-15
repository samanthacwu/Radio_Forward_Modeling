#### this is outdated ###

import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import sigma_T,m_e,c
import argparse
from scipy.interpolate import interp1d
from matplotlib.colors import LogNorm

plt.style.use('plot_styles.mplstyle_new')
parser = argparse.ArgumentParser(description='''Analyze spectrum evolution ''')
parser.add_argument('--model_dir',type=str,default='',help='Path to input model directory')  #e.g. './2.898M_Porb10/'
parser.add_argument('--data_dir',type=str,default='./evolve_spectrum/',help='Path to spectrum output files')
parser.add_argument('--evolve_shock_dir',type=str,default='',help='Path to evolve shock directory') 
parser.add_argument('--dNdgamma_dir',type=str,default='',help='Path to evolve electron spectrum directory') 
parser.add_argument('--B_field_prof',type=str,default='None',help='Path to input magnetic field profile. should be npz file') #e.g. Bfield_vs_t.npz will be in data_dir
parser.add_argument('--Te_csm',type=float,default=1e5,help='Temperature of CSM in K')
parser.add_argument('--f_omega',type=float,default=1,help='covering fraction of CSM')

# data_dir = './evolve_spectrum/E_1.00E+51_Mej_1.00E+00_rhoprof_2.898M_Porb10/'
# python spectrum_evolution.py 'evolve_spectrum/E_1.00E+51_Mej_1.00E+00_rhoprof_2.898M_Porb10/' \n
#'evolve_spectrum/evolve_spectrum_implicit/' --max_step=20000 --print_int=5000 --dt_sc=1e-3 --tf=1e6
args = parser.parse_args()
data_dir = args.data_dir
evolve_shock_dir = args.evolve_shock_dir
model_dir = args.model_dir
dNdgamma_dir = args.dNdgamma_dir
f_omega = args.f_omega

km_s = 1e5 #* cm/s
G = 6.67430e-8
Msun = 1.989e33
Lsun = 3.828e33
Rsun = 6.96e10

secinyear = 3.154e7
secinday = 86400

vels=np.load(evolve_shock_dir + 'shock_vels.npy')
rshocks=np.load(evolve_shock_dir + 'rshocks.npy')
Mencs=np.load(evolve_shock_dir + 'Mencs.npy')
times=np.load(evolve_shock_dir + 'times.npy')
dts=np.load(evolve_shock_dir + 'dts.npy')

rho_test = np.load(model_dir+'density_prof.npz')['rho']
r_test = np.load(model_dir+'density_prof.npz')['r']
density_interp = interp1d(r_test,rho_test,kind='linear')

def B(eps_B,rho,vel): #inputs in CGS
    return np.sqrt(8*np.pi*eps_B*rho*vel**2) #in Gauss
#assuming gamma_max -> infty, normalization of dn/dgamma = n0 * gamma^-p is:
def n0(eps_E,gamma_min,rho,vel): #this is number density for the energy distribution dn/dgamma
    return eps_E*gamma_min*rho*vel**2/(m_e.cgs.value*c.cgs.value**2)
#usually assume gamma_min=1
#adopt p=3 also
dndgamma_func = lambda gamma,n0_const,p: n0_const*gamma**(-p)

shock_data = np.load(data_dir + 'shock_data.npz')
times_orig = shock_data['times']
vsh_of_t = shock_data['vsh_of_t']
rsh_of_t = shock_data['rsh_of_t']
rho_of_t = shock_data['rho_of_t']


if args.B_field_prof == 'None':
    B_of_t = B(eps_B=1e-2,rho=rho_of_t,vel=vsh_of_t)
    times_t = times_orig
    ind_start=0
elif args.B_field_prof != 'None':
    B_data = np.load(data_dir+args.B_field_prof)
    B_of_t = B_data['B_vals']
    B_of_t[0] = B_of_t[1]
    times_t = B_data['times']
    ind_start = int(len(times_orig)-len(times_t))
    
    vsh_of_t = vsh_of_t[ind_start:]
    rsh_of_t = rsh_of_t[ind_start:]
    rho_of_t = rho_of_t[ind_start:]
    print(ind_start,len(times_t),len(vsh_of_t))

n0_of_t = n0(eps_E=1e-1,gamma_min=1,rho=rho_of_t,vel=vsh_of_t)
# B_of_t = B(eps_B=1e-2,rho=rho_of_t,vel=vsh_of_t)
#nondimensionalize to initial values of vsh, rsh, and initial dynamical time t_dyn ~ rsh/vsh
adjust_const = 1 #1e-12
vsh_t0 = vsh_of_t[0]
rsh_t0 = rsh_of_t[0]
tdyn_t0 = rsh_t0/vsh_t0
#also scale N by initial value N_t0 = n0_t0 * rsh_t0**3
n0_t0 = n0_of_t[0]
N_t0 = n0_t0 * rsh_t0**3

vsh_func_ND = interp1d(times_t/tdyn_t0,vsh_of_t/vsh_t0,kind='linear')
rsh_func_ND = interp1d(times_t/tdyn_t0,rsh_of_t/rsh_t0,kind='linear')
n0_func_ND = interp1d(times_t/tdyn_t0,n0_of_t/n0_t0,kind='linear')

B_func = interp1d(times_t/tdyn_t0,B_of_t,kind='linear')
coeff_rad_ND = lambda t: -((sigma_T.cgs.value * B_func(t)**2)/(6 * np.pi * m_e.cgs.value * c.cgs.value))*tdyn_t0

print("initial dynamical time, shock radius, shock velocity, normalization constant", 
      f'{tdyn_t0:1.3E} s, {rsh_t0:1.3E} cm,{vsh_t0:1.3E} cm/s,{n0_t0:1.3E} cm^-3')
print('minimum time', times_t[0]/tdyn_t0,'maximum time', times_t[-1]/tdyn_t0, 'tdyn', times_t[-1]/secinyear, 'yr')

#time dependent coefficients for heating and cooling terms
c1 = lambda t: -vsh_func_ND(t)/rsh_func_ND(t)/adjust_const
c2 = lambda t: coeff_rad_ND(t)/adjust_const
c3 = lambda t: 4 * np.pi * rsh_func_ND(t)**2 * vsh_func_ND(t) * n0_func_ND(t)/adjust_const


print('c1',c1(times_t[0]),'c2',c2(times_t[0]), 'c3',c3(times_t[0]) )

# delta_gamma = 100
delta_gamma = 1e-3 # not using even spacing
gamma_max = 1e8 #actually calculate this
gamma_min = 1
N_g = 256
d_ln_gamma = (np.log(gamma_max) - np.log(gamma_min))/N_g
# gamma_e_vals = np.arange(gamma_min,gamma_max,delta_gamma)
gamma_e_vals = np.zeros(N_g)
gamma_e_vals[0] = gamma_min
for i in np.arange(1,N_g):
    gamma_e_vals[i] = gamma_e_vals[i-1] + (np.exp(d_ln_gamma)-1)*gamma_e_vals[i-1]
delta_gamma_e = gamma_e_vals*(np.exp(d_ln_gamma)-1)

dt_scale_sim = 1e-3

dts = np.load(dNdgamma_dir+'/dts.npy')
tvals=np.load(dNdgamma_dir+'/times.npy')
yvals=np.load(dNdgamma_dir+'/yvals.npz')['arr_0'].reshape((len(tvals),len(gamma_e_vals)))
#convert yvals to the right units
dNdgamma_vals = yvals * n0_t0 * rsh_t0**3
print('normalization',n0_t0 * rsh_t0**3)
tdyn_of_t_sim = dts/dt_scale_sim

tvals_sec = tvals*tdyn_t0

#  Synchrotron emission 
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

from astropy.constants import e,h,m_p

#critical angular frequency
sin_alpha = 2./3. # angle averaged s.t. sin_alpha -> 2/3
omega_c = lambda B,gamma: 3 * gamma**2 * e.gauss.value * B * sin_alpha/(2*m_e.cgs.value*c.cgs.value)

#photon energies
gamma_ph_min = 1.5e-14
gamma_ph_max = 1.5
N_ph = 256
d_ln_gamma_ph = (np.log(gamma_ph_max)-np.log(gamma_ph_min))/(N_ph-1)
gamma_ph_vals = np.zeros(N_ph)
gamma_ph_vals[0] = gamma_ph_min
for i in np.arange(1,N_ph):
    gamma_ph_vals[i] = gamma_ph_vals[i-1] + (np.exp(d_ln_gamma_ph)-1)*gamma_ph_vals[i-1]

#photon frequencies
nu_ph_vals = gamma_ph_vals*m_e.cgs.value*c.cgs.value**2/h.cgs.value

sync_x_vals = lambda B,gamma: nu_ph_vals/(omega_c(B,gamma)/2*np.pi)

#defined above:
# nu_ph_vals = gamma_ph_vals*m_e.cgs.value*c.cgs.value**2/h.cgs.value

# sync_x_vals = lambda B,gamma: nu_ph_vals/(omega_c(B,gamma)/2*np.pi)
from scipy.integrate import quad
rhoshock_func = interp1d(times_t/tdyn_t0,rho_of_t,kind='linear')
#dNdgamma_vals should be yvals returned to correct units
def emission_absorption_at_time(t,dNdgamma_vals,gamma_e_vals,delta_gamma_e,f_omega=1):
    #time should be in units of tdyn_t0
    sin_alpha=2./3. # angle averaged value
    L_nu_integral = 0
    alpha_SSA_integral = 0
    vol_shock_div_dr = f_omega*4*np.pi*(rsh_func_ND(t)*rsh_t0)**2
    for i in np.arange(0,N_g-1):
        gamma_e = gamma_e_vals[i]
        x_vals_tmp = sync_x_vals(B_func(t),gamma_e)
        Pnu = syn_func_fit(x_vals_tmp) * sin_alpha*np.sqrt(3)*e.gauss.value**3*B_func(t)/(m_e.cgs.value*c.cgs.value**2)
        # if i % 50 == 0:
        #     plt.plot(nu_ph_vals,Pnu)
        #     plt.xscale('log')
            # plt.yscale('log')
        L_nu_integral += Pnu*dNdgamma_vals[i]*delta_gamma_e[i]
        # if i % 50 == 0:
        #     plt.plot(nu_ph_vals,L_nu_integral)
        #     plt.xscale('log')
        #     plt.yscale('log')
        #     plt.ylim(0.1,)
        d_dgamma_term = dNdgamma_vals[i+1]/gamma_e_vals[i+1]**2 - dNdgamma_vals[i]/gamma_e_vals[i]**2
        d_dgamma_term = d_dgamma_term/delta_gamma_e[i]
        alpha_SSA_integral += gamma_e**2 * Pnu * d_dgamma_term * delta_gamma_e[i]
    for i in [0,N_g-1]:
        gamma_e = gamma_e_vals[i]
        x_vals_tmp = sync_x_vals(B_func(t),gamma_e)
        Pnu = syn_func_fit(x_vals_tmp) * sin_alpha*np.sqrt(3)*e.gauss.value**3*B_func(t)/(m_e.cgs.value*c.cgs.value**2)
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
T_e_csm = args.Te_csm #K
def tau_ff(t,density_profile,radii,r_out,nu_ph,T_e_csm=1e5,verbose=True): 
    #need to integrate alpha_ff for all densities outside r_shock(t). uints of t should be in t_dyn_t0 (or use other interp fn?)
    Z=2
    # T_e_csm = 1e5 #K
    nu_scale = 10*1e9 #10 GHz
    alpha_ff_arr = 8.4e-28 * Z**2 * density_profile**2/(8*m_p.cgs.value**2) * (T_e_csm/1e4)**(-1.35) * (nu_ph/nu_scale)**(-2.1)
    alpha_ff_func = interp1d(radii,alpha_ff_arr,kind='linear')
    if verbose:
        plt.plot(radii,alpha_ff_arr)
        plt.scatter(rsh_func_ND(t)*rsh_t0,alpha_ff_func(rsh_func_ND(t)*rsh_t0))
        plt.grid()
        plt.yscale('log')
        plt.xscale('log')
        print('alpha_ff(r(t))', alpha_ff_func(rsh_func_ND(t)*rsh_t0))
    ans,err = quad(alpha_ff_func,rsh_func_ND(t)*rsh_t0,r_out,limit=10000)
    if verbose:
        print('error/ans', err/ans)
    print_once = False
    if err/ans > 1e-2:
        if verbose:
            print('warning: error in tau_ff integration is large', err/ans,'nu',nu_ph)

    return ans

plt.figure()
colors = plt.cm.viridis(np.linspace(0,1,15))
count = 0
tau_ff_test_grid = np.zeros((15,len(nu_ph_vals)))
for t in tvals:
#times[np.where((rsh_func_ND(times/tdyn_t0)*rsh_t0>1e14) & (rsh_func_ND(times/tdyn_t0)*rsh_t0<1e15))]/tdyn_t0: #times/tdyn_t0:
    # count+=1
    # if count < 100 and count % 10 ==0:
    if count > 7000 and count % 600 ==0:
        print(count, t)
        tau_ff_test = np.zeros_like(nu_ph_vals)
        for i,nu in enumerate(nu_ph_vals):
            tau_ff_test[i] = tau_ff(t,density_interp(rshocks),rshocks,rshocks[-1],nu,T_e_csm=T_e_csm,verbose=False)
        tau_ff_test_grid[int((count-7000)/600)] = tau_ff_test
        plt.plot(nu_ph_vals,tau_ff_test,color=colors[int((count-7000)/600)])
    count +=1
plt.yscale('log')
plt.xscale('log')
# plt.ylim(1e20,1e30)
# plt.xlim(1e-2,1000)
# plt.axvspan(2,4,color='gray',alpha=0.5)
# plt.axhspan(1e26,1e29,color='gray',alpha=0.5)

fig=plt.figure()

# colors = plt.cm.viridis(np.linspace(0,1,15))
norm = LogNorm(vmin=tvals[7200]*tdyn_t0/secinyear, vmax=tvals[-1]*tdyn_t0/secinyear)
count = 0
tau_ff_test_grid = np.zeros((15,len(nu_ph_vals)))
for t in tvals:
#times[np.where((rsh_func_ND(times/tdyn_t0)*rsh_t0>1e14) & (rsh_func_ND(times/tdyn_t0)*rsh_t0<1e15))]/tdyn_t0: #times/tdyn_t0:
    # count+=1
    # if count < 100 and count % 10 ==0:
    if count > 7000 and count % 600 ==0:
        print(count, t)
        Lnu_test,tau_ssa_test = emission_absorption_at_time(tvals[count],dNdgamma_vals[count],gamma_e_vals,delta_gamma_e,f_omega=f_omega)
        tau_ff_test = tau_ff_test_grid[int((count-7000)/600)]
        
        plt.plot(nu_ph_vals/1e9,Lnu_test*np.exp(-tau_ff_test)*(1-np.exp(-tau_ssa_test))/tau_ssa_test,color=plt.cm.viridis(norm(t*tdyn_t0/secinyear)))
    count +=1
sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
plt.colorbar(sm,label='Time (yr)')
plt.yscale('log')
plt.xscale('log')
plt.axvspan(2,4,color='gray',alpha=0.5)
plt.axhspan(1e26,1e29,color='gray',alpha=0.5)
plt.xlim(0.1,1e6)
plt.ylim(1e20,)
plt.xlabel(r'$\nu$ (GHz)')
plt.ylabel(r'$L_{\nu}$ (erg/s/Hz)')
plt.savefig('Radio_curves_vs_frequency.png',dpi=300,transparent=False,facecolor='white')
# plt.title('Both SSA and FF')

plt.figure()
colors = plt.cm.viridis(np.linspace(0,1,15))
count = 0
tau_ff_3ghz= np.zeros_like(tvals)
Lnu_3ghz = np.zeros_like(tvals)
tau_ssa_3ghz = np.zeros_like(tvals)
nu_in = 3.129e10 # 3 GHz
times_out = np.zeros_like(tvals)
j=0
for i,t in enumerate(tvals):
#times[np.where((rsh_func_ND(times/tdyn_t0)*rsh_t0>1e14) & (rsh_func_ND(times/tdyn_t0)*rsh_t0<1e15))]/tdyn_t0: #times/tdyn_t0:
    # count+=1
    # if count < 100 and count % 10 ==0:
    # if count > 7000 and count % 600 ==0:
    #     print(count, t)
    if i % 100 ==0:
        # print(i)
        # for i,nu in enumerate(nu_ph_vals):
        tau_ff_3ghz[j] = tau_ff(t,density_interp(rshocks),rshocks,rshocks[-1],nu_in,T_e_csm=T_e_csm,verbose=False)

        index_3ghz=np.where(np.abs(nu_ph_vals-3e10)/3e10 < 0.05)[0]
        Lnu_test,tau_ssa_test = emission_absorption_at_time(tvals[i],dNdgamma_vals[i],gamma_e_vals,delta_gamma_e,f_omega=f_omega)
        Lnu_3ghz[j] = Lnu_test[index_3ghz]
        tau_ssa_3ghz[j] = tau_ssa_test[index_3ghz]
        times_out[j] = t
        j+=1
    # count +=1
# print(i)
        
tau_ff_3ghz= np.trim_zeros(tau_ff_3ghz)
Lnu_3ghz = np.trim_zeros(Lnu_3ghz)
tau_ssa_3ghz = np.trim_zeros(tau_ssa_3ghz)
times_out = np.trim_zeros(times_out)
# print(tau_ssa_3ghz)
Lnu_abs_3ghz = Lnu_3ghz*np.exp(-tau_ff_3ghz)*(1-np.exp(-tau_ssa_3ghz))/tau_ssa_3ghz
plt.plot(times_out*tdyn_t0/secinyear,Lnu_3ghz*np.exp(-tau_ff_3ghz)*(1-np.exp(-tau_ssa_3ghz))/tau_ssa_3ghz)
#,color=colors[int((count-7000)/600)])
plt.yscale('log')
plt.xlim(0,100)
plt.ylim(1e20,1e29)
# plt.xscale('log')
plt.xlabel('Time (yr)')
plt.ylabel(r'$L_{\nu}$ (erg/s/Hz)')
plt.title('Emission at 3 GHz')
# plt.ylim(1e20,1e30)
# plt.xlim(1e-2,1000)
# plt.axvspan(2,4,color='gray',alpha=0.5)
plt.axhspan(1e26,1e29,color='gray',alpha=0.5)
plt.axvspan(5,20,color='gray',alpha=0.5)
plt.savefig('Radio_curve_3GHz.png',dpi=300,transparent=False,facecolor='white')

#times_sec should be times_yr oops
np.savez(dNdgamma_dir+f'Lnu_3ghz_sparse_Te_{T_e_csm:1.1E}.npz',times_out=times_out,Lnu_3ghz=Lnu_3ghz,tau_ff_3ghz=tau_ff_3ghz,tau_ssa_3ghz=tau_ssa_3ghz,
         times_sec=times_out*tdyn_t0/secinyear,Lnu_abs_3ghz=Lnu_abs_3ghz)
