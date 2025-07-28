import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import sigma_T,m_e,c
import argparse
from scipy.interpolate import interp1d
from matplotlib.colors import LogNorm

plt.style.use('/Users/samwu/codes/current_projects/RadioTDEFlares/constant_density_evol/plot_styles.mplstyle_new')
parser = argparse.ArgumentParser(description='''Analyze spectrum evolution ''')
parser.add_argument('--model_dir',type=str,default='',help='Path to input model directory')  #e.g. './2.898M_Porb10/'
parser.add_argument('--data_dir',type=str,default='./evolve_spectrum/',help='Path to spectrum output files')
parser.add_argument('--evolve_shock_dir',type=str,default='',help='Path to evolve shock directory') 
parser.add_argument('--dNdgamma_dir',type=str,default='',help='Path to evolve electron spectrum directory') 
parser.add_argument('--B_field_prof',type=str,default='None',help='Path to input magnetic field profile. should be npz file') #e.g. Bfield_vs_t.npz will be in data_dir
parser.add_argument('--eps_B',type=float,default=1e-1,help='value of epsilon_B (magnetic field efficiency factor)')
parser.add_argument('--eps_E',type=float,default=1e-1,help='value of epsilon_E (electron efficiency factor)')
parser.add_argument('--Te_csm',type=float,default=1e5,help='Temperature of CSM in K')
parser.add_argument('--f_omega',type=float,default=1,help='covering fraction of CSM')
parser.add_argument('--p_exp',type=float,default=3,help='exponent of electron power law (p>2)')
parser.add_argument('--epochs',type=int,nargs='+',default=[10,100,1000],help='Epochs to plot SED at, in days') #e.g. [0,1,2,3,4,5] for 6 epochs

# data_dir = './evolve_spectrum/E_1.00E+51_Mej_1.00E+00_rhoprof_2.898M_Porb10/'
# python spectrum_evolution.py 'evolve_spectrum/E_1.00E+51_Mej_1.00E+00_rhoprof_2.898M_Porb10/' \n
#'evolve_spectrum/evolve_spectrum_implicit/' --max_step=20000 --print_int=5000 --dt_sc=1e-3 --tf=1e6
args = parser.parse_args()
data_dir = args.data_dir
evolve_shock_dir = args.evolve_shock_dir
model_dir = args.model_dir
dNdgamma_over_dir = args.dNdgamma_dir
f_omega = args.f_omega
eps_E = args.eps_E
eps_B = args.eps_B
p = args.p_exp
T_e_csm = args.Te_csm #K
epoch_list = args.epochs

km_s = 1e5 #* cm/s
G = 6.67430e-8
Msun = 1.989e33
Lsun = 3.828e33
Rsun = 6.96e10

secinyear = 3.154e7
secinday = 86400

def B(eps_B,rho,vel): #inputs in CGS
    return np.sqrt(8*np.pi*eps_B*rho*vel**2) #in Gauss
#assuming gamma_max -> infty, normalization of dn/dgamma = n0 * gamma^-p is:
def n0(eps_E,gamma_min,rho,vel,p=p): #this is number density for the energy distribution dn/dgamma
    return eps_E*((p-2)*gamma_min**(p-2))*rho*vel**2/(m_e.cgs.value*c.cgs.value**2)
#usually assume gamma_min=1
#adopt p=3 also
dndgamma_func = lambda gamma,n0_const,p: n0_const*gamma**(-p)

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
gamma_ph_min = 1.5e-18 #1.5e-14
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
print("minimum nu_ph",nu_ph_vals[0],'maximum nu_ph',nu_ph_vals[-1])
#defined above:
# nu_ph_vals = gamma_ph_vals*m_e.cgs.value*c.cgs.value**2/h.cgs.value

# sync_x_vals = lambda B,gamma: nu_ph_vals/(omega_c(B,gamma)/2*np.pi)
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
        L_nu_integral += Pnu*dNdgamma_vals[i]*delta_gamma_e[i]
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

def tau_ff(densitysq_integral,nu_ph,T_e_csm=1e4,verbose=True): 
    #need to integrate alpha_ff for all densities outside r_shock(t). uints of t should be in t_dyn_t0 (or use other interp fn?)
    # assuming X=0.7, Y=0.3
    mu_e = (0.7+0.5*0.3)**(-1)
    nu_scale = 10*1e9 #10 GHz
    tau_ff_arr = 8.4e-28 * (T_e_csm/1e4)**(-1.35) * (nu_ph/nu_scale)**(-2.1) * densitysq_integral/(mu_e*m_p.cgs.value**2) 
    return tau_ff_arr #array vs. nu_ph


## Loading shock data and calculating emission for each flare
data_dir = args.data_dir
shock_data = np.load(data_dir + 'shock_data.npz')
times_orig = shock_data['times']
vsh_orig = shock_data['vsh_of_t']
rsh_of_t = shock_data['rsh_of_t']
rho_const = shock_data['rho_const']
rho2_of_t = shock_data['rho2_of_t']
# deltav1_of_t = shock_data['deltav1_of_t']
deltav2_of_t = shock_data['deltav2_of_t']
intrhofl1sqdr = shock_data['int_rhofl_1_sq_dr']

intrhofl1sqdr_vs_t = interp1d(times_orig,intrhofl1sqdr,kind='linear')

# flare_numbers = [0,1]
# for flare in flare_numbers:
#     if flare==0:
#         rho_of_t = rho1_of_t
#         vsh_of_t = deltav1_of_t
#         dNdgamma_dir = dNdgamma_over_dir + '/flare1/'
#     elif flare==1:
flare = 1
# rho_of_t = rho2_of_t
rho_of_t = rho_const
vsh_of_t = deltav2_of_t
dNdgamma_dir = dNdgamma_over_dir + '/flare2/'
print('flare',flare+1,'dNdgamma_dir',dNdgamma_dir)

if args.B_field_prof == 'None':
    B_of_t = B(eps_B=eps_B,rho=rho_of_t,vel=vsh_of_t)
    times = times_orig
    ind_start=0
elif args.B_field_prof != 'None':
    B_data = np.load(data_dir+args.B_field_prof)
    B_of_t = B_data['B_vals']
    print(B_of_t)
    B_of_t[0] = B_of_t[1]
    times = B_data['times']
    ind_start = int(len(times_orig)-len(times))
    print(ind_start)
    vsh_of_t = vsh_of_t[ind_start:]
    rsh_of_t = rsh_of_t[ind_start:]
    rho_of_t = rho_of_t[ind_start:]

n0_of_t = n0(eps_E=eps_E,gamma_min=1,rho=rho_of_t,vel=vsh_of_t,p=p) 
# B_of_t = B(eps_B=1e-2,rho=rho_of_t,vel=vsh_of_t)
#nondimensionalize to initial values of vsh, rsh, and initial dynamical time t_dyn ~ rsh/vsh
adjust_const = 1 #1e-12
vsh_t0 = vsh_orig[0]
rsh_t0 = rsh_of_t[0]
tdyn_t0 = rsh_t0/vsh_t0
#also scale N by initial value N_t0 = n0_t0 * rsh_t0**3
n0_t0 = n0_of_t[0]
N_t0 = n0_t0 * rsh_t0**3
#now new time is t' = t/tdyn_t0, new n0' = n0/n0_t0, new N' = N/(N_t0)
#create non-dimensional interpolating functions vs. time, time is now in units of tdyn_t0
vsh_tot_func_ND = interp1d(times/tdyn_t0,vsh_orig/vsh_t0,kind='linear')
vsh_func_ND = interp1d(times/tdyn_t0,vsh_of_t/vsh_t0,kind='linear')
rsh_func_ND = interp1d(times/tdyn_t0,rsh_of_t/rsh_t0,kind='linear')
n0_func_ND = interp1d(times/tdyn_t0,n0_of_t/n0_t0,kind='linear')
B_func = interp1d(times/tdyn_t0,B_of_t,kind='linear')

B_func = interp1d(times/tdyn_t0,B_of_t,kind='linear')
coeff_rad_ND = lambda t: -((sigma_T.cgs.value * B_func(t)**2)/(6 * np.pi * m_e.cgs.value * c.cgs.value))*tdyn_t0

print("initial dynamical time, shock radius, shock velocity, normalization constant", 
    f'{tdyn_t0:1.3E} s, {rsh_t0:1.3E} cm,{vsh_t0:1.3E} cm/s,{n0_t0:1.3E} cm^-3')
print('minimum time', times[0]/tdyn_t0,'maximum time', times[-1]/tdyn_t0, 'tdyn', times[-1]/secinyear, 'yr')

#time dependent coefficients for heating and cooling terms
c1 = lambda t: -vsh_tot_func_ND(t)/rsh_func_ND(t) #/adjust_const
c2 = lambda t: coeff_rad_ND(t) #/adjust_const
c3 = lambda t: 4 * np.pi * rsh_func_ND(t)**2 * vsh_func_ND(t) * n0_func_ND(t) #/adjust_const


print('c1',c1(times[0]/tdyn_t0),'c2',c2(times[0]/tdyn_t0), 'c3',c3(times[0]/tdyn_t0) )

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

##### Uncomment if want Peak Data ####
fig=plt.figure()

# colors = plt.cm.viridis(np.linspace(0,1,15))
norm = LogNorm(vmin=tvals[1000]*tdyn_t0/secinyear, vmax=tvals[-1]*tdyn_t0/secinyear)
count = 0
Lnu_pkvals = []
t_pkvals = []
nu_pkvals = []
# tau_ff_test_grid = np.zeros((15,len(nu_ph_vals)))
for t in tvals:
#times[np.where((rsh_func_ND(times/tdyn_t0)*rsh_t0>1e14) & (rsh_func_ND(times/tdyn_t0)*rsh_t0<1e15))]/tdyn_t0: #times/tdyn_t0:
    # count+=1
    # if count < 100 and count % 10 ==0:
    if count % 100 ==0:
        print(count,int(count/100), t)
        Lnu_test,tau_ssa_test = emission_absorption_at_time(tvals[count],dNdgamma_vals[count],gamma_e_vals,delta_gamma_e,f_omega=f_omega)
        tau_ff_test = tau_ff(intrhofl1sqdr_vs_t(tvals[count]*tdyn_t0),nu_ph_vals,T_e_csm=T_e_csm,verbose=False)
        Lnu_spectrum = Lnu_test*np.exp(-tau_ff_test)*(1-np.exp(-tau_ssa_test))/tau_ssa_test
        argmax,Lnumax = (np.argmax(Lnu_spectrum),np.amax(Lnu_spectrum))
        # print(nu_ph_vals[argmax]/1e9,Lnumax)
        Lnu_pkvals.append(Lnumax)
        nu_pkvals.append(nu_ph_vals[argmax])
        t_pkvals.append(t*tdyn_t0/secinyear)
        plt.plot(nu_ph_vals/1e9,Lnu_spectrum,color=plt.cm.viridis(norm(t*tdyn_t0/secinyear)))
        # plt.scatter(nu_ph_vals[argmax]/1e9,Lnumax,color=plt.cm.viridis(norm(t*tdyn_t0/secinyear)),s=40)
    count +=1
plt.plot(nu_ph_vals/1e9,1e30*(nu_ph_vals/1e9)**(-1.5),color='black',ls='--',label=r'$\nu^{-1.5}$')
plt.plot(nu_ph_vals/1e9,1e26*(nu_ph_vals/1e9)**(-1),color='grey',ls='--',label=r'$\nu^{-1}$')
sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
plt.colorbar(sm,label='Time (yr)')
plt.yscale('log')
plt.xscale('log')
# plt.axvspan(2,4,color='gray',alpha=0.5)
# plt.axhspan(1e26,1e29,color='gray',alpha=0.5)
plt.legend()
plt.xlim(1e-6,1e8)
plt.ylim(1e10,1e30)
plt.xlabel(r'$\nu$ (GHz)')
plt.ylabel(r'$L_{\nu}$ (erg/s/Hz)')
plt.savefig(f'Radio_curves_vs_frequency_{flare+1}.png',dpi=300,transparent=False,facecolor='white')
# # plt.title('Both SSA and FF')
# print(nu_ph_vals[np.where((nu_ph_vals > 1e10) & (nu_ph_vals < 10e10))])

np.savez(dNdgamma_dir+f'peak_data.npz',Lnu_pk=np.array(Lnu_pkvals),t_pk=np.array(t_pkvals),nu_pk=np.array(nu_pkvals))


#### To plot SED at specific epochs ####
fig=plt.figure()
norm = LogNorm(vmin=tvals[1000]*tdyn_t0/secinday, vmax=tvals[-1]*tdyn_t0/secinday)
print("Plotting SED at epochs (days)", epoch_list)
times_list = []
count_list = []

for epoch in epoch_list: #epochs are in days
    index_at_epoch=np.argmin(np.abs(tvals*tdyn_t0/secinday-epoch)/epoch )
    times_list.append(tvals[index_at_epoch]) #in units of tdyn_t0
    count_list.append(index_at_epoch)

SEDs_to_save_dict = {}
SEDs_to_save_dict['nu_ph_vals_GHz'] = nu_ph_vals/1e9 #in GHz

for t,count,epoch in zip(times_list,count_list,epoch_list):
#times[np.where((rsh_func_ND(times/tdyn_t0)*rsh_t0>1e14) & (rsh_func_ND(times/tdyn_t0)*rsh_t0<1e15))]/tdyn_t0: #times/tdyn_t0:
    # count+=1
    # if count < 100 and count % 10 ==0:

    print(count, t*tdyn_t0/secinday, 'days')
    Lnu_test,tau_ssa_test = emission_absorption_at_time(tvals[count],dNdgamma_vals[count],gamma_e_vals,delta_gamma_e,f_omega=f_omega)
    tau_ff_test = tau_ff(intrhofl1sqdr_vs_t(tvals[count]*tdyn_t0),nu_ph_vals,T_e_csm=T_e_csm,verbose=False)
    Lnu_spectrum = Lnu_test*np.exp(-tau_ff_test)*(1-np.exp(-tau_ssa_test))/tau_ssa_test
    SEDs_to_save_dict[f'epoch_{epoch}'] = Lnu_spectrum
    plt.plot(nu_ph_vals/1e9,Lnu_spectrum,color=plt.cm.viridis(norm(t*tdyn_t0/secinday)),label=f'{epoch:.1f} days'.format(epoch=epoch))
        # plt.scatter(nu_ph_vals[argmax]/1e9,Lnumax,color=plt.cm.viridis(norm(t*tdyn_t0/secinyear)),s=40)
    count +=1
plt.plot(nu_ph_vals/1e9,1e30*(nu_ph_vals/1e9)**(-1.5),color='black',ls='--',label=r'$\nu^{-1.5}$')
plt.plot(nu_ph_vals/1e9,1e26*(nu_ph_vals/1e9)**(-1),color='grey',ls='--',label=r'$\nu^{-1}$')
sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
plt.colorbar(sm,label='Time (day)')
plt.yscale('log')
plt.xscale('log')
# plt.axvspan(2,4,color='gray',alpha=0.5)
# plt.axhspan(1e26,1e29,color='gray',alpha=0.5)
plt.legend()
plt.xlim(1e-6,1e8)
plt.ylim(1e10,1e30)
plt.xlabel(r'$\nu$ (GHz)')
plt.ylabel(r'$L_{\nu}$ (erg/s/Hz)')
plt.savefig(f'SED_at_epochs_{flare+1}.png',dpi=300,transparent=False,facecolor='white')
np.savez(dNdgamma_dir+f'SED_data.npz',SED_vs_epoch=SEDs_to_save_dict)

