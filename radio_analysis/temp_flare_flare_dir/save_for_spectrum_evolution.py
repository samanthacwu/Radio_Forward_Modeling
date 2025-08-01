# plot some quantities to analyze the evolution.py and Bfield_vs_t.py results, then save into combined npz file in the directory where we want to analyze spectrum evolution 
#loading in text file with these headers: 
# Mfl=0.01Msun, vmin=0.05c, vmax=0.408248c, power-law-index=0.5
# time [s], rsh [cm], vsh [cm/s], Deltav_1 [cm/s], Deltav_2 [cm/s], rhofl_1 [g/cm3], rhofl_2 [g/cm3], int_rhosq_1dr [g^2/cm^5]

import numpy as np
import matplotlib.pyplot as plt
import mesa_reader as mr
import os
import sys
import argparse
from scipy.interpolate import interp1d


parser = argparse.ArgumentParser(description='''Save and analyze shock properties. ''')
parser.add_argument('--evolve_shock_dir',type=str,default='',help='Path to evolve shock directory')  #e.g. shell_propagation/
parser.add_argument('--filename',type=str,default='',help='Name of file to load') #e.g. shell_evolution_deltat_1yr.txt
parser.add_argument('--data_dir',type=str,default='./evolve_spectrum/',help='Path to save output files')
parser.add_argument('--Bfield',type=bool,default=False,help='Whether to plot Bfield')
parser.add_argument('--eps_B',type=float,default=1e-2,help='value of epsilon_B (B field efficiency)')

args = parser.parse_args()
evolve_shock_dir = args.evolve_shock_dir
data_dir = args.data_dir+'/'
datafile_name = '/' + args.filename 

print(evolve_shock_dir.split('/')[-2])
print('loading',evolve_shock_dir+datafile_name)

if evolve_shock_dir=='':
    print('Please provide evolve shock directory')
    sys.exit()

G = 6.67430e-8
Msun = 1.989e33
Lsun = 3.828e33
Rsun = 6.96e10

secinyear = 3.154e7
secinday = 86400
plt.style.use('plot_styles.mplstyle_new')

shell_data=np.loadtxt(evolve_shock_dir+datafile_name)
times = shell_data[:,0] #s
rsh = shell_data[:,1] #cm
vsh = shell_data[:,2] #cm/s
deltav_1 = shell_data[:,3] #cm/s
deltav_2 = shell_data[:,4] #cm/s
rhofl_1 = shell_data[:,5] #g/cm^3
rhofl_2 = shell_data[:,6] #g/cm^3
int_rhofl_1_sq_dr = shell_data[:,7] #g^2/cm^5 


#python evolution.py './2.898M_Porb10/density_prof.npz' --dt_in=1e-2 --E=1e51 --max_step=50000 --print_int=1000 --tf=250
plt.plot(rsh/215/Rsun,vsh/1e5,color='black',ls=':',label=r'$v_{\rm sh}$')
plt.plot(rsh/215/Rsun, deltav_1/1e5,color='crimson',label=r'$\Delta v_{1}$')
plt.plot(rsh/215/Rsun, deltav_2/1e5,color='dodgerblue',label=r'$\Delta v_{2}$')
# plt.scatter(rsh/215/Rsun,vsh/1e5,s=10,marker='x',c='tab:olive')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('radius (AU)')
plt.ylabel('shell velocity (km/s)')
plt.legend()
plt.show()
plt.close()
print('Initial shell velocity',vsh[0]/1e5,'initial time (s)',times[0:2])

# plt.title(r'SN: $E=10^{51}$ erg, $M_{\rm ej}=1\, M_{\odot}$, Progenitor: $M_{\rm He}=2.898\,M_{\odot}$, $P_{\rm orb}=10$d')
# plt.plot(r_test/Rsun/215,rho_test,label='full density profile')
plt.plot(rsh/215/Rsun,rhofl_1,label='density traversed by shock 1')
plt.plot(rsh/215/Rsun,rhofl_2,label='density traversed by shock 2')
plt.plot(rsh/Rsun/215,6e11*rsh**(-2),label=r'$r^{-2}$ wind',color='black',ls=':')
# plt.axvline(5e2)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('radius (AU)')
plt.ylabel(r'density (g/cm$^3$)')
plt.ylim(1e-30,1e-12)
plt.legend()
plt.show()
plt.close()

# plt.title(r'SN: $E=10^{51}$ erg, $M_{\rm ej}=1\, M_{\odot}$, Progenitor: $M_{\rm He}=2.898\,M_{\odot}$, $P_{\rm orb}=10$d')
plt.plot(times/secinyear,rhofl_1,label='density traversed by shock 1')
plt.plot(times/secinyear,rhofl_2,label='density traversed by shock 2')
# plt.plot(times/secinyear,6e11*r_test**(-2),label=r'$r^{-2}$ wind')
# plt.axvline(5e2)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Time (yr)')
plt.ylabel(r'density (g/cm$^3$)')
plt.ylim(1e-30,1e-12)
plt.legend()
plt.show()
plt.close()

plt.plot(times/secinyear, vsh/1e5,color='black',ls=':',label=r'$v_{\rm sh}$')
plt.plot(times/secinyear, (vsh-deltav_1)/1e5,color='crimson',label=r'$v_{1}$')
plt.plot(times/secinyear, (deltav_2+vsh)/1e5,color='dodgerblue',label=r'$v_{2}$')
plt.xlabel('Time (yr)')
plt.ylabel('Shell Velocity (km/s)')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

# edit this to be done for each flare

if not args.Bfield:
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    else:
        print('Data dir exists')
    save_dir = data_dir+datafile_name.split('.txt')[0]
    print('Creating save dir',save_dir)
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    else:
        print('Directory exists')

    np.savez(save_dir+'/shock_data.npz',times=times,rsh_of_t=rsh,vsh_of_t=vsh,deltav1_of_t=deltav_1,rho1_of_t=rhofl_1,deltav2_of_t=deltav_2,rho2_of_t=rhofl_2,int_rhofl_1_sq_dr=int_rhofl_1_sq_dr) 

# if args.Bfield:
#     save_dir = data_dir+evolve_shock_dir.split('/')[-2]
#     from astropy.constants import m_e,c,sigma_T
#     km_s = 1e5 #* cm/s
#     eps_B=args.eps_B
#     #eps_B ~1e-2, eps_E ~1e-2-1e-1 from radio modeling
#     def B(eps_B,rho,vel): #inputs in CGS
#         return np.sqrt(8*np.pi*eps_B*rho*vel**2) #in Gauss
#     #assuming gamma_max -> infty, normalization of dn/dgamma = n0 * gamma^-p is:
#     def n0(eps_E,gamma_min,rho,vel): #this is number density for the energy distribution dn/dgamma
#         return eps_E*gamma_min*rho*vel**2/(m_e.cgs.value*c.cgs.value**2)
#     #usually assume gamma_min=1
#     #adopt p=3 also
#     dndgamma_func = lambda gamma,n0_const,p: n0_const*gamma**(-p)
#     B_func = interp1d(times,B(eps_B,density_interp(rsh),vsh),kind='cubic')
#     Bfield_test=np.load(save_dir+'/Bfield_vs_t.npz')

#     ind_start = int(len(times)-len(Bfield_test['times']))

#     plt.plot(Bfield_test['times']/secinyear,Bfield_test['B_vals'],label='Integrated B field')
#     plt.plot(times[ind_start:]/secinyear,B_func(times[ind_start:]),label='Original def')
#     # plt.plot(times,np.sqrt(Bfield_int*6*eps_B/(rsh**3)))
#     plt.yscale('log')
#     plt.xscale('log')
#     plt.ylabel('B field (G)')
#     plt.xlabel('Time (yr)')
#     plt.legend()
#     plt.show()
#     plt.close()
