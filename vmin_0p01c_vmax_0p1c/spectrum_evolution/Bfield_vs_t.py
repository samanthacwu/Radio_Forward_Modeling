import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import sigma_T,m_e,c
from functools import partial
import argparse
from copy import copy
import os
from scipy.interpolate import interp1d
from scipy.integrate import simpson,quad
parser = argparse.ArgumentParser(description=''' Run evolution of shock radius given input density profile, energy, ejecta mass. ''')
parser.add_argument('data_dir',type=str,default='',help='Path to input properties of shock vs time. should contain npz files with standard names')
# parser.add_argument('save_dir',type=str,default='',help='Path to save output files')
parser.add_argument('--eps_B',type=float,default=1e-2,help='value of epsilon_B (B field efficiency)')
parser.add_argument('--t0',type=float,default=0,help='initial time (years)')
parser.add_argument('--R0',type=float,default=1.1e-1,help='initial radius (AU)')
parser.add_argument('--f_omega',type=float,default=1,help='covering fraction of CSM')
parser.add_argument('--dt_in',type=float,default=1e-5,help='dt (years)')
parser.add_argument('--tf',type=float,default=1e2,help='final time (years)')
parser.add_argument('--max_step',type=float,default=1e2,help='max no. of steps')
parser.add_argument('--print_int',type=float,default=10,help='interval for plotting')
#will load rsh_of_t, vsh_of_t, and times from npz files for each model/evolution results
args = parser.parse_args()
data_dir = args.data_dir
# save_dir = args.save_dir

km_s = 1e5 #* cm/s
G = 6.67430e-8
Msun = 1.989e33
Lsun = 3.828e33
Rsun = 6.96e10
AU_cm = 1.496e13
AU_Rsun = 215
secinyear = 3.154e7
secinday = 86400

t0_in = args.t0
R0_in = args.R0
dt_in = args.dt_in
f_omega = args.f_omega
# dt_scale_in = args.dt_sc
tf_in = args.tf
max_step_in = args.max_step

shock_data = np.load(data_dir + 'shock_data.npz')
times = shock_data['times']
dts = shock_data['dts']
vsh_of_t = shock_data['vsh_of_t']
rsh_of_t = shock_data['rsh_of_t']
rho_of_t = shock_data['rho_of_t']
rhosh_of_r = interp1d(rsh_of_t,rho_of_t,kind='linear')

#calculate initial radius to start integration.
kappa=0.2 # cm^2/g for H-poor fully ionized gas
integrand = lambda r: kappa*rhosh_of_r(r)
RHS = c.cgs.value/vsh_of_t[0]
R0 = R0_in
i_start = 0
for i in np.arange(len(rsh_of_t)): #should be pretty close to the interior
    ans = simpson(integrand(rsh_of_t[i:]),rsh_of_t[i:])-RHS
    if ans < 0:
        R0 = rsh_of_t[i]
        i_start = i
        print(i, ans)
        break

# subtract t0 from times s.t. start integration at t=0
t0 = times[i_start]
times_int = times[i_start:] - t0
dts_int = dts[i_start:]
#create interpolating functions vs. time
vsh_func = interp1d(times_int,vsh_of_t[i_start:],kind='linear')
rsh_func = interp1d(times_int,rsh_of_t[i_start:],kind='linear')
rhosh_func = interp1d(times_int,rho_of_t[i_start:],kind='linear')

#
#eps_B ~1e-2, eps_E ~1e-2-1e-1 from radio modeling
eps_B=args.eps_B #1e-2 default


def Bsq_RHS(t,time,f_omega=1): #inputs in CGS
    return f_omega*rhosh_func(t)*vsh_func(t)**3*rsh_func(t)**2 * rsh_func(t) #adiabatic expansion factor
def coeff_of_t(time,f_omega=1):
    vol = f_omega*(4./3.)*np.pi*rsh_func(time)**3
    return 8.*np.pi*eps_B*4*np.pi/vol/rsh_func(time)
# B_vals = np.zeros(int(len(times)/10)+1)
# t_saved = np.zeros(int(len(times)/10)+1)
coeff_vals = np.zeros(len(times_int))
dBsq_vals = np.zeros(len(times_int))
B_vals = np.zeros(len(times_int))
print('length of array',len(times_int))
# t_saved = np.zeros(len(times_int))
if not os.path.exists(data_dir + 'Bfield_vs_t.npz'):
    print('Calculating B field vs time for',data_dir)
    for i, (dt,time) in enumerate(zip(dts_int,times_int)):
        if i==0:
            dBsq_vals[i] = 0
            continue
        if i % 500 == 0:
            print('______i______',i)
            print('time left',times_int[-1]-time)
            print('vsh',vsh_func(time))
            print('rsh',rsh_func(time))
            print('rho',rhosh_func(time))
            print('dBsq_vals',dBsq_vals[i-1])
            # print('Bsq_scaled',Bsq_RHS(time)*coeff_of_t)
        # if i % 10 == 0: 
        Bsq_RHS_scaled = lambda t: Bsq_RHS(t,time,f_omega=f_omega)
        ans, err = quad(Bsq_RHS_scaled,times_int[i-1],time,limit=1000)
            # if ans !=0:
            #     if err/ans > 1e-1:
            #         print('Warning: error in integration is large:',err/ans)
        dBsq_vals[i] = ans
        coeff_vals[i]=coeff_of_t(time,f_omega=f_omega)
            # t_saved[int(i/10)] = time
    B_vals = np.sqrt(np.cumsum(dBsq_vals)*coeff_vals)
    np.savez(data_dir + 'Bfield_vs_t.npz',times=times_int,B_vals=B_vals)
elif os.path.exists(data_dir + 'Bfield_vs_t.npz'):
    print('Loading B field vs time')
    B_data = np.load(data_dir + 'Bfield_vs_t.npz')
    B_vals = B_data['B_vals']
    times_int = B_data['times']
    eps_B_ratio = eps_B/1e-2 # saved Bfields should all be 1e-2 for now
    if eps_B_ratio != 1:
        print('Scaling B field values by ratio of',eps_B_ratio, 'w/sqrt', np.sqrt(eps_B_ratio))
    B_vals = B_vals*np.sqrt(eps_B_ratio) #this is to save time instead of integrating 
    np.savez(data_dir + f'Bfield_vs_t_{eps_B:1.0E}.npz',times=times_int,B_vals=B_vals)
def B(eps_B,rho,vel): #inputs in CGS
    return np.sqrt(8*np.pi*eps_B*rho*vel**2) #in Gauss
# print(B_vals)
plt.plot(times_int,B_vals)
plt.plot(times,B(eps_B,rho_of_t,vsh_of_t))
plt.yscale('log')
plt.xscale('log')
# plt.show()
plt.savefig(data_dir + 'Bfield_vs_t.png')
print('Completed B field vs time')



