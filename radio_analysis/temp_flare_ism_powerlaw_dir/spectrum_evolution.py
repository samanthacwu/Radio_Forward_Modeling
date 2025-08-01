import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import sigma_T,m_e,c
from functools import partial
import argparse
from copy import copy
import os
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser(description=''' Run evolution of shock radius given input density profile, energy, ejecta mass. ''')
parser.add_argument('data_dir',type=str,default='',help='Path to input properties of shock vs time. should contain npz files with standard names')
parser.add_argument('save_dir',type=str,default='',help='Path to save output files')
parser.add_argument('--B_field_prof',type=str,default='None',help='Path to input magnetic field profile. should be npz file')
#e.g. Bfield_vs_t.npz will be in data_dir
parser.add_argument('--t0',type=float,default=0,help='initial time (years)')
parser.add_argument('--dt_sc',type=float,default=1e-5,help='scaling of dt (t_dyn_0)')
parser.add_argument('--tf',type=float,default=1e2,help='final time (years)')
parser.add_argument('--max_step',type=float,default=1e2,help='max no. of steps')
parser.add_argument('--print_int',type=float,default=1e2,help='interval for plotting')
parser.add_argument('--eps_B',type=float,default=1e-2,help='value of epsilon_B (B field efficiency)')
parser.add_argument('--eps_E',type=float,default=1e-1,help='value of epsilon_E (electron efficiency factor)')
parser.add_argument('--f_omega',type=float,default=1,help='covering fraction of CSM')
parser.add_argument('--p_exp',type=float,default=3,help='exponent of electron power law (p>2)')

#will load rsh_of_t, vsh_of_t, and times from npz files for each model/evolution results
args = parser.parse_args()
data_dir = args.data_dir
save_dir = args.save_dir

km_s = 1e5 #* cm/s
G = 6.67430e-8
Msun = 1.989e33
Lsun = 3.828e33
Rsun = 6.96e10

secinyear = 3.154e7
secinday = 86400

t0_in = args.t0
# dt_in = args.dt_in
dt_scale_in = args.dt_sc
tf_in = args.tf
max_step_in = args.max_step
eps_B = args.eps_B
eps_E = args.eps_E
f_omega = args.f_omega
p = args.p_exp

#eps_B ~1e-2, eps_E ~1e-2-1e-1 from radio modeling
def B(eps_B,rho,vel): #inputs in CGS
    return np.sqrt(8*np.pi*eps_B*rho*vel**2) #in Gauss
#assuming gamma_max -> infty, normalization of dn/dgamma = n0 * gamma^-p is:
def n0(eps_E,gamma_min,rho,vel,p=p): #this is number density for the energy distribution dn/dgamma
    return eps_E*((p-2)*gamma_min**(p-2))*rho*vel**2/(m_e.cgs.value*c.cgs.value**2)
#usually assume gamma_min=1
#adopt p=3 also
dndgamma_func = lambda gamma,n0_const,p: n0_const*gamma**(-p)

data_dir = args.data_dir
shock_data = np.load(data_dir + 'shock_data.npz')
times_orig = shock_data['times']
vsh_orig = shock_data['vsh_of_t']
rsh_of_t = shock_data['rsh_of_t']
rho_const = shock_data['rho_const']
rho2_of_t = shock_data['rho2_of_t']
# deltav1_of_t = shock_data['deltav1_of_t']
deltav2_of_t = shock_data['deltav2_of_t']

# flare_numbers = [0,1]
# for flare in flare_numbers:
flare = 1 #just do one flare

# rho_of_t = rho2_of_t
rho_of_t= rho_const # use forward shock
vsh_of_t = vsh_orig
# print('flare',flare+1)

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

#create interpolating functions vs. time
vsh_tot_func = interp1d(times,vsh_orig,kind='linear')
vsh_func = interp1d(times,vsh_of_t,kind='linear')
rsh_func = interp1d(times,rsh_of_t,kind='linear')
B_func = interp1d(times,B_of_t,kind='linear')
n0_func = interp1d(times,n0_of_t,kind='linear')

print(len(times),len(B_of_t),len(vsh_of_t))
def dot_gamma_e_ad(gamma_e,t,vsh_tot_func,rsh_func):
    return -(vsh_tot_func(t)/rsh_func(t)) * gamma_e
def dot_gamma_e_rad(gamma_e,t,B_func):
    return -((sigma_T.cgs.value * B_func(t)**2)/(6 * np.pi * m_e.cgs.value * c.cgs.value)) * gamma_e**2
def q_e(gamma_e,t,vsh_func,rsh_func,n0_func,p=p,f_omega=1):
    # return dndgamma_func(gamma_e,n0_func(t),p)
    return f_omega*4 * np.pi * rsh_func(t)**2 * vsh_func(t) * dndgamma_func(gamma_e,n0_func(t),p)

#set gamma_e values
#set spacing (evenly spaced) of gamma_e values
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

# adjust_const = 1e30 #1e17


dot_gamma_e_ad_of_t = partial(dot_gamma_e_ad,vsh_tot_func=vsh_tot_func,rsh_func=rsh_func)
dot_gamma_e_rad_of_t = partial(dot_gamma_e_rad,B_func=B_func)
q_e_of_t = partial(q_e,vsh_func=vsh_func,rsh_func=rsh_func,n0_func=n0_func,p=3,f_omega=f_omega)

#nondimensionalize to initial values of vsh, rsh, and initial dynamical time t_dyn ~ rsh/vsh
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
coeff_rad_ND = lambda t: -((sigma_T.cgs.value * B_func(t)**2)/(6 * np.pi * m_e.cgs.value * c.cgs.value))*tdyn_t0

print("initial dynamical time, shock radius, shock velocity, normalization constant", 
    f'{tdyn_t0:1.3E} s, {rsh_t0:1.3E} cm,{vsh_t0:1.3E} cm/s,{n0_t0:1.3E} cm^-3')
print('minimum time', times[0]/tdyn_t0,'initial B field',B_func(times[0]/tdyn_t0),B_of_t[0],'max time',times[-1]/tdyn_t0)

#time dependent coefficients for heating and cooling terms
c1 = lambda t: -vsh_tot_func_ND(t)/rsh_func_ND(t) #/adjust_const
c2 = lambda t: coeff_rad_ND(t) #/adjust_const
c3 = lambda t: 4 * np.pi * rsh_func_ND(t)**2 * vsh_func_ND(t) * n0_func_ND(t) #/adjust_const

afunc = lambda x,t: c1(t) * x + c2(t) * x**2.
qfunc = lambda x,t: c3(t) * x**(-3.)

#define new functions for new equation formulation
#y = x^3 * u
Afunc = lambda x,t: -2*c1(t)-c2(t)*x
Qfunc = lambda x,t: c3(t)*np.ones(len(x))
print('c1(0)',c1(times[0]/tdyn_t0),'c2(0)',c2(times[0]/tdyn_t0),'c3(0)', c3(times[0]/tdyn_t0) )
# plt.figure()
# # plt.plot(gamma_e_vals,-afunc(gamma_e_vals,t=times[0]))
# # plt.plot(gamma_e_vals,Afunc(gamma_e_vals,t=times[0]))
# # plt.plot(gamma_e_vals,Qfunc(gamma_e_vals,t=times[0]))
# plt.plot(gamma_e_vals,-afunc(gamma_e_vals,t=times[0]))
# plt.plot(gamma_e_vals,qfunc(gamma_e_vals,t=times[0]))
# plt.yscale('log')
# plt.xscale('log')
# plt.xlabel('Gamma')
# plt.ylabel('A, Q')
# plt.show()
# plt.close()

plt.figure()
plt.plot(times/secinyear,B_func(times/tdyn_t0))
plt.plot(times_orig[ind_start:]/secinyear,B(eps_B=eps_B,rho=rho_of_t,vel=vsh_of_t),ls=':')
plt.xlabel('Time (yr)')
plt.ylabel('B field')
plt.xscale('log')
plt.yscale('log')
plt.show()
plt.close()

plt.figure()
plt.plot(times/secinyear,vsh_func_ND(times/tdyn_t0))
plt.plot(times_orig[ind_start:]/secinyear,vsh_of_t/vsh_t0,label=f'delta_v{flare+1}',ls=':')
plt.plot(times_orig[ind_start:]/secinyear,vsh_orig/vsh_t0,label='shock velocity',ls=':',color='black')
plt.xlabel('Time (yr)')
plt.ylabel('vel')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
plt.close()

# first, establish RHS function and the finite difference scheme for RHS
from solvers import elec_time_evol, euler, RK4, df_dx
from calcs import evolve_ODE
plt.style.use('plot_styles.mplstyle_new')

#quickly test the FD scheme - it works :)
# a_prime=df_dx2(a_func(gamma_e_vals,t=times[0]), delta_gamma)

# da_dgamma = -vsh_func(times[0])/rsh_func(times[0]) - 2*((sigma_T.cgs.value * B_func(times[0])**2)/(6 * np.pi * m_e.cgs.value * c.cgs.value)) * gamma_e_vals

# plt.plot(gamma_e_vals,a_prime)
# plt.plot(gamma_e_vals,da_dgamma,ls=':')
# plt.show()

#_________set system parameters____________

# print('rho_prof',args.rho_prof.split("/")[1])
rundir = args.save_dir+'/flare'+str(flare+1)+'/'


#*secinyear # 1e-7*secinyear #1e4*secinyear
t0 = t0_in*secinyear/tdyn_t0 #secinyear converted to dynamical time
dt0 =dt_scale_in/np.abs(c1(t0)) #should be dt_scale * initial dynamical time
final_age = tf_in*secinyear/tdyn_t0 #*secinyear #1.6e7*secinyear #1.6e9*secinyear #ages_array[-1]

print('times/tdyn_t0',times/tdyn_t0,t0,final_age)

max_step = int(max_step_in) #2000
verbose = True
print_interval = args.print_int #1e-2*max_step #print every 1% of steps
save_interval = 1e-1*max_step #saves every timestep, but this is when to save in case need to restart
# plot_interval = int(plotint_in) #10000
# save_interval = int(saveint_in) #1e6 
print("print interval",print_interval,max_step)
# time.sleep(2)
print(f'Creating save dir for flare{flare+1}',rundir)
if not os.path.exists(rundir): os.mkdir(rundir)
if not os.path.exists(rundir+'/plots/'): os.mkdir(rundir+'/plots/')

print('------------initial values: -------------')

count = 0
# dt_scale = 0.0003
min_dt = 1e-10 #*secinyear
dt_limit = 1e2 #*secinyear # choose an appropriate limit?
stop = False
dt_i = copy(dt0)
t_i = copy(t0)
#save values of times, dNdgamma values, etc.
y_saved = np.array([],dtype=np.float64)
t_saved = np.array([],dtype=np.float64)
dt_saved = np.array([],dtype=np.float64)

print("initial shock radius, shock velocity, normalization constant", f'{rsh_func(t_i*tdyn_t0):1.3E},{vsh_func(t_i*tdyn_t0):1.3E},{n0_func(t_i*tdyn_t0):1.3E}')

#this version of y should be like dN/dgamma
while (t_i < final_age) and (count < max_step) and (dt_i > min_dt) and (stop == False):
    redo_flag = False
    reduce_dt_flag = False
    if count == 0:
        t_i = copy(t0) #sec
        y_i = np.zeros_like(gamma_e_vals) #initial condition. y_i is dNdgamma here
    else:
        t_i = copy(t)
        y_i = copy(y_next) #array of length gamma_e_vals
    
    # if t_i < 1e-1*secinyear:
    #     # dt_i = 1e-7*secinyear
    #     dt_i = 1e-6*secinyear
    # else:
    dt_i = copy(dt0)
    #just fix dt for now
    if count % print_interval == 0:
        # RHS_func = lambda y_n,t_n: ( Afunc(gamma_e_vals,t_n)*y_n + afunc(gamma_e_vals,t_n)*df_dx(y_n,delta_gamma) + Qfunc(gamma_e_vals,t_n) )

        print("----------------Beginning of Step-----------------")
        print('count',count,'time (tdyn)',f'{t_i:1.5E}')
        print('c1(t)',c1(t_i),'c2(t)',c2(t_i),'c3(t)', c3(t_i) )
        # print('Q',Qfunc(gamma_e_vals,t_i)[0],Qfunc(gamma_e_vals,t_i)[-1])


    try:
        dt_i = dt_scale_in/np.abs(c1(t_i)) #should be inverse of dynamical time
    except:
        print("error in dt_i",dt_i,t_i,times/tdyn_t0)
    #change to the smaller of dynamical time, injection time
    # print('dt',dt_i)
    #include their dt adaptation
    #### time evolution
    ###### RHS: du/dt = (d/dx)f(u,t) + q_e(u,t)
    # RHS_func = lambda y_n,t_n: ( Afunc(gamma_e_vals,t_n)*y_n + afunc(gamma_e_vals,t_n)*df_dx(y_n,delta_gamma) + Qfunc(gamma_e_vals,t_n) )
    Energy = m_e.cgs.value * c.cgs.value**2
    delta_energy = Energy*gamma_e_vals*(np.exp(d_ln_gamma)-1)
    Pcool = -Energy*(c1(t_i)*gamma_e_vals + c2(t_i)*gamma_e_vals**2)
    q_e_inj = c3(t_i)*gamma_e_vals**(-p)
    y_next = elec_time_evol(dt_i,gamma_e_vals,delta_energy, y_i, q_e_inj, Pcool)
    t = t_i + dt_i
    # evolve_ODE(y_i,t_i,dt_i, RHS_func,solver=RK4) #RHS_func and RK4 are functions
    
    t_old = copy(t_i) #for next time, in case need to check
    y_old = copy(y_i)

    if count % print_interval == 0:
        print(f'dt (t_dyn): {dt_i:1.2E}','dt limit (t_dyn)',f'{dt_limit:1.5E}')
        print(f'y_next:{y_next[0]:1.2E},{y_next[-1]:1.2E}')
        print('y_i',y_i[0],y_i[-3:])

    #save values-from end of step
    y_saved = np.append(y_saved,y_next)
    t_saved = np.append(t_saved,t)
    dt_saved = np.append(dt_saved,dt_i)

    if count % save_interval == 0:
        np.savez(rundir+'/yvals.npz',y_saved)
        np.save(rundir+'/times.npy',t_saved)
        np.save(rundir+'/dts.npy',dt_saved)
    
    if (dt_i <= min_dt):
        print("stopping for min dt limit")
    count+=1

#save final values

print("number of steps",count)
np.set_printoptions(formatter={'float':'{:1.8E}'.format})
if count < 10:
    print("dts (yr)",dt_saved*tdyn_t0/secinyear)
    print("times (yr)",t_saved*tdyn_t0/secinyear)
else:
    inds = np.arange(0,count,int(np.floor(count/10)))
    print(inds)
    print("dts (yr)",dt_saved[inds]*tdyn_t0/secinyear,dt_saved[-1]*tdyn_t0/secinyear)
    print("times (yr)",t_saved[inds]*tdyn_t0/secinyear,t_saved[-1]*tdyn_t0/secinyear)
# if t_i < final_age:
    np.savez(rundir+'/yvals.npz',y_saved)
    np.save(rundir+'/times.npy',t_saved)
    np.save(rundir+'/dts.npy',dt_saved)
