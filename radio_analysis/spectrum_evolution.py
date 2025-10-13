import numpy as np
import matplotlib.pyplot as plt
from functools import partial
import argparse
from copy import copy
import os
from scipy.interpolate import interp1d

# first, establish RHS function and the finite difference scheme for RHS
from .solvers import elec_time_evol #, euler, RK4, df_dx
# from calcs import evolve_ODE

script_dir = os.path.dirname(os.path.abspath(__file__))
style_path = os.path.join(script_dir, '../plot_styles.mplstyle_new')
plt.style.use(style_path)

from hydro_evol.model import Model
from hydro_evol.constants_list import *

def gamma_e_func(gamma_min,N_g,gamma_max=1e8):

    #set gamma_e values
    # gamma_max = 1e8 #actually calculated this at some point
    # gamma_min = m.gamma_min
    # N_g = m.N_g
    d_ln_gamma = (np.log(gamma_max) - np.log(gamma_min))/N_g
    # gamma_e_vals = np.arange(gamma_min,gamma_max,delta_gamma)
    gamma_e_vals = np.zeros(N_g)
    gamma_e_vals[0] = gamma_min
    for i in np.arange(1,N_g):
        gamma_e_vals[i] = gamma_e_vals[i-1] + (np.exp(d_ln_gamma)-1)*gamma_e_vals[i-1]
    delta_gamma_e = gamma_e_vals*(np.exp(d_ln_gamma)-1)
    return gamma_e_vals, delta_gamma_e, d_ln_gamma

#define some functions
dndgamma_func = lambda gamma,n0_const,p: n0_const*gamma**(-p)
def dot_gamma_e_ad(gamma_e,t,vsh_tot_func,rsh_func):
    return -(vsh_tot_func(t)/rsh_func(t)) * gamma_e
def dot_gamma_e_rad(gamma_e,t,B_func):
    return -((sigma_T * B_func(t)**2)/(6 * np.pi * m_e * c)) * gamma_e**2
def q_e(gamma_e,t,vsh_func,rsh_func,n0_func,p,f_omega=1):
    # return dndgamma_func(gamma_e,n0_func(t),p)
    return f_omega*4 * np.pi * rsh_func(t)**2 * vsh_func(t) * dndgamma_func(gamma_e,n0_func(t),p)
coeff_rad_ND_func = lambda t,B_func,tdyn_t0: -((sigma_T * B_func(t)**2)/(6 * np.pi * m_e * c))*tdyn_t0

def evolve_spectrum(simtype, data_dir, save_dir_in='', t0_in=0, dt_scale_in=1e-3, tf_in=1e2, max_step_in=1e2, print_int=1e2,plotting=True):
    """
    Run evolution of spectrum given hydrodynamical shock info.
    
    Parameters:
    simtype (str): Type of simulation ('flare_flare' or 'flare_ism').
    data_dir (str): Path to input properties of shock vs time.
    save_dir (str): Path to save output files of spectrum evolution.
    t0_in (float): Initial time in years.
    dt_scale_in (float): Scaling of dt (in units of t_dyn_0).
    tf_in (float): Final time in years.
    max_step_in (int): Maximum number of steps.
    print_int (int): Interval for plotting.
    """

    if save_dir_in=='':
        save_dir = data_dir
    else:
        save_dir = save_dir_in

    m = Model(data_dir+'shock_data.npz',simtype=simtype)
    m.generate_ND_interp_funcs(simtype)
    m.generate_interp_funcs(simtype)

    #set gamma_e values
    gamma_e_vals, delta_gamma_e, d_ln_gamma = gamma_e_func(m.gamma_min,m.N_g,gamma_max=1e8)

    if simtype=='flare_flare':
        flare_list=['fwd','bwd']
    elif simtype=='flare_ism':
        flare_list=['fwd']


    for flarenum,flare in enumerate(flare_list):
        if flare == 'fwd':
            coeff_rad_ND = partial(coeff_rad_ND_func,B_func=m.B_fwd_ND_func,tdyn_t0=m.tdyn_t0)
            vsh_ND_func = m.vsh_fwd_ND_func
            n0_func_ND = m.n0_fwd_ND_func
            n0_func = m.n0_fwd_func
            B_ND_func = m.B_fwd_ND_func
            B_of_t = m.B_fwd
            vsh_of_t = m.vsh_fwd
            
        elif flare == 'bwd':
            coeff_rad_ND = partial(coeff_rad_ND_func,B_func=m.B_bwd_ND_func,tdyn_t0=m.tdyn_t0)
            vsh_ND_func = m.vsh_bwd_ND_func
            n0_func_ND = m.n0_bwd_ND_func
            n0_func = m.n0_bwd_func
            B_ND_func = m.B_bwd_ND_func
            B_of_t = m.B_bwd
            vsh_of_t = m.vsh_bwd

        print(f"initial dynamical time: {m.tdyn_t0:1.3E} s, shock radius: {m.rsh_t0:1.3E} cm,")
        print(f"initial {flare} shock velocity: {m.vsh_t0:1.3E} cm/s, normalization constant for {flare} shock: {m.n0_fwd[0]:1.3E} cm^-3, initial B field for {flare} shock: {m.B_fwd[0]}") 
        print('minimum time (tdyn0)', m.times[0]/m.tdyn_t0,'max time',m.times[-1]/m.tdyn_t0,t0_in*secinyear/m.tdyn_t0)
        print('minimum time (yr)', m.times[0]/secinyear,'max time',m.times[-1]/secinyear,t0_in)

        #time dependent coefficients for heating and cooling terms
        c1 = lambda t: -m.vsh_tot_ND_func(t)/m.rsh_ND_func(t) #/adjust_const
        c2 = lambda t: coeff_rad_ND(t) #/adjust_const
        c3 = lambda t: 4 * np.pi * m.rsh_ND_func(t)**2 * vsh_ND_func(t) * n0_func_ND(t) #/adjust_const

        afunc = lambda x,t: c1(t) * x + c2(t) * x**2.
        qfunc = lambda x,t: c3(t) * x**(-3.)

        Afunc = lambda x,t: -2*c1(t)-c2(t)*x
        Qfunc = lambda x,t: c3(t)*np.ones(len(x))
        print('c1(0)',c1(m.times[0]/m.tdyn_t0),'c2(0)',c2(m.times[0]/m.tdyn_t0),'c3(0)', c3(m.times[0]/m.tdyn_t0) )

        if plotting:
            plt.figure()
            plt.plot(m.times/secinyear,B_ND_func(m.times/m.tdyn_t0))
            plt.plot(m.times/secinyear,B_of_t,ls=':')
            plt.xlabel('Time (yr)')
            plt.ylabel('B field')
            plt.xscale('log')
            plt.yscale('log')
            plt.show()
            plt.close()

            plt.figure()
            plt.plot(m.times/secinyear,vsh_ND_func(m.times/m.tdyn_t0))
            plt.plot(m.times/secinyear,vsh_of_t/m.vsh_t0,label=f'{flare} shock velocity',ls='--')
            plt.plot(m.times/secinyear,m.vsh/m.vsh_t0,label='shock velocity',ls=':',color='black')
            plt.xlabel('Time (yr)')
            plt.ylabel('vel')
            plt.xscale('log')
            plt.yscale('log')
            plt.legend()
            plt.show()
            plt.close()

        #_________set system parameters____________

        rundir = save_dir+f'/{flare}_flare'+str(flarenum+1)+'/'
        

        t0 = t0_in*secinyear/m.tdyn_t0 #secinyear converted to dynamical time
        dt0 = 1e-2*dt_scale_in/np.abs(c1(t0)) #should be dt_scale * initial dynamical time
        final_age = tf_in*secinyear/m.tdyn_t0  #final time in dynamical time units

        print('times/tdyn_t0',m.times/m.tdyn_t0,t0,final_age)

        max_step = int(max_step_in) #2000
        verbose = True
        print_interval = print_int #1e-2*max_step #print every 1% of steps
        save_interval = 1e-1*max_step #saves every timestep, but this is when to save in case need to restart
        print("print interval",print_interval,max_step)
        print(f'Creating save dir for {flare} flare{flarenum+1}',rundir)
        if not os.path.exists(rundir): os.mkdir(rundir)
        if not os.path.exists(rundir+'/plots/'): os.mkdir(rundir+'/plots/')

        
        print('------------initial values: -------------')

        count = 0
        min_dt = 1e-10
        dt_limit = 1e2 # choose an appropriate limit?
        stop = False
        dt_i = copy(dt0)
        t_i = copy(t0)
        #save values of times, dNdgamma values, etc.
        y_saved = np.array([],dtype=np.float64)
        t_saved = np.array([],dtype=np.float64)
        dt_saved = np.array([0],dtype=np.float64)
        
        print("initial shock radius, shock velocity, normalization constant", f'{m.rsh_func(t_i*m.tdyn_t0):1.3E},{m.vsh_tot_func(t_i*m.tdyn_t0):1.3E},{n0_func(t_i*m.tdyn_t0):1.3E}')

        #this version of y should be like dN/dgamma
        while (t_i < final_age) and (count < max_step) and (dt_i > min_dt) and (stop == False):
            redo_flag = False
            reduce_dt_flag = False
            if count == 0:
                t_i = copy(t0) #sec
                y_i = np.zeros_like(gamma_e_vals) #initial condition. y_i is dNdgamma here

                y_saved = np.append(y_saved,y_i)
                t_saved = np.append(t_saved,t_i)
            else:
                t_i = copy(t)
                y_i = copy(y_next) #array of length gamma_e_vals
            
            dt_i = copy(dt0)
            #just fix dt for now
            if count % print_interval == 0:
                # RHS_func = lambda y_n,t_n: ( Afunc(gamma_e_vals,t_n)*y_n + afunc(gamma_e_vals,t_n)*df_dx(y_n,delta_gamma) + Qfunc(gamma_e_vals,t_n) )

                print("----------------Beginning of Step-----------------")
                print('count',count,'time (tdyn)',f'{t_i:1.5E}')
                print('c1(t)',c1(t_i),'c2(t)',c2(t_i),'c3(t)', c3(t_i) )
                # print('Q',Qfunc(gamma_e_vals,t_i)[0],Qfunc(gamma_e_vals,t_i)[-1])


            try:
                # if t_i*m.tdyn_t0-t0_in*secinyear < 1e-2*secinyear:
                #     dt_i = 5e-2*dt_scale_in/np.abs(c1(t_i)) #should be inverse of dynamical time
                #     # print('reduced timestep',t_i,dt_i)
                # else:
                dt_i = dt_scale_in/np.abs(c1(t_i)) #should be inverse of dynamical time
            except:
                print("error in dt_i",dt_i,t_i,m.times/m.tdyn_t0)
            #### time evolution
                
            Energy = m_e * c**2
            delta_energy = Energy*gamma_e_vals*(np.exp(d_ln_gamma)-1)
            Pcool = -Energy*(c1(t_i)*gamma_e_vals + c2(t_i)*gamma_e_vals**2)
            q_e_inj = c3(t_i)*gamma_e_vals**(-m.p)
            y_next = elec_time_evol(dt_i,gamma_e_vals,delta_energy, y_i, q_e_inj, Pcool)
            t = t_i + dt_i
            
            
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
            print("dts (yr)",dt_saved*m.tdyn_t0/secinyear)
            print("times (yr)",t_saved*m.tdyn_t0/secinyear)
        else:
            inds = np.arange(0,count,int(np.floor(count/10)))
            print(inds)
            print("dts (yr)",dt_saved[inds]*m.tdyn_t0/secinyear,dt_saved[-1]*m.tdyn_t0/secinyear)
            print("times (yr)",t_saved[inds]*m.tdyn_t0/secinyear,t_saved[-1]*m.tdyn_t0/secinyear)
        # if t_i < final_age:
            np.savez(rundir+'/yvals.npz',y_saved)
            np.save(rundir+'/times.npy',t_saved)
            np.save(rundir+'/dts.npy',dt_saved)
            np.save(rundir+'/gamma_e_vals.npy',gamma_e_vals)
            
