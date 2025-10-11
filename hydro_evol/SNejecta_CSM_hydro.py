import numpy as np
from copy import copy
import os
from functools import partial
from scipy.interpolate import CubicSpline,interp1d
from constants_list import *
from radio_analysis.solvers import euler, RK4
from .calcs import evolve_ODE, calc_dMdt_shell, calc_dvdt_shell, calc_dEintdt_shell, calc_shell_values
try:
	from scipy.integrate import simpson
except:
	from scipy.integrate import simps as simpson
    
def SNejecta_profile(E,Mej,delta=1,n=10):
    g_to_n = 1./(4.*np.pi*(n-delta)) * (2.*(5.-delta)*(n-5.)*E)**(n/2-1.5) / ((3.-delta)*(n-3)*Mej*Msun)**(n/2-2.5)
    v_t = np.sqrt(2.*(5.-delta)*(n-5.)*E / ((3.-delta)*(n-3)*Mej*Msun))
    return g_to_n,v_t

def Mej(g_to_n,shock_velocity,delta=1,n=10):
    return 4.*np.pi*g_to_n*(shock_velocity)**(3.-n) / (n-3.)

def evolve_ejectaCSM_shock(rho_prof_path,rundir,E,Mej,f_omega=1,R0_in=0.11,dt_in=1e-5,tf_in=1e2,t0_in=1.6e-5,max_step=50000,print_int=10000):
    # rho_prof_path: path to the file containing the density profile of the CSM
    # datadir: place to save output. should be a directory like 
    # f'E_{E:1.2E}_Mej_{Mej_calc:1.2E}_'+'rhoprof_'+f'{he_mass:1.3f}'+'M_'+model_dir.split('/')[-2].split('_')[1]+'_velfac_'+f'{vel_factor:0.2f}'+'_fomega_'+f'{f_omega:0.2f}''/'
    # could be where the density profile is stored, or a new directory if want to separate.
    # E: explosion energy in erg
    # Mej: ejecta mass in Msun
    # R0_in: inner radius of the CSM in AU 
    # f_omega: fraction of the sphere covered by the CSM (1 for spherical, <1 for aspherical)
    # dt_in: initial timestep in years
    # tf_in: final time in years
    # t0_in: initial time in years. change this if initial condition is too large of a velocity
    # max_step: maximum number of steps to take
    # print_int: print progress every print_int steps

    rho_prof_vs_r = np.load(rho_prof_path)
    print(f"Shock Energy = {E:1.2E} (erg), Ejecta mass = {Mej:1.2E} (Msun), t0 = {t0_in:1.2E} (yr), R0 = {R0_in:1.3E} (AU),tf = {tf_in:1.2E} (yr)","f_omega",f_omega)
    print('rho_prof',rho_prof_path)
    v_t, g_to_n = SNejecta_profile(E,Mej)
    print('v_t (km/s)',v_t/1e5,'g_to_n',g_to_n)

    dt0 = dt_in*secinyear 
    t0 = t0_in*secinyear
    R0 = R0_in*AU_cm #initial radius in cm
    final_age = tf_in*secinyear
    max_step = int(max_step)
    print_interval = int(print_int)
    save_interval = 1e-1*max_step

    if not os.path.exists(rundir): os.mkdir(rundir)
    if not os.path.exists(rundir+'/plots/'): os.mkdir(rundir+'/plots/')

    print('------------initial values: -------------')

    count = 0
    min_dt = 1e-30*secinyear #1e-10*secinyear
    dt_limit = 1e2*secinyear # choose an appropriate limit?
    stop = False
    dt_i = copy(dt0)
    t_i = copy(t0)

    #interpolate density profile
    radius_array = rho_prof_vs_r['r'] #cm
    rho_array = rho_prof_vs_r['rho'] #g/cm^3
    rho_vs_r = interp1d(radius_array,rho_array)

    print('Minimum radius (AU)', f'{np.amin(radius_array)/AU_cm:1.5E}', 'Maximum radius (AU)',f'{np.amax(radius_array)/AU_cm:1.5E}')


    #save values of times, shock parameters, etc.
    t_saved = np.array([],dtype=np.float64)
    dt_saved = np.array([],dtype=np.float64)
    rshock_saved = np.array([],dtype=np.float64)
    vshock_saved = np.array([],dtype=np.float64)
    Menc_saved = np.array([],dtype=np.float64)
    Eint_saved = np.array([],dtype=np.float64)

    t_saved = np.append(t_saved,t0) #s
    dt_saved = np.append(dt_saved,dt_i) #s
    rshock_saved = np.append(rshock_saved,R0) #cm

    rhoej_saved = np.array([],dtype=np.float64)
    rhoCSM_saved = np.array([],dtype=np.float64)
    dv_ejsh_saved = np.array([],dtype=np.float64)
    dv_shcsm_saved = np.array([],dtype=np.float64)

    while (t_i < final_age) and (count < max_step) and (dt_i > min_dt) and (stop == False):

        if count == 0:
            t_i = copy(t0) #sec
            # initial shock velocity is set to ejecta velocity at initial radius (should be ~0.5c for R0_in=0.11AU, t0~10^-5 yr)
            shock_velocity_i = R0 / t0 # cm/s
            # initial shell mass is set to mass occupied by power-law ejecta from initial velocity outwards
            # \int_(v_sh)^(+inf) 4pi*(t*v)**2*rho_ej*(t*dv) = 4pi*g^n*v^(1-n)/(n-1)
            Menc_i = Mej(g_to_n,shock_velocity_i) #g
            rshock_i = copy(R0) #cm
            Eint_i = 0.5*Menc_i*shock_velocity_i**2 #erg
            print('Menc_i',Menc_i,'E_int_i',Eint_i)

            vshock_saved = np.append(vshock_saved,shock_velocity_i) #cm/s
            Menc_saved = np.append(Menc_saved,Menc_i) #g
            Eint_saved = np.append(Eint_saved,Eint_i) #erg
        else:
            t_i = copy(t)
            Menc_i = copy(Menc)
            rshock_i = copy(rshock)
            shock_velocity_i = copy(shock_velocity)
            Eint_i = copy(Eint)

            
        #### need this for E=1.5e51/E=0.5e51. this works with t0_in = 1.6e-5
        if shock_velocity_i > 2e4*1e5: #km/s
            dt_i = 2.5e3*secinyear/shock_velocity_i 
        elif shock_velocity_i > 1e4*1e5:  #km/s
            if rshock_i > 1*AU_cm:
                dt_i = 5e4*secinyear/shock_velocity_i
            else:
                dt_i = 1e4*secinyear/shock_velocity_i
        #relax dt at large distances
        if rshock_i > 1e3*AU_cm:
            dt_i = 5e-2*secinyear

        if count % print_interval == 0:
            print("----------------Beginning of Step-----------------")
            print('count',count,'time (yr)',f'{t_i/secinyear:1.5E}','mass swept up (Msun)',Menc_i/Msun,'shock distance (AU)',rshock_i/AU_cm,
                'Eint',Eint_i, 'rho', rho_vs_r(rshock_i),
                    'shock velocity(km/s)',shock_velocity_i/1e5,'dt (yr)',dt_i/secinyear)

        rshock_of_t = lambda t: rshock_i + shock_velocity_i*(t-t_i) #linear approximation
        f_of_t_dvdt = partial(calc_dvdt_shell,Menc_of_t=Menc_i,Eint_of_t=Eint_i,radius_func=rshock_of_t,g_to_n=g_to_n,v_t=v_t,rho_func=rho_vs_r,f_omega=f_omega)
        shock_velocity,t = evolve_ODE(shock_velocity_i,t_i,dt_i,f_of_t_dvdt,solver=RK4)
        f_of_t_dEdt = partial(calc_dEintdt_shell,shock_velocity=shock_velocity_i, radius_func=rshock_of_t, rho_func=rho_vs_r,f_omega=f_omega)
        Eint_out, t = evolve_ODE(Eint_i,t_i,dt_i,f_of_t_dEdt,solver=RK4)
        Eint = max(Eint_out,0)
        f_of_t_dMdt = partial(calc_dMdt_shell,shock_velocity=shock_velocity_i,radius_func=rshock_of_t,g_to_n=g_to_n,v_t=v_t,rho_func=rho_vs_r,f_omega=f_omega) 
        Menc,t = evolve_ODE(Menc_i,t_i,dt_i,f_of_t_dMdt,solver=RK4) #calc_dMdt and RK4 are functions
        rshock = rshock_of_t(t)
        
        rhoej_of_t, dv_ejsh, rhoCSM_of_t, dv_shcsm = calc_shell_values(t, shock_velocity, rshock_of_t, g_to_n, v_t, rho_vs_r,f_omega=f_omega)

        ## calculate the density squared integral ahead of the shock for later
        r_gtr_rsh = np.logspace(np.log10(rshock), np.log10(np.amax(radius_array)), 100)
        rho1_gtr_rsh = [rho_vs_r(r)**2 for r in r_gtr_rsh]
        rho1sq_int_arr = np.append(rho1sq_int_arr, simpson(rho1_gtr_rsh, r_gtr_rsh))

        if rshock > 0.99*radius_array[0]:
            print("Reached outer radius")
            stop = True
        t_old = copy(t_i) #for next time, in case need to back up
        rshock_old = copy(rshock_i)

        if count % print_interval == 0:
            print(f'dt (yr): {dt_i/secinyear:1.2E}','dt limit (yr)',f'{dt_limit/secinyear:1.5E}')

        #save values-from end of step
        t_saved = np.append(t_saved,t)
        dt_saved = np.append(dt_saved,dt_i)
        rshock_saved = np.append(rshock_saved,rshock)
        Menc_saved = np.append(Menc_saved,Menc)
        vshock_saved = np.append(vshock_saved,shock_velocity) #cm/s
        Eint_saved = np.append(Eint_saved,Eint) #erg

        rhoej_saved = np.append(rhoej_saved,rhoej_of_t)
        rhoCSM_saved = np.append(rhoCSM_saved,rhoCSM_of_t)
        dv_ejsh_saved = np.append(dv_ejsh_saved,dv_ejsh)
        dv_shcsm_saved = np.append(dv_shcsm_saved,dv_shcsm)

        if count % save_interval == 0:
            np.savez(rundir+'/shock_data.npz',dts=dt_saved,times=t_saved,rsh_of_t=rshock_saved,vsh_of_t=vshock_saved,
                     Menc_of_t=Menc_saved,Eint_of_t=Eint_saved,
                     rhoCSM_of_t=rhoCSM_saved,rhoej_of_t=rhoej_saved,dv_ejsh=dv_ejsh_saved,dv_shcsm=dv_shcsm_saved)
        
        if (dt_i <= min_dt):
            print("stopping for min dt limit",dt_i)
        count+=1

    #save final values
    print("number of steps",count)
    np.set_printoptions(formatter={'float':'{:1.8E}'.format})
    if count < 10:
        print("dts (yr)",dt_saved/secinyear)
        print("times (yr)",t_saved/secinyear)
    else:
        inds = np.arange(0,count,int(np.floor(count/10)))
        print(inds)
        print("dts (yr)",dt_saved[inds]/secinyear)
        print("times (yr)",t_saved[inds]/secinyear)
    np.savez(rundir+'/shock_data.npz',dts=dt_saved,times=t_saved,rsh_of_t=rshock_saved,vsh_of_t=vshock_saved,
                     Menc_of_t=Menc_saved,Eint_of_t=Eint_saved,
                     rhoCSM_of_t=rhoCSM_saved,rhoej_of_t=rhoej_saved,dv_ejsh=dv_ejsh_saved,dv_shcsm=dv_shcsm_saved)
    return