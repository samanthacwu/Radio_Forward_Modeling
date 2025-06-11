import numpy as np

Msun = 1.3271244e26/6.67430e-8 #g
Rsun = 6.957e10 #cm
GG = 6.67430e-8 # grav. constant (cm^3sg^-1)
Lsun = 3.839e33 #erg/s
secinyear=3.1536*10**7
secinday=86400

def evolve_ODE(un,tn,dt,f,solver): #f is the function that evolves semimajor axis. solver is RK4
    #this is a general function that can evolve the ode using a solver (probably using RK4)
    u_next= solver(un,tn,dt,f) 
    t_next = tn+dt #update to be t(n) + dt
    return u_next, t_next

    
#when input into solver, use partial(calc_dMdt,Menc,t,radius=rshock_i,velocity=shock_velocity_i,rho_func=rho_vs_r)
def calc_dMdt(Menc_of_t,t,radius_func,energy,rho_func):#, #M1, q, and arrays are constants. t in seconds, a in cm
    velocity =  np.sqrt(2*energy/Menc_of_t)
    radius = radius_func(t)
    rho_of_t = rho_func(radius)
    # print("density",rho_of_t)
    # print("velocity",velocity)
    # print("radius",radius)
    # print(4*np.pi*rho_of_t*radius**2*velocity)
    return 4*np.pi*rho_of_t*radius**2*velocity


# (Added by DT) assuming ejecta-CSM interface as thin shell but resolving the ejecta,
# the bottom two functions solve the increase of velocity and mass of the shell per unit time
def calc_dMdt_shell(Menc_of_t, t, shock_velocity, radius_func, g_to_n, v_t, rho_func,f_omega=1):
    n = 10.
    delta = 1.
    radius = radius_func(t)
    rhoCSM_of_t = rho_func(radius)
    # NOTE: v_CSM is assumed to be zero. 
    # It's good approximation for our binary case, but we'll need to set a finite velocity (profile) for the TDE flare calculation.
    v_CSM = 0.0
    # calculate ejecta density at r=radius
    v_ej = radius/t
    rhoej_at_vt = g_to_n/v_t**n/t**3
    if v_ej > v_t:
        rhoej_of_t = rhoej_at_vt * (v_ej/v_t)**(-n)
    else:
        rhoej_of_t = rhoej_at_vt * (v_ej/v_t)**(-delta)
    return f_omega*4*np.pi*radius**2*(rhoej_of_t*(v_ej-shock_velocity) + rhoCSM_of_t*(shock_velocity-v_CSM))

def calc_dvdt_shell(shock_velocity, t, Menc_of_t, radius_func, g_to_n, v_t, rho_func,f_omega=1):
    n = 10.
    delta = 1.
    radius = radius_func(t)
    # print('radius',radius,'t',t,'v',shock_velocity)
    rhoCSM_of_t = rho_func(radius)
    # NOTE: v_CSM is assumed to be zero. 
    # It's good approximation for our binary case, but we'll need to set a finite velocity (profile) for the TDE flare calculation.
    v_CSM = 0.0
    # calculate ejecta density at r=radius
    v_ej = radius/t
    rhoej_at_vt = g_to_n/v_t**n/t**3
    if v_ej > v_t:
        rhoej_of_t = rhoej_at_vt * (v_ej/v_t)**(-n)
    else:
        rhoej_of_t = rhoej_at_vt * (v_ej/v_t)**(-delta)
    # print('result',4*np.pi*radius**2*(rhoej_of_t*(v_ej-shock_velocity)**2 - rhoCSM_of_t*(shock_velocity-v_CSM)**2) / Menc_of_t)
    return f_omega*4*np.pi*radius**2*(rhoej_of_t*(v_ej-shock_velocity)**2 - rhoCSM_of_t*(shock_velocity-v_CSM)**2) / Menc_of_t
