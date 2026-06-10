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

def calc_shell_values(t, shock_velocity, radius_func, g_to_n, v_t, rho_func,f_omega=1):
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
    return rhoej_of_t,(v_ej-shock_velocity), rhoCSM_of_t, (shock_velocity-v_CSM)

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

def calc_dvdt_shell(shock_velocity, t, Menc_of_t, Eint_of_t, radius_func, g_to_n, v_t, rho_func,f_omega=1):
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
    return f_omega*4*np.pi*radius**2*(rhoej_of_t*(v_ej-shock_velocity)**2 - rhoCSM_of_t*(shock_velocity-v_CSM)**2) / Menc_of_t + 2*Eint_of_t/radius/Menc_of_t

## add energy conservation to account for adiabatic PdV work and radiative cooling losses
def calc_dEintdt_shell(Eint, t, shock_velocity, radius_func, rho_func,f_omega=1):
    # n = 10.
    # delta = 1.
    gamma = 5./3.
    radius = radius_func(t)
    # print('radius',radius,'t',t,'v',shock_velocity)
    rhoCSM_of_t = rho_func(radius)
    # NOTE: v_CSM is assumed to be zero. 
    # It's good approximation for our binary case, but we'll need to set a finite velocity (profile) for the TDE flare calculation.
    v_CSM = 0.0

    m_p = 1.67e-24 # g
    k_b = 1.38e-16 # erg/K
    mu = 4./3. # helium rich gas values (assuming fully ionized)
    T_down = (3.*mu*m_p/16./k_b)*(shock_velocity-v_CSM)**2
    rho_down = 4*rhoCSM_of_t
    Z = 2 #helium rich gas values
    n_e = rho_down/m_p/2. #helium rich gas values
    n_i = rho_down/m_p/4. #helium rich gas values
    g_B = 1.2
    eps_ff = 1.4*10**(-27.) * T_down**0.5 * Z**2 * n_e * n_i * g_B
    vol_fac = ((gamma-1)/(gamma+1)) * f_omega * 4*np.pi*radius**3/3
    # print('T_down',T_down,'eps_ff',eps_ff,'n_e',n_e,'n_i',n_i,'rho_down',rho_down)
    # print('E_int',Eint,'v',shock_velocity)
    # print('term 1', f_omega*4*np.pi*radius**2*(2*(gamma+1)**(-2.))*rhoCSM_of_t*(shock_velocity-v_CSM)**3, 'term2', 2*Eint*shock_velocity/radius, 'term3', eps_ff*vol_fac)
    return f_omega*4*np.pi*radius**2*(2*(gamma+1)**(-2.))*rhoCSM_of_t*(shock_velocity-v_CSM)**3 - 2*Eint*shock_velocity/radius - eps_ff*vol_fac

#for flares:
def calc_dMdt_flare(Menc_of_t, t, shock_velocity, radius, v1, v2, rho1, rho2,f_omega=1):
    #rho1, v1 are for forward shock; rho2,v2 are for backward shock
    return f_omega*4*np.pi*radius**2*(rho2*(v2-shock_velocity) + rho1*(shock_velocity-v1))

def calc_dvdt_flare(shock_velocity, t, Menc_of_t, Eint_of_t, radius, v1, v2, rho1, rho2,f_omega=1):
    return f_omega*4*np.pi*radius**2*(rho2*(v2-shock_velocity)**2 - rho1*(shock_velocity-v1)**2) / Menc_of_t + 2*Eint_of_t/radius/Menc_of_t


## add energy conservation to account for adiabatic PdV work and radiative cooling losses
def calc_dEintdt_flare(Eint, t, shock_velocity, radius, v1, rho1,f_omega=1):

    gamma = 5./3.

    m_p = 1.67e-24 # g
    k_b = 1.38e-16 # erg/K
    mu = 0.5 #hydrogen rich gas values
    T_down = (3.*mu*m_p/16./k_b)*(shock_velocity-v1)**2
    rho_down = 4*rho1
    Z = 1 #hydrogen rich gas values
    n_e = rho_down/m_p #hydrogen rich gas values
    n_i = rho_down/m_p #hydrogen rich gas values
    g_B = 1.2
    eps_ff = 1.4*10**(-27.) * T_down**0.5 * Z**2 * n_e * n_i * g_B
    vol_fac = ((gamma-1)/(gamma+1)) * f_omega * 4*np.pi*radius**3/3
    # print('T_down',T_down,'eps_ff',eps_ff,'n_e',n_e,'n_i',n_i,'rho_down',rho_down)
    # print('E_int',Eint,'v',shock_velocity)
    # print('term 1', f_omega*4*np.pi*radius**2*(2*(gamma+1)**(-2.))*rho1*(shock_velocity-v1)**3, 'term2', 2*Eint*shock_velocity/radius, 'term3', eps_ff*vol_fac)
    return f_omega*4*np.pi*radius**2*(2*(gamma+1)**(-2.))*rho1*(shock_velocity-v1)**3 - 2*Eint*shock_velocity/radius - eps_ff*vol_fac
