import numpy as np
#time integration
def RK4(u_n,t_n,h,f): #f is the input function during evolution
    k1 = h*f(u_n,t_n)
    k2 = h*f(u_n+k1/2,t_n+h/2)
    k3 = h*f(u_n+k2/2,t_n+h/2)
    k4 = h*f(u_n+k3,t_n+h) 
    return u_n+(k1+2*k2+2*k3+k4)/6.
def euler(u_n,t_n,h,f):
    return u_n + h*f(u_n,t_n) #using forward euler

#spatial differentiation

def df_dx(f,dx): #second order finite difference
    df_array = (np.roll(f,-1) - np.roll(f,1))/2
    
    df_array[0] = -1.5*f[0] + 2*f[1] - 0.5*f[2]
    df_array[-1] = 1.5*f[-1] - 2*f[-2] + 0.5*f[-3]
    # #do fourth order at boundaries
    # df_array[0] = -(25./12.)*f[0] + 4.*f[1] - 3.*f[2] + (4./3.)*f[3] - (1./4.)*f[4]
    # df_array[-1] = (25./12.)*f[-1] - 4*f[-2] + 3*f[-3] - (4./3.)*f[-4] + (1./4.)*f[-5]
    return df_array/dx

def elec_time_evol(dt,gamma_e,delta_energy,u_old,q_e_inj,P_cool):
    # u_new and u_old are arrays over the grid of gamma_e (length N_g)
    #gamma_e is spaced logarithmically from gamma_min to gamma_max
    #the grid spacing is given by delta_energy
    #u_new is calculated implicitly

    N_g = len(gamma_e)
    u_new = np.zeros(N_g)
    u_new[1:N_g-1] = (u_old[1:N_g-1] + q_e_inj[1:N_g-1]*dt + u_old[2:N_g]*P_cool[2:N_g]*dt/delta_energy[1:N_g-1])/(1 + P_cool[1:N_g-1]*dt/delta_energy[1:N_g-1])
    #do boundaries
    u_new[0] = (u_old[0] + q_e_inj[0]*dt + u_old[1]*P_cool[1]*dt/delta_energy[0])/(1 + P_cool[0]*dt/delta_energy[0])
    u_new[N_g-1] = (u_old[N_g-1] + q_e_inj[N_g-1]*dt)/(1 + P_cool[N_g-1]*dt/delta_energy[N_g-1])

    return u_new