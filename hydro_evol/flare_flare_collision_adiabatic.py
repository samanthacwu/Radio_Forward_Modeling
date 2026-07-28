# import argparse
import os
import math
import numpy as np
from matplotlib import pyplot as plt
from functools import partial
from radio_analysis.solvers import euler, RK4
from .calcs import calc_dvdt_flare, calc_dMdt_flare, calc_dEintdt_flare, evolve_ODE
try:
	from scipy.integrate import simpson
except:
	from scipy.integrate import simps as simpson
plt.style.use('tableau-colorblind10')
# color arrays for plotting
color_array = plt.rcParams['axes.prop_cycle'].by_key()['color']

# constants
c = 2.9989e10
Msun = 1.989e33
yr_to_sec = 3.156e7

# flare density profile
def rho_flare(r, t, A_norm, v_min, v_max,s=0.5): 
	# only exists for v_min < v=r/t < v_max
	if r/t < v_min or r/t > v_max:
		return 0.
	else:
		return A_norm/t**3 * (r/v_max/t)**(-2.*s-3.)

def evolve_flares(M_flares,delta_ts,v_min_c,v_max_c,s=0.5,data_dir='./evolve_spectrum/',dt_scale=1e-3):
	"""
		M_flares: (list of) flare masses in solar masses
		delta_ts: (list of) time between flares in years
		v_min_c: minimum flare velocity in units of c
		v_max_c: maximum flare velocity in units of c
		s: power-law index of flare density profile
		data_dir: directory to save output data
		dt_scale: set dt = dt_scale*(rsh/vsh) at each timestep to make sure shell doesn't expand too much at one timestep.
	"""
		
	M_flare_arr = np.array(M_flares) * Msun # flare masses in g
	deltat_arr = np.array(delta_ts) * yr_to_sec # time between flares in seconds
	v_min = v_min_c * c
	v_max = v_max_c * c

	# plot shell radius and mass as function of time
	plt.rcParams['font.size'] = 15
	fig, axes = plt.subplots(2, 2, figsize=(11, 10))
	ax1 = axes[0,0]
	ax2 = axes[0,1]
	ax3 = axes[1,0]
	ax4 = axes[1,1]
	ax1.set_title(r'$v_{\rm min}=%gc$' % (v_min/c))
	ax1.set_ylabel(r'shell radius [cm]')
	ax2.set_ylabel(r'shell velocity [$c$]')
	ax3.set_ylabel(r'shell mass [$M_\odot$]')
	ax4.set_ylabel(r'$dM_{\rm sh}/dt$ [$M_\odot$/yr]')
	for ax in [ax1,ax2,ax3,ax4]:
		ax.set_xlabel(r'time from collision [yr]')
		ax.set_xlim(1e-3, 1e2)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.grid(linestyle=':')
	ax2.set_yscale('linear')

	t_arr = np.linspace(0,100000)
	# ax2.plot(t_arr, np.ones(len(t_arr))*v_min/c, color='gray', linestyle='dashdot', label=r'$v_{\rm min}$')
	# ax2.plot(t_arr, np.ones(len(t_arr))*v_max/c, color='black', linestyle='dashdot', label=r'$v_{\rm max}$')
	for M_flare in M_flare_arr:
		ax3.plot(t_arr, np.ones(len(t_arr))*M_flare/Msun, linestyle='dashdot', label=r'$M_{\rm flare}=$'+f'{M_flare/Msun:1.2E} $M_\odot$')
	# ls_arr = ['solid', 'dashed', 'dotted']

	############# MAIN #################

	# sanity checks
	assert s>0 # s=0 will lead to no mass ejection. s<0 is unphysical.
	assert v_max > v_min
	for i,delta_t in enumerate(deltat_arr):
		for j,M_flare in enumerate(M_flare_arr):
			
			print(delta_t,M_flare)

			# obtain normalization of flare density profile
			A = M_flare*(2.*s)/(4.*math.pi*v_max**3) / ((v_max/v_min)**(2.*s) - 1.)

			# output kinetic energy. when s=1 the integration is a bit different
			if s == 1:
				print('kinetic energy: %g erg' % (2.*math.pi*A*v_max**5*math.log(v_max/v_min)))
			else:
				print('kinetic energy: %g erg' % (2.*math.pi*A*v_max**5/(2.-2.*s) * (1.-(v_max/v_min)**(2.*s-2.))))
			
			# initial conditions at collision
			t0 = v_min*delta_t / (v_max - v_min) # t defined as time from launch of later flare 
			rsh_0 = v_min*v_max*delta_t / (v_max - v_min) # collision radii
			# initial vsh set by pressure equilibrium
			rho_flare_1_init = rho_flare(rsh_0, t0+delta_t, A, v_min, v_max,s=s)
			rho_flare_2_init = rho_flare(rsh_0, t0, A, v_min, v_max,s=s)
			rho12_ratio = rho_flare_2_init/rho_flare_1_init
			vsh_0 = (v_min + v_max * math.sqrt(rho12_ratio)) / (1. + math.sqrt(rho12_ratio))
			v_fl_1_init = rsh_0/(t0+delta_t)
			v_fl_2_init = rsh_0/t0
			Msh_0 = 4.*math.pi*rsh_0**2*(rho_flare_2_init*(v_fl_2_init-vsh_0) + rho_flare_1_init*(vsh_0-v_fl_1_init)) * dt_scale*rsh_0/vsh_0
			#internal energy of shock
			Eint_0 = 0.5*Msh_0*vsh_0**2 # assumes KE=Eint

			# print('A', A, 'rho12_ratio', rho12_ratio, 'rho2_0', rho_flare_2_init, 'rho1_0', rho_flare_1_init, 'vsh_0 (c)', vsh_0/c)
			#vsh_0 = v_max 

			# initialize
			t = t0
			Msh = Msh_0
			rsh = rsh_0
			vsh = vsh_0
			Eint = Eint_0
			# initialize array
			t_arr = np.array([])
			dt_arr = np.array([])
			Msh_arr = np.array([]) 
			rsh_arr = np.array([])
			vsh_arr = np.array([])
			v1_arr = np.array([])
			v2_arr = np.array([])
			rho1_arr = np.array([])
			rho2_arr = np.array([])
			rho1sq_int_arr = np.array([])
			dMdt_arr = np.array([]) 
			Eint_arr = np.array([])
			print("initial t: %s sec" % t, 'initial Msh',Msh_0/Msun,'initial Eint',Eint_0)
			
			# solve shock propagation
			while t < 1000*t0:
				# make sure shell doesn't expand too much at one timestep
				dt = dt_scale*(rsh/vsh)
				# get current rho_flare
				rho_fl_1 = rho_flare(rsh, t+delta_t, A, v_min, v_max,s=s)
				rho_fl_2 = rho_flare(rsh, t, A, v_min, v_max,s=s)
				v_fl_1 = rsh/(t+delta_t)
				v_fl_2 = rsh/t

				# setup 
				rsh_old = rsh # backup of rsh
				t_old = t
				t_arr = np.append(t_arr, t)
				dt_arr = np.append(dt_arr, dt)
				#evolve
				dMdt_arr = np.append(dMdt_arr, 4.*math.pi*rsh**2*(rho_fl_2*(v_fl_2-vsh) + rho_fl_1*(vsh-v_fl_1)))
				f_of_t_dMdt = partial(calc_dMdt_flare,shock_velocity=vsh,radius=rsh,v1=v_fl_1,v2=v_fl_2,rho1=rho_fl_1,rho2=rho_fl_2)
				Msh, t = evolve_ODE(Msh,t_old,dt,f_of_t_dMdt,solver=RK4)
				# Msh += 4.*math.pi*rsh**2*(rho_fl_2*(v_fl_2-vsh) + rho_fl_1*(vsh-v_fl_1)) * (t-t_old)

				#calculate internal energy in forward shock, assuming dissipation is dominated by forward shock interacting with next flare
				f_of_t_dEdt = partial(calc_dEintdt_flare,shock_velocity=vsh,radius=rsh,v1=v_fl_1,rho1=rho_fl_1)
				Eint, t = evolve_ODE(Eint,t_old,dt,f_of_t_dEdt,solver=RK4)
				Eint = max(Eint,0) # make sure internal energy doesn't go negative
				Eint_arr = np.append(Eint_arr, Eint) 
				# Eint = 0

				f_of_t_dvdt = partial(calc_dvdt_flare,Menc_of_t=Msh,Eint_of_t=Eint,radius=rsh,v1=v_fl_1,v2=v_fl_2,rho1=rho_fl_1,rho2=rho_fl_2)
				vsh, t = evolve_ODE(vsh,t_old,dt,f_of_t_dvdt,solver=RK4)

				rsh += vsh*(t-t_old)
				
				# append
				# t_arr = np.append(t_arr, t)
				Msh_arr = np.append(Msh_arr, Msh)
				rsh_arr = np.append(rsh_arr, rsh)
				vsh_arr = np.append(vsh_arr, vsh)
				v1_arr = np.append(v1_arr, v_fl_1)
				v2_arr = np.append(v2_arr, v_fl_2)
				rho1_arr = np.append(rho1_arr, rho_fl_1)
				rho2_arr = np.append(rho2_arr, rho_fl_2)
				# integrated rho^2 for flare 1 outside rsh (used for free-free absorption)
				r_gtr_rsh = np.logspace(math.log10(rsh), math.log10(v_max*t+delta_t), 100)
				rho1_gtr_rsh = [rho_flare(r, t+delta_t, A, v_min, v_max,s=s)**2 for r in r_gtr_rsh]
				rho1sq_int_arr = np.append(rho1sq_int_arr, simpson(rho1_gtr_rsh, r_gtr_rsh))
			# plot evolution of shell parameters
			print('Msh initial',Msh_arr[0]/Msun)
			ax1.plot((t_arr-t0)/yr_to_sec, rsh_arr, color=color_array[i], label=r'$M_{\rm flare}=$'+f'{M_flare/Msun:1.2E} $M_\odot$, $\Delta t={delta_t/yr_to_sec:0.2f}$ yr')
			ax2.plot((t_arr-t0)/yr_to_sec, vsh_arr/c, color=color_array[i])
			ax2.plot((t_arr-t0)/yr_to_sec, (vsh_arr-v1_arr)/c, color=color_array[i],ls='--')
			ax2.plot((t_arr-t0)/yr_to_sec, (v2_arr-vsh_arr)/c, color=color_array[i],ls=':')
			ax3.plot((t_arr-t0)/yr_to_sec, Msh_arr/Msun, color=color_array[i])
			ax4.plot((t_arr-t0)/yr_to_sec, dMdt_arr*yr_to_sec/Msun, color=color_array[i])
			# save shell parameters
			save_dir = data_dir+f'shell_evolution_Mflare_{M_flare/Msun:1.0E}_deltat_{delta_t/yr_to_sec:0.2f}yr/'
			if not os.path.exists(save_dir):
				os.mkdir(save_dir)
			else:
				print('Directory exists')
			# np.savetxt('shell_evolution_deltat_%gyr.txt' % (delta_t/yr_to_sec), np.c_[t_arr, rsh_arr, vsh_arr, vsh_arr-v1_arr, v2_arr-vsh_arr, rho1_arr, rho2_arr,rho1sq_int_arr], header='Mfl=%gMsun, vmin=%gc, vmax=%gc, power-law-index=%g\ntime [s], rsh [cm], vsh [cm/s], Deltav_1 [cm/s], Deltav_2 [cm/s], rhofl_1 [g/cm3], rhofl_2 [g/cm3], int_rhosq_1dr [g^2/cm^5]' % (M_flare/Msun, v_min/c, v_max/c, p), fmt='%.8g')
			np.savez(save_dir+'/shock_data.npz',
				times=t_arr,
				rsh_of_t=rsh_arr,
				vsh_of_t=vsh_arr,
				deltav1_of_t=vsh_arr-v1_arr,
				deltav2_of_t=v2_arr-vsh_arr,
				rho1_of_t=rho1_arr,
				rho2_of_t=rho2_arr,
				int_rhofl_1_sq_dr=rho1sq_int_arr,
				Msh_of_t=Msh_arr,
				dMshdt_of_t=dMdt_arr,
				dt_arr=dt_arr
			)
			
			# plot v1, v2, vshell

			# fig0,ax0 = plt.subplots()
			# ax0.set_xlabel(r'time from collision [yr]')
			# ax0.set_ylabel(r'velocity [$c$]')
			# ax0.set_xlim(1e-3, 1e2)
			# ax0.set_xscale('log')
			# ax0.set_title('$\Delta t=%g$ yr' % (delta_t/yr_to_sec))
			# ax0.grid(linestyle=':')
			# ax0.plot((t_arr-t0)/yr_to_sec, vsh_arr/c, color=color_array[i], label=r'$v_{\rm sh}$')
			# ax0.plot((t_arr-t0)/yr_to_sec, v1_arr/c, linestyle='dashed', color=color_array[i+1], label=r'$v_1=r_{\rm sh}/(t+\Delta t)$')
			# ax0.plot((t_arr-t0)/yr_to_sec, v2_arr/c, linestyle='dotted', color=color_array[i+2], label=r'$v_2=r_{\rm sh}/t$')
			# ax0.legend()
			# fig0.tight_layout()
			# fig0.savefig('velocity_evolution_%gyr.pdf' % (delta_t/yr_to_sec))
	# finish
	ax1.legend()
	ax2.legend()
	ax3.legend()
	fig.tight_layout()
	fig.savefig('shell_evolution.pdf')

	return



