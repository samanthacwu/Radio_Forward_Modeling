import math
import numpy as np
from matplotlib import pyplot as plt
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

# model parameters
# current assumption is that the two flares are identical
p = 0.5
v_max = 0.1 * c
v_min = 0.01 * c
M_flare = 1e-1 * Msun

# deltat_arr = np.array([0.1, 0.3, 1.]) * yr_to_sec # vary delta_t
# deltat_arr = np.array([1.]) * yr_to_sec # vary delta_t
rho_ism_arr = np.array([1.,10.,100.,1000.])*1e-24 # g/cm^3
#np.array([1,10,100])* 1e-24 # g/cm^3

# flare density profile
def rho_flare(r, t, A_norm, v_min, v_max): 
	# only exists for v_min < v=r/t < v_max
	if r/t < v_min or r/t > v_max:
		return 0.
	else:
		return A_norm/t**3 * (r/v_max/t)**(-2.*p-3.)


# plot shell radius and mass as function of time
plt.rcParams['font.size'] = 15
fig, axes = plt.subplots(2, 2, figsize=(11, 10))
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]
ax1.set_title(r'$M_{\rm flare}=%gM_\odot$, $v_{\rm min}=%gc$' % (M_flare/Msun, v_min/c))
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
ax2.plot(t_arr, np.ones(len(t_arr))*v_min/c, color='gray', linestyle='dashdot', label=r'$v_{\rm min}$')
ax2.plot(t_arr, np.ones(len(t_arr))*v_max/c, color='black', linestyle='dashdot', label=r'$v_{\rm max}$')
ax3.plot(t_arr, np.ones(len(t_arr))*M_flare/Msun, linestyle='dashdot', label=r'$M_{\rm flare}$')
ls_arr = ['solid', 'dashed', 'dotted', '-.']

############# MAIN #################

# sanity checks
assert p>0 # p=0 will lead to no mass ejection. p<0 is unphysical.
assert v_max > v_min

# obtain normalization of flare density profile
A = M_flare*(2.*p)/(4.*math.pi*v_max**3) / ((v_max/v_min)**(2.*p) - 1.)

# output kinetic energy. when p=1 the integration is a bit different
if p == 1:
	print('kinetic energy: %g erg' % (2.*math.pi*A*v_max**5*math.log(v_max/v_min)))
else:
	print('kinetic energy: %g erg' % (2.*math.pi*A*v_max**5/(2.-2.*p) * (1.-(v_max/v_min)**(2.*p-2.))))

for i, rho_ism in enumerate(rho_ism_arr):
	
	# initial conditions at collision
	t0 = 0.01 * yr_to_sec #v_min*delta_t / (v_max - v_min) # t defined as time from launch of later flare 
	Msh_0 = 0.0
	rsh_0 = v_max * t0 
	#v_min*v_max*delta_t / (v_max - v_min) # collision radii
	# initial vsh set by pressure equilibrium
	rho_flare_1_init = rho_ism #rho_flare(rsh_0, t0+delta_t, A, v_min, v_max)
	rho_flare_2_init = rho_flare(rsh_0, t0, A, v_min, v_max)
	rho12_ratio = rho_flare_2_init/rho_ism
	vsh_0 = (v_min + v_max * math.sqrt(rho12_ratio)) / (1. + math.sqrt(rho12_ratio))
	#vsh_0 = v_max 

	# initialize
	t = t0
	Msh = Msh_0
	rsh = rsh_0
	vsh = vsh_0
	# initialize array
	t_arr = np.array([])
	Msh_arr = np.array([]) 
	rsh_arr = np.array([])
	vsh_arr = np.array([])
	# v1_arr = np.array([])
	v2_arr = np.array([])
	rho1_arr = np.array([])
	rho2_arr = np.array([])
	rho1sq_int_arr = np.array([])
	dMdt_arr = np.array([]) 
	print("initial t: %s sec" % t)
	# solve shock propagation
	while t < 1000*t0:
		# make sure shell doesn't expand too much at one timestep
		dt = 1e-3*(rsh/vsh)
		# get current rho_flare
		# rho_fl_1 = rho_flare(rsh, t+delta_t, A, v_min, v_max)
		rho_fl_2 = rho_flare(rsh, t, A, v_min, v_max)
		# v_fl_1 = rsh/(t+delta_t)
		v_fl_2 = rsh/t
		# evolve
		rsh_old = rsh # backup of rsh
		dMdt_arr = np.append(dMdt_arr, 4.*math.pi*rsh**2*(rho_fl_2*(v_fl_2-vsh) + rho_ism*vsh))
		Msh += 4.*math.pi*rsh**2*(rho_fl_2*(v_fl_2-vsh) + rho_ism*vsh) * dt
		rsh += vsh * dt
		vsh += (4.*math.pi*rsh_old**2/Msh) * (rho_fl_2*(v_fl_2-vsh)**2 - rho_ism*vsh**2) * dt
		t += dt
		# append
		t_arr = np.append(t_arr, t)
		Msh_arr = np.append(Msh_arr, Msh)
		rsh_arr = np.append(rsh_arr, rsh)
		vsh_arr = np.append(vsh_arr, vsh)
		# v1_arr = np.append(v1_arr, v_fl_1)
		v2_arr = np.append(v2_arr, v_fl_2)
		rho1_arr = np.append(rho1_arr, rho_ism)
		rho2_arr = np.append(rho2_arr, rho_fl_2)
		# integrated rho^2 for flare 1 outside rsh (used for free-free absorption)
		r_gtr_rsh = np.logspace(math.log10(rsh), math.log10(v_max*t), 100)
		rho1_gtr_rsh = [rho_ism**2 for r in r_gtr_rsh]
		# rho1_gtr_rsh = [rho_flare(r, t, A, v_min, v_max)**2 for r in r_gtr_rsh]
		rho1sq_int_arr = np.append(rho1sq_int_arr, simpson(rho1_gtr_rsh, r_gtr_rsh))
	# plot evolution of shell parameters
	ax1.plot((t_arr-t0)/yr_to_sec, rsh_arr, ls=ls_arr[i], color=color_array[i], label=r'$\rho_{\rm ISM} =%g \times 10^{-24}$ g/cm$^3$ ' % (rho_ism/1e-24))
	ax2.plot((t_arr-t0)/yr_to_sec, vsh_arr/c, ls=ls_arr[i], color=color_array[i])
	ax3.plot((t_arr-t0)/yr_to_sec, Msh_arr/Msun, ls=ls_arr[i], color=color_array[i])
	ax4.plot((t_arr-t0)/yr_to_sec, dMdt_arr*yr_to_sec/Msun, ls=ls_arr[i], color=color_array[i])
	# save shell parameters
	np.savetxt('shell_evolution_rhoISM_%ge-24_cgs.txt' % (rho_ism/1e-24), np.c_[t_arr, rsh_arr, vsh_arr, v2_arr-vsh_arr, rho1_arr, rho2_arr, rho1sq_int_arr], header='Mfl=%gMsun, vmin=%gc, vmax=%gc, power-law-index=%g\ntime [s], rsh [cm], vsh [cm/s], Deltav_2 [cm/s], rho_ISM [g/cm3], rhofl_2 [g/cm3], int_rhosq_1dr [g^2/cm^5]' % (M_flare/Msun, v_min/c, v_max/c, p), fmt='%.8g')
	# plot v1, v2, vshell

	fig0,ax0 = plt.subplots()
	ax0.set_xlabel(r'time from collision [yr]')
	ax0.set_ylabel(r'velocity [$c$]')
	ax0.set_xlim(1e-3, 1e2)
	ax0.set_xscale('log')
	ax0.set_title(r'$\rho_{\rm ISM} =%g \times 10^{-24}$ g/cm$^3$' % (rho_ism/1e-24))
	ax0.grid(linestyle=':')
	ax0.plot((t_arr-t0)/yr_to_sec, vsh_arr/c, color=color_array[i], label=r'$v_{\rm sh}$')
	# ax0.plot((t_arr-t0)/yr_to_sec, v1_arr/c, linestyle='dashed', color=color_array[i+1], label=r'$v_1=r_{\rm sh}/(t+\Delta t)$')
	ax0.plot((t_arr-t0)/yr_to_sec, v2_arr/c, linestyle='dotted', color=color_array[i+2], label=r'$v_2=r_{\rm sh}/t$')
	ax0.legend()
	fig0.tight_layout()
	fig0.savefig('velocity_evolution_rhoISM_%ge-24_cgs.pdf' % (rho_ism/1e-24))
# finish
ax1.legend()
ax2.legend()
ax3.legend()
fig.tight_layout()
fig.savefig('shell_evolution.pdf')



