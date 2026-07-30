import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/samwu/codes/current_projects/Radio_Forward_Modeling/')
from scipy.interpolate import interp1d
# constants
c = 2.99792e10
Msun = 1.989e33
yr_to_sec = 3.154e7
plt.style.use('/Users/samwu/codes/current_projects/Radio_Forward_Modeling/plot_styles.mplstyle_new')

path_to_data = './evolve_spectrum_adiabatic/shell_evolution_Mflare_1E-02_pISM_0_rhoISM0_1E+03m_H/shock_data.npz'
path_to_data_old = './evolve_spectrum/shell_evolution_Mflare_1E-02_pISM_0_rhoISM0_1E+03m_H/shock_data.npz'

fig, ax = plt.subplots(3,1, figsize=(5, 10))

path_list = [path_to_data,path_to_data_old]
color_list = ['dodgerblue', 'crimson']

dMdt_funcs = []
dM_arrs = []
time_arrs = []
for path_to_data,color in zip(path_list,color_list):

    shock_data = np.load(path_to_data)
    times_yr = shock_data['times']/yr_to_sec
    t0 = times_yr[0]
    print(t0)
    ax[0].plot(times_yr, shock_data['rsh_of_t'],color=color)

    ax[1].plot(times_yr, shock_data['vsh_of_t']/c, label='shock velocity',color=color)
    # ax[1].plot(times_yr, shock_data['deltav1_of_t']/c, label='dv1',ls=':',color=color)
    ax[1].plot(times_yr, shock_data['deltav2_of_t']/c, label='dv2',ls='--',color=color)

    # dM = shock_data['dMshdt_of_t']*shock_data['dt_arr']/Msun
    # # print(dM[20:40])

    # ax[2].plot(times_yr, dM,color=color)
    ax[2].plot(times_yr, shock_data['rho2_of_t'], label='rho2',ls='--',color=color)
    ax[2].plot(times_yr, shock_data['rho_ism'], label='rho_ism',ls='-',color=color)

    # ax[3].plot(times_yr, shock_data['Msh_of_t']/Msun, label='shell mass',color=color)
    # ax[3].plot(times_yr,np.cumsum(dM),color=color,ls='--')
    # print(shock_data['Msh_of_t'][:20]/Msun)
    # print(times_yr[:10]-t0)
    # dMdt_funcs.append(interp1d(times_yr,dM))
    # ax[3].plot(times_yr,np.cumsum(interp1d(times_yr,dM)(times_yr)),color=color)
    # dM_arrs.append(np.cumsum(dM))
    # time_arrs.append(times_yr)

ax[0].set_xlabel('time (yr)')
ax[0].set_ylabel('shock radius (cm)')
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].set_xlim(1e-3,100)

ax[1].set_xlabel('time (yr)')
ax[1].set_ylabel('shock velocity (c)')
ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_xlim(1e-3,100)

ax[2].set_xlabel('time (yr)')
ax[2].set_ylabel('density (g/cm^3)')
# ax[2].set_ylabel(r'$dM_{\rm sh}$ (Msun)')
ax[2].set_yscale('log')
ax[2].set_xscale('log')
ax[2].set_xlim(1e-3,100)
# ax[2].set_ylim()

# ax[3].set_xlabel('time (yr)')
# ax[3].set_ylabel('enclosed mass (Msun)')
# ax[3].set_yscale('log')
# ax[3].set_xscale('log')
# ax[3].set_xlim(1e-1,100)

ax[1].legend()
ax[2].legend()


fig.tight_layout()
plt.savefig('./hydro_props.pdf',bbox_inches='tight')
plt.close()

# fig2,ax2 = plt.subplots()
# times_yr = np.logspace(-0.95,2,1000)
# print(time_arrs[0][0],time_arrs[0][-1],len(time_arrs[0]))
# print(time_arrs[1][0],time_arrs[1][-1],len(time_arrs[1]))
# ax2.plot(times_yr, np.cumsum(dMdt_funcs[0](times_yr)),color='dodgerblue')
# ax2.plot(times_yr, np.cumsum(dMdt_funcs[1](times_yr)),color='crimson')
# # ax2.plot(time_arrs[0], np.cumsum(dMdt_funcs[0](time_arrs[0])),color='dodgerblue',ls=':')
# # ax2.plot(time_arrs[1], np.cumsum(dMdt_funcs[1](time_arrs[1])),color='crimson',ls=':')

# # print(np.cumsum(dMdt_funcs[0](times_yr)))
# # print(np.cumsum(dMdt_funcs[0](time_arrs[0])))
# ax2.plot(time_arrs[0], dM_arrs[0],color='dodgerblue',ls='--')
# ax2.plot(time_arrs[1], dM_arrs[1],color='crimson',ls='--')
# ax2.set_yscale('log')
# ax2.set_xscale('log')
# plt.show()