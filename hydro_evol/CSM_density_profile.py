import numpy as np
import matplotlib.pyplot as plt
import mesa_reader as mr
import os
import sys
# import argparse
from scipy.interpolate import interp1d

script_dir = os.path.dirname(os.path.abspath(__file__))
style_path = os.path.join(script_dir, '../plot_styles.mplstyle_new')
plt.style.use(style_path)

def CSM_density_profile(model_dir,data_dir,E,vel_factor=1,f_omega=1,wind_scaling_factor=1):
    # model_dir: path to input model directory
    # data_dir: path to save output files
    # E: explosion energy in erg
    # vel_factor: factor by which to scale the orbital velocity of the CSM 
    # f_omega: fraction of the sphere covered by the CSM (1 for spherical, <1 for aspherical)
    # wind_scaling_factor: factor by which to scale the wind mass loss rate

    if model_dir=='':
        print('Please provide model directory')
        sys.exit()

    if np.any(np.array(os.listdir(model_dir))=='history2.data'):
        print('history2.data exists')
        # bh = mr. MesaData(model_dir+'binary_history.data')
        h1 = mr. MesaData(model_dir+'history.data')
        h2 = mr. MesaData(model_dir+'history2.data')
        h3 = mr. MesaData(model_dir+'history3.data')
        histories = {}
        histories['h1'] = h1
        histories['h2'] = h2    
        histories['h3'] = h3
        final_age=histories['h3'].star_age[-1]
        final_radius = np.average(histories['h3'].star_1_radius[-10:])
        #calculate Mej
        print('Final mass',histories['h3'].star_1_mass[-1],'Final radius',final_radius)
        Mej_calc = histories['h3'].star_1_mass[-1]-1.4
        print('Ejecta mass',Mej_calc)
    else:
        print('history2.data does not exist')
        # sys.exit()
        # bh = mr. MesaData(model_dir+'binary_history.data')
        h1 = mr. MesaData(model_dir+'history.data')
        histories = {}
        histories['h1'] = h1
        # histories['bh'] = bh
        final_age=histories['h1'].star_age[-1]
        final_radius = np.average(histories['h1'].star_1_radius[-10:])
        #calculate Mej
        print('Final mass',histories['h1'].star_1_mass[-1],'Final radius',final_radius)
        Mej_calc = histories['h1'].star_1_mass[-1]-1.4
        print('Ejecta mass',Mej_calc)


    he_mass=h1.star_mass[0]
    print(f'{he_mass:1.3f}',model_dir.split('/')[-2].split('_')[1])

    G = 6.67430e-8
    Msun = 1.989e33
    Lsun = 3.828e33
    Rsun = 6.96e10

    secinyear = 3.154e7
    secinday = 86400

    def log_Mdot_wind(L,Y,Z):
        return -11 + 1.29*np.log10(L) + 1.73*np.log10(Y) + 0.47*np.log10(Z) #Nugis and Lamers 2000 for He rich

    def rho_CSM(Mdot_CSM,v_CSM,ttil_cc,Rstar,f_omega=1):
        r_CSM = v_CSM*ttil_cc + Rstar
    #     print(r_CSM)
    #     print(Mdot_CSM)
        return Mdot_CSM/(f_omega*4*np.pi*r_CSM**2*v_CSM)


    ttil_cc_arr = np.array([])
    v_csm_arr = np.array([])
    v_esc_arr = np.array([])
    Mdot_csm = np.array([])
    Mdot_wind = np.array([])
    Rstar = np.array([])


    initial_mass=histories['h1'].star_1_mass[0]
    Yval = (histories['h1'].total_mass_he4/histories['h1'].star_mass)[0]
    Zval = (1 - (histories['h1'].total_mass_h1+histories['h1'].total_mass_he4)/histories['h1'].star_mass)[0]
    for i,key in enumerate(histories.keys()):
        indices_skip =np.where(histories[key].lg_mstar_dot_1!=1e-99)
        ttil_cc_arr = np.append(ttil_cc_arr,(final_age-histories[key].star_age)[indices_skip]*secinyear)
        v_csm_arr = np.append(v_csm_arr,histories[key].v_orb_2[indices_skip]*10**5) #convert km/s to cm/s. should be for the NS
        Mdot_csm = np.append(Mdot_csm,(10**histories[key].lg_mstar_dot_1[indices_skip])*Msun/secinyear)
        Mdot_wind = np.append(Mdot_wind, 10**log_Mdot_wind(10**histories[key].log_L,Yval,Zval)[indices_skip]*Msun/secinyear )
        v_esc_arr = np.append(v_esc_arr,
                        np.sqrt(2*G*histories[key].star_1_mass[indices_skip]*Msun/(histories[key].star_1_radius[indices_skip]*Rsun))) 
    # Mdot_csm += Mdot_wind
    Rstar = np.full(len(ttil_cc_arr),final_radius*Rsun)

    ttil_cc_arr = np.copy(ttil_cc_arr)[np.where(ttil_cc_arr>0)]
    v_csm_arr = vel_factor*np.copy(v_csm_arr)[np.where(ttil_cc_arr>0)]
    Mdot_csm = np.copy(Mdot_csm)[np.where(ttil_cc_arr>0)]
    Mdot_wind = wind_scaling_factor * np.copy(Mdot_wind)[np.where(ttil_cc_arr>0)]
    v_esc_arr = np.copy(v_esc_arr)[np.where(ttil_cc_arr>0)]
    Rstar = np.copy(Rstar)[np.where(ttil_cc_arr>0)]

    vel_shock = np.sqrt(2*E/Mej_calc/Msun) #in  cm/s
    print('Shock velocity [km/s]',vel_shock/1e5)

    save_density_dir =  data_dir +'/'+ f'E_{E:1.2E}_Mej_{Mej_calc:1.2E}_'+'rhoprof_'+f'{he_mass:1.3f}'+'M_'+model_dir.split('/')[-2].split('_')[1]+'_velfac_'+f'{vel_factor:0.2f}'+'_fomega_'+f'{f_omega:0.2f}''/'
    print('Creating save dir',save_density_dir)
    if not os.path.exists(save_density_dir):
        os.mkdir(save_density_dir)
    else:
        print('Directory exists')
        # sys.exit()

    plt.figure(figsize=(10,6))
    # final_age=histories['h3'].star_age[-1]
    initial_mass=histories['h1'].star_1_mass[0]

    for i,key in enumerate(histories.keys()):
        if i==0:
            label1=r'$\dot{M}$'
            label2=r'$\Delta M$'
        else:
            label1=label2=None
        indices_skip =np.where(histories[key].lg_mstar_dot_1!=1e-99)
        # print(indices_skip)
        plt.plot((final_age-histories[key].star_age)[indices_skip],10**histories[key].lg_mstar_dot_1[indices_skip],label=label1,c='tab:blue')
        plt.plot((final_age-histories[key].star_age)[indices_skip],initial_mass-histories[key].star_1_mass[indices_skip],label=label2,c='tab:orange')
        
    plt.yscale('log')
    plt.ylim(1e-7,)
    plt.xlim(1e6,1e-3)
    plt.axhline(1,ls='--',c='k')
    plt.xscale('log')
    plt.xlabel('Time til O burn (yr)')
    plt.ylabel('Mdot (Msun/yr)/Mass lost (Msun)')
    plt.legend()
    plt.savefig(save_density_dir+'Mdot_Massloss.png')

    print('CSM velocity (km/s)',np.average(v_csm_arr)/1e5,'without vel factor',np.average(v_csm_arr)/vel_factor/1e5,'covering fraction f_omega',f_omega)
    print('Escape velocity (km/s)', np.average(v_esc_arr)/1e5, 'Mdot wind (Msun/yr)',np.average(Mdot_wind)/Msun*secinyear )
    # v_csm_final = (f_omega*Mdot_wind*v_esc_arr + Mdot_csm*v_csm_arr)/(f_omega*Mdot_wind+Mdot_csm) # only f_omega of the wind contributes to momentum
    dist_CSM = v_csm_arr*ttil_cc_arr + Rstar
    density_CSM = rho_CSM(Mdot_csm,v_csm_arr,ttil_cc_arr,Rstar,f_omega=f_omega)

    dist_wind_orig = v_esc_arr*ttil_cc_arr + Rstar
    dist_wind = np.flip(np.logspace(np.log10(1.01*np.min(dist_wind_orig)),np.log10(0.99*np.max(dist_wind_orig)),5000))
    Mdot_wind_interp = interp1d(dist_wind_orig,Mdot_wind,kind='linear')
    v_esc_interp = interp1d(dist_wind_orig,v_esc_arr,kind='linear')
    ttil_cc_interp = interp1d(dist_wind_orig,ttil_cc_arr,kind='linear')
    Rstar_interp = interp1d(dist_wind_orig,Rstar,kind='linear')
    density_wind = rho_CSM(Mdot_wind_interp(dist_wind),v_esc_interp(dist_wind),ttil_cc_interp(dist_wind),Rstar_interp(dist_wind),f_omega=1) #just assume wind is isotropic but ejected with same velocity as CSM

    density_CSM_func = interp1d(dist_CSM,density_CSM)
    density_wind_func = interp1d(dist_wind,density_wind)

    dist_arr = np.flip(np.logspace(np.log10(1.01*max(np.min(dist_CSM),np.min(dist_wind))),np.log10(0.99*min(np.max(dist_CSM),np.max(dist_wind))),5000))
    print('distance array max/min',dist_arr[0],dist_arr[-1])
    density = density_CSM_func(dist_arr) + density_wind_func(dist_arr)
    # plt.figure()
    # plt.plot(ttil_cc_arr,v_csm_arr,label='CSM')
    # plt.plot(ttil_cc_arr,v_esc_arr,label='Wind')
    # # plt.plot(ttil_cc_arr,v_csm_final,label='Final')
    # plt.yscale('log')
    # plt.xscale('log')
    # plt.xlabel('Time til CC (sec)')
    # plt.legend()
    # # plt.xticks([1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3])
    # plt.ylabel(r'Velocity (cm/s)')
    # # plt.ylim(1e-41,)
    # plt.show()

    # plt.figure()
    # plt.plot(ttil_cc_arr/secinyear,density_CSM,label='CSM')
    # plt.plot(ttil_cc_interp(dist_wind)/secinyear,density_wind,label='Wind')
    # plt.plot(ttil_cc_arr/secinyear,density_CSM_func(dist_CSM),label='CSM',ls=':')
    # plt.plot(ttil_cc_interp(dist_wind)/secinyear,density_wind_func(dist_wind),label='Wind',ls=':')
    # plt.yscale('log')
    # plt.xscale('log')
    # plt.xlabel('Time til CC (yr)')
    # plt.legend()
    # # plt.xticks([1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3])
    # plt.ylabel(r'Density (g/cm$^3$)')
    # plt.ylim(1e-31,)
    # plt.show()


    plt.figure()

    plt.plot(dist_arr,density_CSM_func(dist_arr),label='CSM')
    plt.plot(dist_CSM,density_CSM,ls=':',color='grey')
    plt.plot(dist_arr,density_wind_func(dist_arr),label='Wind')
    plt.plot(dist_wind,density_wind,ls=':',color='black')
    plt.plot(dist_arr,density,label='Total')

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Distance (cm)')
    plt.legend()
    # plt.xticks([1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3])
    plt.ylabel(r'Density (g/cm$^3$)')
    plt.ylim(1e-30,1e-5)
    plt.savefig(save_density_dir+'wind+CSM_density_profile.png')
    # plt.show()


    plt.figure()
    plt.plot(dist_arr/vel_shock/secinyear,density)

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Approx. time since explosion (yr)')
    # plt.legend()
    plt.xticks([1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4])
    plt.ylabel(r'Density (g/cm$^3$)')
    plt.ylim(1e-30,1e-5)

    plt.savefig(save_density_dir+'density_profile.png')
    # plt.show()

    np.savez(save_density_dir+'density_prof.npz',rho=density,r=dist_arr) #save the density profile. radius is in cm