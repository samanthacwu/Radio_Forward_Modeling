import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from hydro_evol.model import Model
from hydro_evol.constants_list import *
from .spectrum_evolution import gamma_e_func
from .radio_tools import *

def analyze_multiwavelength_spectrum(simtype,data_dir,freq_in,T_e_csm,dNdgamma_dir_in='',f_omega=1):
    """
    Analyze multiwavelength spectrum evolution from simulation data.
        
    Parameters:
    simtype (str): Type of simulation ('flare_flare' or 'flare_ism').
    data_dir: Directory containing shock evolution data.
    freq_in: Frequency to plot in GHz.
    Te_csm (float): Temperature of CSM in Kelvin.
    dNdgamma_dir: Directory containing electron spectrum data.
    f_omega (float): Covering fraction of CSM (default is 1).

    """
    if dNdgamma_dir_in=='':
        dNdgamma_over_dir = data_dir
    else:
        dNdgamma_over_dir = dNdgamma_dir_in

    if simtype=='flare_flare':
        flare_list=['fwd','bwd']
    elif simtype=='flare_ism':
        flare_list=['fwd']

    m = Model(data_dir+'shock_data.npz',simtype=simtype)
    m.generate_ND_interp_funcs(simtype)
    m.generate_interp_funcs(simtype)

    #set gamma_e values
    gamma_e_vals, delta_gamma_e, d_ln_gamma = gamma_e_func(m.gamma_min,m.N_g,gamma_max=1e8)
    #set gamma_ph values
    gamma_ph_vals = gamma_ph_func(gamma_ph_min=1.5e-18,gamma_ph_max=1.5,N_ph=256)
    nu_ph_vals = nu_ph_func(gamma_ph_vals)  # in Hz
    freq_to_plot = freq_in*1e9 # Convert to Hz from GHz

    for flarenum,flare in enumerate(flare_list):
        if flare == 'fwd':
            B_ND_func = m.B_fwd_ND_func
            n0_t0 = m.n0_fwd[0]
            
        elif flare == 'bwd':
            B_ND_func = m.B_bwd_ND_func
            n0_t0 = m.n0_bwd[0]
            
        dNdgamma_dir = dNdgamma_over_dir + f'/{flare}_flare'+str(flarenum+1)+'/'
        dts = np.load(dNdgamma_dir+'/dts.npy')
        tvals=np.load(dNdgamma_dir+'/times.npy')
        yvals=np.load(dNdgamma_dir+'/yvals.npz')['arr_0'].reshape((len(tvals),len(gamma_e_vals)))
        dNdgamma_vals = yvals * n0_t0 * m.rsh_t0**3
        print('normalization',n0_t0 * m.rsh_t0**3)
    
        # find closest frequency index to freq_to_plot
        index_atfreq=np.argmin(np.abs(nu_ph_vals-freq_to_plot)/freq_to_plot )
        nu_in = nu_ph_vals[index_atfreq] # 3 GHz
        print('nu_in',nu_in,'given frequency',freq_to_plot)

        plt.figure()
        colors = plt.cm.viridis(np.linspace(0,1,15))
        count = 0
        tau_ff_atfreq= np.zeros_like(tvals)
        Lnu_atfreq = np.zeros_like(tvals)
        tau_ssa_atfreq = np.zeros_like(tvals)
        times_out = np.zeros_like(tvals)
        j=0
        for i,t in enumerate(tvals): #in units of tdyn_t0
            if i % 10 ==0:
                # print(i)
                tau_ff_test = tau_ff(m.int_rhofwd_sq_dr_func(t*m.tdyn_t0),nu_ph_vals,T_e_csm=T_e_csm)
                tau_ff_atfreq[j] = tau_ff_test[index_atfreq]
                Lnu_test,tau_ssa_test = emission_absorption_at_time(tvals[i],dNdgamma_vals[i],gamma_e_vals,delta_gamma_e,m.N_g,nu_ph_vals,
                                                                    m.rsh_ND_func,m.rsh_t0,B_ND_func,f_omega=f_omega)
                Lnu_atfreq[j] = Lnu_test[index_atfreq]
                tau_ssa_atfreq[j] = tau_ssa_test[index_atfreq]
                times_out[j] = t*m.tdyn_t0
                j+=1

        tau_ff_atfreq= np.trim_zeros(tau_ff_atfreq)
        Lnu_atfreq = np.trim_zeros(Lnu_atfreq)
        tau_ssa_atfreq = np.trim_zeros(tau_ssa_atfreq)
        times_out = np.trim_zeros(times_out)  #in seconds now
        # print(tau_ssa_atfreq)
        Lnu_abs_atfreq = Lnu_atfreq*np.exp(-tau_ff_atfreq)*(1-np.exp(-tau_ssa_atfreq))/tau_ssa_atfreq
        print(np.amax(Lnu_atfreq),np.amin(Lnu_atfreq))
        plt.plot(times_out/secinyear,Lnu_atfreq*np.exp(-tau_ff_atfreq)*(1-np.exp(-tau_ssa_atfreq))/tau_ssa_atfreq)
        #,color=colors[int((count-7000)/600)])
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e-2,2e2)
        plt.ylim(1e10,1e30)
        # plt.xscale('log')
        plt.xlabel('Time (yr)')
        plt.ylabel(r'$L_{\nu}$ (erg/s/Hz)')
        plt.title(f'Emission at {freq_to_plot/1e9:1.2f} GHz')
        # plt.ylim(1e20,1e30)
        # plt.xlim(1e-2,1000)
        # plt.axvspan(2,4,color='gray',alpha=0.5)
        # plt.axhspan(1e26,1e29,color='gray',alpha=0.5)
        # plt.axvspan(5,20,color='gray',alpha=0.5)
        plt.savefig(f'Radio_curve_{freq_to_plot/1e9:1.2f}GHz_{flarenum+1}.png',dpi=300,transparent=False,facecolor='white')

        np.savez(dNdgamma_dir+f'Lnu_{freq_to_plot/1e9:1.2f}GHz_sparse_Te_{T_e_csm:1.1E}.npz',times_out=times_out, #in seconds
                Lnu_atfreq=Lnu_atfreq,tau_ff_atfreq=tau_ff_atfreq,tau_ssa_atfreq=tau_ssa_atfreq,
                times_yr=times_out/secinyear,Lnu_abs_atfreq=Lnu_abs_atfreq)

def analyze_SED(simtype,data_dir,freq_in,T_e_csm,
                epoch_list,dNdgamma_dir_in='',f_omega=1,tval_ind=1000,SED_interval=100):
    """
    Analyze multiwavelength spectrum evolution from simulation data.
        
    Parameters:
    simtype (str): Type of simulation ('flare_flare' or 'flare_ism').
    data_dir: Directory containing shock evolution data.
    freq_in: Frequency to plot in GHz.
    Te_csm (float): Temperature of CSM in Kelvin.
    epoch_list: list of epochs to save SED data
    dNdgamma_dir: Directory containing electron spectrum data.
    f_omega (float): Covering fraction of CSM (default is 1).
    tval_ind, SED_interval control the time step for plotting SEDs.

    """
    from matplotlib.colors import LogNorm

    if dNdgamma_dir_in=='':
        dNdgamma_over_dir = data_dir
    else:
        dNdgamma_over_dir = dNdgamma_dir_in

    if simtype=='flare_flare':
        flare_list=['fwd','bwd']
    elif simtype=='flare_ism':
        flare_list=['fwd']

    m = Model(data_dir+'shock_data.npz',simtype=simtype)
    m.generate_ND_interp_funcs(simtype)
    m.generate_interp_funcs(simtype)

    #set gamma_e values
    gamma_e_vals, delta_gamma_e, d_ln_gamma = gamma_e_func(m.gamma_min,m.N_g,gamma_max=1e8)
    #set gamma_ph values
    gamma_ph_vals = gamma_ph_func(gamma_ph_min=1.5e-18,gamma_ph_max=1.5,N_ph=256)
    nu_ph_vals = nu_ph_func(gamma_ph_vals)  # in Hz
    freq_to_plot = freq_in*1e9 # Convert to Hz from GHz

    for flarenum,flare in enumerate(flare_list):
        if flare == 'fwd':
            B_ND_func = m.B_fwd_ND_func
            n0_t0 = m.n0_fwd[0]
            
        elif flare == 'bwd':
            B_ND_func = m.B_bwd_ND_func
            n0_t0 = m.n0_bwd[0]
            
        dNdgamma_dir = dNdgamma_over_dir + f'/{flare}_flare'+str(flarenum+1)+'/'
        dts = np.load(dNdgamma_dir+'/dts.npy')
        tvals=np.load(dNdgamma_dir+'/times.npy')
        yvals=np.load(dNdgamma_dir+'/yvals.npz')['arr_0'].reshape((len(tvals),len(gamma_e_vals)))
        dNdgamma_vals = yvals * n0_t0 * m.rsh_t0**3
        print('normalization',n0_t0 * m.rsh_t0**3)
    
        # find closest frequency index to freq_to_plot
        index_atfreq=np.argmin(np.abs(nu_ph_vals-freq_to_plot)/freq_to_plot )
        nu_in = nu_ph_vals[index_atfreq] # 3 GHz
        print('nu_in',nu_in,'given frequency',freq_to_plot)

##### Uncomment if want Peak Data ####
        fig=plt.figure()

        norm = LogNorm(vmin=tvals[tval_ind]*m.tdyn_t0/secinyear, vmax=tvals[-1]*m.tdyn_t0/secinyear)
        count = 0
        Lnu_pkvals = []
        t_pkvals = []
        nu_pkvals = []
        
        for t in tvals:
            if count % SED_interval ==0:
                print(count,int(count/SED_interval), t)
                Lnu_test,tau_ssa_test = emission_absorption_at_time(tvals[count],dNdgamma_vals[count],gamma_e_vals,delta_gamma_e,m.N_g,nu_ph_vals,
                                                                    m.rsh_ND_func,m.rsh_t0,B_ND_func,f_omega=f_omega)
                tau_ff_test = tau_ff(m.int_rhofwd_sq_dr_func(tvals[count]*m.tdyn_t0),nu_ph_vals,T_e_csm=T_e_csm)
                Lnu_spectrum = Lnu_test*np.exp(-tau_ff_test)*(1-np.exp(-tau_ssa_test))/tau_ssa_test
                argmax,Lnumax = (np.argmax(Lnu_spectrum),np.amax(Lnu_spectrum))
                
                Lnu_pkvals.append(Lnumax)
                nu_pkvals.append(nu_ph_vals[argmax])
                t_pkvals.append(t*m.tdyn_t0/secinyear)
                plt.plot(nu_ph_vals/1e9,Lnu_spectrum,color=plt.cm.viridis(norm(t*m.tdyn_t0/secinyear)))
                
            count +=1
        plt.plot(nu_ph_vals/1e9,1e30*(nu_ph_vals/1e9)**(-1.5),color='black',ls='--',label=r'$\nu^{-1.5}$')
        plt.plot(nu_ph_vals/1e9,1e26*(nu_ph_vals/1e9)**(-1),color='grey',ls='--',label=r'$\nu^{-1}$')
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
        plt.colorbar(sm,label='Time (yr)')
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        plt.xlim(1e-6,1e8)
        plt.ylim(1e10,1e30)
        plt.xlabel(r'$\nu$ (GHz)')
        plt.ylabel(r'$L_{\nu}$ (erg/s/Hz)')
        plt.savefig(f'Radio_curves_vs_frequency_{flarenum+1}.png',dpi=300,transparent=False,facecolor='white')

        np.savez(dNdgamma_dir+f'peak_data.npz',Lnu_pk=np.array(Lnu_pkvals),t_pk=np.array(t_pkvals),nu_pk=np.array(nu_pkvals))


    #### To plot SED at specific epochs ####
        fig=plt.figure()
        norm = LogNorm(vmin=epoch_list[0], vmax=epoch_list[-1])
        print("Plotting SED at epochs (days)", epoch_list)
        times_list = []
        count_list = []

        for epoch in epoch_list: #epochs are in days
            index_at_epoch=np.argmin(np.abs(tvals*m.tdyn_t0/secinday-epoch)/epoch )
            times_list.append(tvals[index_at_epoch]) #in units of tdyn_t0
            count_list.append(index_at_epoch)

        SEDs_to_save_dict = {}
        SEDs_to_save_dict['nu_ph_vals_GHz'] = nu_ph_vals/1e9 #in GHz

        for t,count,epoch in zip(times_list,count_list,epoch_list):
            print(count, t*m.tdyn_t0/secinday, 'days')
            Lnu_test,tau_ssa_test = emission_absorption_at_time(tvals[count],dNdgamma_vals[count],gamma_e_vals,delta_gamma_e,m.N_g,nu_ph_vals,
                                                                    m.rsh_ND_func,m.rsh_t0,B_ND_func,f_omega=f_omega)
            tau_ff_test = tau_ff(m.int_rhofwd_sq_dr_func(tvals[count]*m.tdyn_t0),nu_ph_vals,T_e_csm=T_e_csm)
            Lnu_spectrum = Lnu_test*np.exp(-tau_ff_test)*(1-np.exp(-tau_ssa_test))/tau_ssa_test
            SEDs_to_save_dict[f'epoch_{epoch}'] = Lnu_spectrum
            plt.plot(nu_ph_vals/1e9,Lnu_spectrum,color=plt.cm.viridis(norm(t*m.tdyn_t0/secinday)),label=f'{epoch:.1f} days'.format(epoch=epoch))
            count +=1
        plt.plot(nu_ph_vals/1e9,1e30*(nu_ph_vals/1e9)**(-1.5),color='black',ls='--',label=r'$\nu^{-1.5}$')
        plt.plot(nu_ph_vals/1e9,1e26*(nu_ph_vals/1e9)**(-1),color='grey',ls='--',label=r'$\nu^{-1}$')
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
        plt.colorbar(sm,label='Time (day)')
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        plt.xlim(1e-6,1e8)
        plt.ylim(1e10,1e30)
        plt.xlabel(r'$\nu$ (GHz)')
        plt.ylabel(r'$L_{\nu}$ (erg/s/Hz)')
        plt.savefig(f'SED_at_epochs_{flarenum+1}.png',dpi=300,transparent=False,facecolor='white')
        np.savez(dNdgamma_dir+f'SED_data.npz',SED_vs_epoch=SEDs_to_save_dict)
