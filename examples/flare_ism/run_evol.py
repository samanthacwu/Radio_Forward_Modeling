import sys
sys.path.append('/Users/samwu/codes/current_projects/Radio_Forward_Modeling/')
import numpy as np
import os
import hydro_evol.flare_ism_collision as fic
import radio_analysis.spectrum_evolution as spec_evol
from radio_analysis.analysis import analyze_multiwavelength_spectrum, analyze_SED
from hydro_evol.model import Model
from hydro_evol.constants_list import *

def r_infl(MBH): #MBH in Msun units
    
    return 13 * pc_cm *(MBH/1e8)**(1./2.)

r_out = r_infl(1e7) # r_out in pc_cm. currently using the sphere of influence of the black hole to set the outer radius for ff absorption
print('using r_out =', r_out/pc_cm,'pc')

#implicitly assuming eps_B=eps_E=0.1, electron spectrum power law p=3
fic.evolve_flares(M_flare=0.01,rho_ism0s=[1e0,1e1,1e2,1e3],v_min_c=0.04,v_max_c=0.4,t0_in=0.001,stop_ratio=110000,p=0.5,p_ism=0,r_out=r_out,data_dir='./evolve_spectrum/')
# fic.evolve_flares(M_flare=0.01,rho_ism0s=[1e0,3e0,1e1,1e2,1e3,1e4],v_min_c=0.04,v_max_c=0.4,t0_in=0.001,stop_ratio=110000,p=0.5,p_ism=0,r_out=r_out,data_dir='./evolve_spectrum/')

# for i,dirname in enumerate(os.listdir('./evolve_spectrum/')):
# # for i,dirname in enumerate(os.listdir('./evolve_spectrum_pflare1/')):
    
# #     # if dirname.split('_')[-2][0]=='1' and dirname.split('_')[3]!='7E-02':
# #     # if dirname.split('_')[3]!='7E-03':
# #     #     continue
# #     # print(dirname)
#     # if dirname != 'shell_evolution_Mflare_1E-02_pISM_0_rhoISM0_3E+01m_H' and dirname != 'shell_evolution_Mflare_1E-02_pISM_0_rhoISM0_3E+02m_H' and dirname != 'shell_evolution_Mflare_1E-01_pISM_0_rhoISM0_3E+02m_H':
#     #     continue
#     # if dirname != 'shell_evolution_Mflare_1E-02_pISM_0_rhoISM0_1E+02m_H':
#     #     continue
#     pathname = './evolve_spectrum/'+dirname + '/'
#     print(dirname)
#     # pathname = './evolve_spectrum_pflare1/'+dirname + '/'
    
#     m = Model(pathname+'shock_data.npz',simtype='flare_ism')
#     print(m.times[0]/secinyear,m.times[-1]/secinyear)
#     # t0_in = 1.01*m.times[0]/secinyear
#     # tf_in = 0.99*m.times[-1]/secinyear
#     # spec_evol.evolve_spectrum('flare_ism',pathname,t0_in=t0_in,tf_in=tf_in,max_step_in=100000,print_int=10000,plotting=False)

#     # analyze_multiwavelength_spectrum(m,'flare_ism',pathname,freq_in=15.5,T_e_csm=1e4)

#     analyze_SED(m,'flare_ism',pathname,freq_in=6,T_e_csm=1e4,
#                 epoch_list=list(np.arange(25,2500,25)),SED_interval=500)

