import sys
sys.path.append('/Users/samwu/codes/current_projects/Radio_Forward_Modeling/')
import numpy as np
import os
import hydro_evol.flare_flare_collision_adiabatic as ffc
import radio_analysis.spectrum_evolution as spec_evol
from radio_analysis.analysis import analyze_multiwavelength_spectrum, analyze_SED
from hydro_evol.model import Model
from hydro_evol.constants_list import *


#implicitly assuming eps_B=eps_E=0.1, electron spectrum power law p=3
ffc.evolve_flares(M_flares=[0.01],delta_ts=[0.1,0.3,1,2],v_min_c=0.04,v_max_c=0.4,data_dir='./evolve_spectrum_adiabatic/',dt_scale=1e-3)
# ffc.evolve_flares(M_flares=[0.1],delta_ts=[0.1,0.3,1,2],v_min_c=0.04,v_max_c=0.4)

# ffc.evolve_flares(M_flares=[0.1],delta_ts=[0.3],v_min_c=0.04,v_max_c=0.4,p=1,data_dir='./evolve_spectrum_pflare1/')
# for i,dirname in enumerate(os.listdir('./evolve_spectrum/')):
# # for i,dirname in enumerate(os.listdir('./evolve_spectrum_pflare1/')):
    
#     # if dirname.split('_')[-1]!='0.20yr' and dirname.split('_')[-1]!='0.10yr':
#     #     continue
#     # if dirname.split('_')[-1]!='1.50yr' and dirname.split('_')[-1]!='0.05yr' and dirname.split('_')[-1]!='0.60yr' and dirname.split('_')[-1]!='0.20yr':
#     #     continue
#     # if dirname != 'shell_evolution_Mflare_1E-01_deltat_0.05yr':
#     #     continue
#     # print(dirname)
#     pathname = './evolve_spectrum/'+dirname + '/'
# #     pathname = './evolve_spectrum_pflare1/'+dirname + '/'
    
#     m = Model(pathname+'shock_data.npz',simtype='flare_flare')
#     print(m.times[0]/secinyear,m.times[-1]/secinyear)
#     t0_in = 1.01*m.times[0]/secinyear
#     tf_in = 0.99*m.times[-1]/secinyear
#     # spec_evol.evolve_spectrum(m,'flare_flare',pathname,t0_in=t0_in,tf_in=tf_in,max_step_in=100000,print_int=10000,plotting=False)

#     # analyze_multiwavelength_spectrum(m,'flare_flare',pathname,freq_in=15,T_e_csm=1e4)

#     analyze_SED(m,'flare_flare',pathname,freq_in=6,T_e_csm=1e4,
#                 epoch_list=list(np.arange(25,2000,25)),SED_interval=500)


