import sys
sys.path.append('/Users/samwu/codes/current_projects/Radio_Forward_Modeling/')
import numpy as np
import os
import hydro_evol.flare_flare_collision_adiabatic as ffc
import radio_analysis.spectrum_evolution as spec_evol
from radio_analysis.analysis import analyze_multiwavelength_spectrum, analyze_SED
from hydro_evol.model import Model
from hydro_evol.constants_list import *


#### run this line to get the file shock_data.npz in each directory, which contains the properties of the shock vs. time
#### remember to create the directory input to data_dir if it doesn't exist yet.

# ffc.evolve_flares(M_flares=[0.01],delta_ts=[0.1,0.3,1,2],v_min_c=0.04,v_max_c=0.4,data_dir='./evolve_spectrum_adiabatic/',dt_scale=1e-3)
#implicitly assuming eps_B=eps_E=0.1, electron spectrum power law p=3

#### can run again for larger flare mass if necessary
# ffc.evolve_flares(M_flares=[0.1],delta_ts=[0.1,0.3,1,2],v_min_c=0.04,v_max_c=0.4)

#### this will loop through all the data directories created under data_dir. 
top_dir = './evolve_spectrum_adiabatic/'
for i,dirname in enumerate(os.listdir(top_dir)):
    print(dirname)

    #### set name of the path to the data directory that is being loaded during this step.
    pathname = top_dir+dirname + '/'
    if not os.path.isdir(pathname):
        continue
    #### the next line sets up the Model class for the radio forward modeling
    m = Model(pathname+'shock_data.npz',simtype='flare_flare')
    print(m.times[0]/secinyear,m.times[-1]/secinyear) #### this tells you the initial and final time of the shock evolution in years.

    #### set the evolution of the radio spectrum to the same time interval as the shock evolution, 
    #### with slightly buffered start and end times to prevent interpolation errors.
    t0_in = 1.01*m.times[0]/secinyear
    tf_in = 0.99*m.times[-1]/secinyear

    # #### the next line calculates the electron spectrum at each time
    # spec_evol.evolve_spectrum(m,'flare_flare',pathname,t0_in=t0_in,tf_in=tf_in,max_step_in=100000,print_int=10000,plotting=False)

    # #### the above line created the electron spectrum (as a function of frequency), one spectrum per timestep. 
    # #### the next line generates the radio light curve at a given frequency, freq_in
    # analyze_multiwavelength_spectrum(m,'flare_flare',pathname,freq_in=6,T_e_csm=1e4)

    ### you can also generate the SEDs at different epochs, given in epoch_list.

    analyze_SED(m,'flare_flare',pathname,T_e_csm=1e4,
                epoch_list=list(np.arange(25,2000,250)),SED_interval=1000)


