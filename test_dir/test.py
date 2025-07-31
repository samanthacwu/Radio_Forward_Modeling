import sys
sys.path.append('/Users/samwu/codes/current_projects/RadioTDEFlares_module_branch/')
import numpy as np

import hydro_evol.flare_flare_collision as ffc
import hydro_evol.flare_ism_collision as fic
import radio_analysis.spectrum_evolution as spec_evol
from radio_analysis.analysis import analyze_multiwavelength_spectrum, analyze_SED
from hydro_evol.model import Model
from hydro_evol.constants_list import *
# ffc.evolve_flares(M_flares=[0.01,0.1],delta_ts=[0.1,1],v_min_c=0.1,v_max_c=0.4)

# m = Model('./evolve_spectrum/shell_evolution_Mflare_1E-01_deltat_1.00yr/'+'shock_data.npz',simtype='flare_flare')
# print(m.times,m.times[0]/secinyear,m.times[-1]/secinyear)
# spec_evol.evolve_spectrum('flare_flare','.//evolve_spectrum/shell_evolution_Mflare_1E-01_deltat_1.00yr/',t0_in=0.35,tf_in=333,max_step_in=100000,print_int=10000)

analyze_multiwavelength_spectrum('flare_flare','./evolve_spectrum/shell_evolution_Mflare_1E-01_deltat_1.00yr/',freq_in=6,T_e_csm=1e4)

analyze_SED('flare_flare','./evolve_spectrum/shell_evolution_Mflare_1E-01_deltat_1.00yr/',freq_in=6,T_e_csm=1e4,
            epoch_list=[50, 100, 250, 500, 750, 1000, 1250, 1500],SED_interval=500)

# fic.evolve_flares(M_flare=0.01,rho_ism0s=[10,100],v_min_c=0.1,v_max_c=0.4,t0_in=0.01,p=0.5,p_ism=0,data_dir='./evolve_spectrum/')

# spec_evol.evolve_spectrum('flare_ism','./evolve_spectrum/shell_evolution_Mflare_1E-02_rhoISM0_1E+01yr/',t0_in=0.0101,tf_in=99.5,max_step_in=100000,print_int=10000)

# analyze_multiwavelength_spectrum('flare_ism','./evolve_spectrum/shell_evolution_Mflare_1E-02_rhoISM0_1E+01yr/',freq_in=6,T_e_csm=1e4)

# analyze_SED('flare_ism','./evolve_spectrum/shell_evolution_Mflare_1E-02_rhoISM0_1E+01yr/',freq_in=6,T_e_csm=1e4,
#             epoch_list=[50, 100, 250, 500, 750, 1000, 1250, 1500],SED_interval=500)