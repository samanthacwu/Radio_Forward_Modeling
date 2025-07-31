import sys
sys.path.append('/Users/samwu/codes/current_projects/RadioTDEFlares_module_branch/')
import numpy as np

import hydro_evol.flare_flare_collision as ffc
import hydro_evol.flare_ism_collision as fic
import radio_analysis.spectrum_evolution as spec_evol

# ffc.evolve_flares(M_flares=[0.01,0.1],delta_ts=[0.1,1],v_min_c=0.1,v_max_c=0.4)

# fic.evolve_flares(M_flare=0.01,rho_ism0s=[10,100],v_min_c=0.1,v_max_c=0.4,t0_in=0.01,p=0.5,p_ism=0,data_dir='./evolve_spectrum/')

spec_evol.evolve_spectrum('flare_ism','./evolve_spectrum/shell_evolution_Mflare_1E-02_rhoISM0_1E+01yr/',t0_in=0.0101,max_step_in=2,print_int=1)
