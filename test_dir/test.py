import sys
sys.path.append('/Users/samwu/codes/current_projects/RadioTDEFlares_module_branch/')
import numpy as np

import hydro_evol.flare_flare_collision as ffc
import radio_analysis.spectrum_evolution as spec_evol

ffc.evolve_flares(M_flares_list=[0.01,0.1],delta_t_list=[0.1,1],v_min_c=0.1,v_max_c=0.4)



