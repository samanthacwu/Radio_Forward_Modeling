import sys
#put your path to where the Radio_Forward_Modeling directory is located 
sys.path.append('/Users/samwu/codes/current_projects/Radio_Forward_Modeling/')
import numpy as np

import hydro_evol.CSM_density_profile as csm_density
import hydro_evol.SNejecta_CSM_hydro as hydro_evol
import radio_analysis.spectrum_evolution as spec_evol
from radio_analysis.analysis import analyze_multiwavelength_spectrum, analyze_SED
from hydro_evol.model import Model
from hydro_evol.constants_list import *


#### replace './3.68M_Porb1d/' with the directory where your MESA model is located, if you are creating a CSM density profile from a MESA model, and uncomment the next line.

#csm_density.CSM_density_profile('./3.68M_Porb1d/','./evolve_shock/',1e51,vel_factor=0.3,f_omega=1,wind_scaling_factor=1)

#### if you did the above step, it creates a directory with the below path. 
#### if you are starting with your own input density profile, save it as an npz file called density_prof.npz (or whatever you prefer) in a directory and put the path to that directory below for path_to_output
#### the contents of density_prof.npz are just the density and radius arrays, e.g.
#### np.savez(save_density_dir+'density_prof.npz',rho=density,r=dist_arr)

path_to_output = './evolve_shock/E_1.00E+51_Mej_2.10E+00_rhoprof_3.680M_Porb1d_velfac_0.30_fomega_1.00/'

#### run this line to get the file shock_data.npz, which contains the properties of the shock vs. time
#### P.S.: if this is ever slow, go into the hydro_evol file and play with the initial conditions

hydro_evol.evolve_ejectaCSM_shock(path_to_output+'density_prof.npz',path_to_output,E=1e51,Mej=2.10,
                               f_omega=1,R0_in=0.11,dt_in=1e-2,tf_in=1e2,t0_in=1.6e-5,max_step=50000,print_int=10000)

##after you have created shock_data.npz, run the next line to set up the Model class for the radio forward modeling

m = Model(path_to_output+'shock_data.npz',simtype='SNejecta_CSM',X_H=0,X_He=1,calculate_SBO=True,integrated_Bfield=True) #here can specify eps_B, eps_E, f_omega.
print(m.times,m.times[0]/secinyear,m.times[-1]/secinyear)

##the next line calculates the electron spectrum at each time time

spec_evol.evolve_spectrum(m,'SNejecta_CSM',path_to_output,t0_in=1.01*m.times[0]/secinyear,tf_in=0.99*m.times[-1]/secinyear,max_step_in=100000,print_int=10000,plotting=False)

##once you have the electron spectrum, you can generate the light curve at input frequency freq_in

analyze_multiwavelength_spectrum(m,'SNejecta_CSM',path_to_output,freq_in=3,T_e_csm=1e5)

##and you can generate the SEDs at different epochs, given in epoch_list

analyze_SED(m,'SNejecta_CSM',path_to_output,T_e_csm=1e5,
            epoch_list=[10,50],SED_interval=500)
