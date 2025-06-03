
# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_rhoISM_1000e-24_cgs.txt' --data_dir './evolve_spectrum_epsB_1e-1/'

# python spectrum_evolution.py './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/' './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/' --p_exp=3 --eps_B=0.1 --eps_E=0.1 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.101 --tf=99.5

# python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_1000e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_1000e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/'
# mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/

# python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_100e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_100e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/'
# mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/

# python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_10e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_10e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/'
# mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/

python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 526 769 945 1092 1103 1325 1343 --model_dir './evolve_shock/shell_evolution_rhoISM_10e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_10e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/

python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 526 769 945 1092 1103 1325 1343 --model_dir './evolve_shock/shell_evolution_rhoISM_100e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_100e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/

python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 526 769 945 1092 1103 1325 1343 --model_dir './evolve_shock/shell_evolution_rhoISM_1000e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_1000e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/