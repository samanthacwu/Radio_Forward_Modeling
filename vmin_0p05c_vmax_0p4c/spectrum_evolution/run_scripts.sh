
# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_deltat_0.3yr.txt' --data_dir './evolve_spectrum_epsB_1e-1/'

# python spectrum_evolution.py './evolve_spectrum_epsB_1e-1/shell_evolution_deltat_0.3yr/' './evolve_spectrum_epsB_1e-1/shell_evolution_deltat_0.3yr/' --p_exp=3 --eps_B=0.1 --eps_E=0.1 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.05 --tf=41.9

# python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_deltat_0.3yr/' --evolve_shock_dir './evolve_shock/shell_evolution_deltat_0.3yr/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_deltat_0.3yr/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_deltat_0.3yr/'
# mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_deltat_0.3yr/

python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 50 100 250 500 750 1000 1250 1500 --model_dir './evolve_shock/shell_evolution_deltat_0.3yr/' --evolve_shock_dir './evolve_shock/shell_evolution_deltat_0.3yr/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_deltat_0.3yr/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_deltat_0.3yr/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_deltat_0.3yr/