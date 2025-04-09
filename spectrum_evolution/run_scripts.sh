
# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_deltat_1yr.txt' --data_dir './evolve_spectrum_epsB_1e-2/'

python spectrum_evolution.py './evolve_spectrum_epsB_1e-2/shell_evolution_deltat_1yr/' './evolve_spectrum_epsB_1e-2/shell_evolution_deltat_1yr/' --p_exp=2.1 --eps_B=0.01 --eps_E=0.01 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.15 --tf=135

# python analyze_multiwavelength_spectrum.py --p_exp=2.1 --f_omega=1.00 --eps_E=0.001 --eps_B=0.01 --model_dir './evolve_shock/shell_evolution_deltat_1yr/' --evolve_shock_dir './evolve_shock/shell_evolution_deltat_1yr/' --data_dir './evolve_spectrum_epsB_1e-2/shell_evolution_deltat_1yr/' --dNdgamma_dir './evolve_spectrum_epsB_1e-2/shell_evolution_deltat_1yr/lower_epsE/'
# mv Radio_curve*.png ./evolve_spectrum_epsB_1e-2/shell_evolution_deltat_1yr/lower_epsE
