
python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_deltat_1yr.txt' --data_dir './evolve_spectrum_epsB_1e-2/'

# python spectrum_evolution.py './evolve_spectrum_p2p1_epsB_1e-2/E_1.00E+50_Mej_1.50E+00_fomega_0.75/' './evolve_spectrum_p2p1_epsB_1e-2/E_1.00E+50_Mej_1.50E+00_fomega_0.75/lower_epsE/' --p_exp=2.1 --eps_B=0.01 --eps_E=0.001 --f_omega=1.00 --max_step=20000 --print_int=1000 --dt_sc=1e-3 --t0=1 --tf=1263

# python analyze_multiwavelength_spectrum.py --p_exp=2.1 --f_omega=1.00 --eps_E=0.001 --eps_B=0.01 --model_dir './evolve_shock/E_1.00E+50_Mej_1.50E+00_fomega_0.75/' --evolve_shock_dir './evolve_shock/E_1.00E+50_Mej_1.50E+00_fomega_0.75/' --data_dir './evolve_spectrum_p2p1_epsB_1e-2/E_1.00E+50_Mej_1.50E+00_fomega_0.75/' --dNdgamma_dir './evolve_spectrum_p2p1_epsB_1e-2/E_1.00E+50_Mej_1.50E+00_fomega_0.75/lower_epsE/'
# mv Radio_curve*.png ./evolve_spectrum_p2p1_epsB_1e-2/E_1.00E+50_Mej_1.50E+00_fomega_0.75/lower_epsE
