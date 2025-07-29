
# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_rhoISM_250e-24_cgs.txt' --data_dir './evolve_spectrum_epsB_1e-1/'
# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_rhoISM_500e-24_cgs.txt' --data_dir './evolve_spectrum_epsB_1e-1/'
# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_rhoISM_750e-24_cgs.txt' --data_dir './evolve_spectrum_epsB_1e-1/'

python spectrum_evolution.py './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_250e-24_cgs/' './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_250e-24_cgs/' --p_exp=3 --eps_B=0.1 --eps_E=0.1 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.0101 --tf=99.5
python spectrum_evolution.py './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_500e-24_cgs/' './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_500e-24_cgs/' --p_exp=3 --eps_B=0.1 --eps_E=0.1 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.0101 --tf=99.5
python spectrum_evolution.py './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_750e-24_cgs/' './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_750e-24_cgs/' --p_exp=3 --eps_B=0.1 --eps_E=0.1 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.0101 --tf=99.5

python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_250e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_250e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_250e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_250e-24_cgs/'
mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_250e-24_cgs/
python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_500e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_500e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_500e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_500e-24_cgs/'
mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_500e-24_cgs/
python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_750e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_750e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_750e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_750e-24_cgs/'
mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_750e-24_cgs/

python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 50 100 250 300 400 450 500 750 1000 1250 1500 2000 --model_dir './evolve_shock/shell_evolution_rhoISM_250e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_250e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_250e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_250e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_250e-24_cgs/
python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 50 100 250 300 400 450 500 750 1000 1250 1500 2000 --model_dir './evolve_shock/shell_evolution_rhoISM_500e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_500e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_500e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_500e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_500e-24_cgs/
python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 50 100 250 500 750 1000 1250 1500 2000 --model_dir './evolve_shock/shell_evolution_rhoISM_750e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_750e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_750e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_750e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_750e-24_cgs/

# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_rhoISM_1e-24_cgs.txt' --data_dir './evolve_spectrum_epsB_1e-1/'
# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_rhoISM_10e-24_cgs.txt' --data_dir './evolve_spectrum_epsB_1e-1/'
# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_rhoISM_100e-24_cgs.txt' --data_dir './evolve_spectrum_epsB_1e-1/'
# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_rhoISM_1000e-24_cgs.txt' --data_dir './evolve_spectrum_epsB_1e-1/'
# python save_for_spectrum_evolution.py  --evolve_shock_dir '../shell_propagation/' --filename 'shell_evolution_rhoISM_10000e-24_cgs.txt' --data_dir './evolve_spectrum_epsB_1e-1/'

python spectrum_evolution.py './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10000e-24_cgs/' './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10000e-24_cgs/' --p_exp=3 --eps_B=0.1 --eps_E=0.1 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.0101 --tf=99.5
python spectrum_evolution.py './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/' './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/' --p_exp=3 --eps_B=0.1 --eps_E=0.1 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.0101 --tf=99.5
python spectrum_evolution.py './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/' './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/' --p_exp=3 --eps_B=0.1 --eps_E=0.1 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.0101 --tf=99.5
python spectrum_evolution.py './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/' './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/' --p_exp=3 --eps_B=0.1 --eps_E=0.1 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.0101 --tf=99.5
python spectrum_evolution.py './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1e-24_cgs/' './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1e-24_cgs/' --p_exp=3 --eps_B=0.1 --eps_E=0.1 --f_omega=1.00 --max_step=100000 --print_int=10000 --dt_sc=1e-3 --t0=0.0101 --tf=99.5

python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_10000e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_10000e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10000e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10000e-24_cgs/'
mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10000e-24_cgs/

python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_1000e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_1000e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/'
mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/

python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_100e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_100e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/'
mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/

python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_10e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_10e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/'
mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/

python analyze_multiwavelength_spectrum.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --model_dir './evolve_shock/shell_evolution_rhoISM_1e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_1e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1e-24_cgs/'
mv Radio_curve*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1e-24_cgs/

python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 50 100 250 500 750 1000 1250 1500 --model_dir './evolve_shock/shell_evolution_rhoISM_10000e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_10000e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10000e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10000e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10000e-24_cgs/

python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 50 100 250 500 750 1000 1250 1500 --model_dir './evolve_shock/shell_evolution_rhoISM_1000e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_1000e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1000e-24_cgs/

python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 50 100 250 500 750 1000 1250 1500 2000 --model_dir './evolve_shock/shell_evolution_rhoISM_100e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_100e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_100e-24_cgs/

python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 50 100 250 500 750 1000 1250 1500 --model_dir './evolve_shock/shell_evolution_rhoISM_10e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_10e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_10e-24_cgs/

python analyze_SED.py --p_exp=3 --f_omega=1.00 --eps_E=0.1 --eps_B=0.1 --Te_csm=1e4 --epochs 50 100 250 500 750 1000 1250 1500 --model_dir './evolve_shock/shell_evolution_rhoISM_1e-24_cgs/' --evolve_shock_dir './evolve_shock/shell_evolution_rhoISM_1e-24_cgs/' --data_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1e-24_cgs/' --dNdgamma_dir './evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1e-24_cgs/'
mv SED*.png ./evolve_spectrum_epsB_1e-1/shell_evolution_rhoISM_1e-24_cgs/
