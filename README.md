# exomoon_project

Here is an explanation of the steps that were taken to generate the simulated lightcurves in the mixed-SNR runs-

GENERATING MODEL LIGHTCURVES

1) Using the code 'best_transit_fraction_systems.py', filter out the systems (out of 1000) for which each moon shows transits in at least half of the total transit epochs. 
2) Using the code 'most_detectable_systems.py' further find the systems with the top 10 highest median moon-transit-depth for both resonant and non-resonant systems and for each moon-number = 1,2,3,4,5. This gives total 100 systems and these were used to generate 100 simulated lightcurves. 
3) Using 'save_unpacked_systems_details.py', unpack and save system architecture details of the 100 systems. (Note: the values of parameters 'w_bary' and 'Omega_moon' were randonmly chosen in the 'best_transit_fraction_systems.py' code and the values are saved in 'best_transit_fraction_systems.txt'.)
4) Using 'save_model_LCs.py', generate and save the model lightcurves using Pandora for the 100 systems. 
