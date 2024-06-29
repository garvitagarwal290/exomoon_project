# exomoon_project

Paper: https://academic.oup.com/mnras/article/529/2/1232/7617719

Here is an explanation of the steps that were taken to generate the simulated lightcurves in the mixed-SNR runs-


## GENERATING MODEL LIGHTCURVES

1) Using the code 'best_transit_fraction_systems.py', filter out the systems (out of 1000) for which each moon shows transits in at least half of the total transit epochs. Output file: 'best_transit_fraction_systems_{nores/res}.txt'
2) Using the code 'most_detectable_systems.py' further find the systems with the top 10 highest median moon-transit-depth for both resonant and non-resonant systems and for each moon-number = 1,2,3,4,5. This gives total 100 systems and these were used to generate 100 simulated lightcurves. Output file: 'most_detectable_systems_final.txt'
3) Using 'save_unpacked_systems_details.py', unpack and save system architecture details of the 100 systems. (Note: the values of parameters 'w_bary' and 'Omega_moon' were randonmly chosen in the 'best_transit_fraction_systems.py' code and the values are saved in 'best_transit_fraction_systems.txt'.)
4) Using 'save_model_LCs.py', generate and save the model lightcurves using Pandora for the 100 systems. 


## SELECTING GOOD KEPLER LIGHTCURVES

1) Given list of Kepler IDs of lightcurves without any planet, using 'dw_statistic_filter.py' filter out the lightcurves for which the Durbin-Watson statistic after detrending is between 1.5 and 2.5. Output file: 'good_kepler_data.txt'
2) Using 'sort_kepler_bynoise.py', sort the lightcurves in 'good_kepler_data.txt' in the increasing order of flux uncertainties.


## PREPARING FOR INJECTION 

The code 'pre-injection.py' decides which system is to be injected into which Kepler lightcurve. For every pair of system-kepler_lightcurve it makes sure that the 4 data quality requirements are met: (i) SNR of each moon is > 2.0 (ii) actual number of transits is not < 3 (iii) transit-fraction for each moon is not < 0.5 (iv) DW statistic is between 1.5 and 2.5. The code summarises the details of these 100 pairs in the file 'catalogue.txt'

## INJECTION

The code 'injection.py' finally injects the model LC of each of the 100 systems in the Kepler LC as mentioned in the catalogue.txt file. The 100 final lightcurves are saved.

## MODELING

The code 'modeling.py' runs the model fitting on the simulated lightcurves using the nested sampler Ultranest. 


Following steps were taken to generate the simulated lightcurves for the fixed-SNR runs-


## SAVING THE LIGHTCURVE TRENDS

Since I didnt save the trends of the lightcurves when I made the mixed-SNR lightcurves, I had to do that separately later using the script 'save_trends.py'. The script detrends the 100 kepler lightcurves using cofiam and saves them in files.

## GENERATING THE FIXED SNR LIGHTCURVES

The script 'modify_SNR.py' first modifies the noise in the kepler lightcurves using the trends saved earlier and the target fixed SNR value. Then it performs the injection of the model lightcurve on to the modified Kepler lightcurve exactly like 'injeciton.py'. The new 100 lightcurves are then ready to be model-fitted by an appropriately modified 'modeling.py' script. 

