#!/usr/bin/python2.7

import numpy as np#para cargar numpy
from pylab import *#cargar mtplotlib
import os#cargar herramientas de archivos
import re # use regular patterns
import sys, getopt # system commands
import string # string functions4
import math
#from scipy import linalg as scipylinalg
#from scipy import stats
#from scipy import ndimage
#from scipy import signal as scipysignal
import scipy
from pylab import *
from numpy import nan
import random
import subprocess
import pandas as pd
import glob
from astropy.io.votable import parse
from scipy.spatial import distance
from astropy.io import ascii

import timeit

start = timeit.default_timer()

catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_27-05-19_dupRem_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'


#df_highprob = pd.read_csv( 'highProbTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs_wNames' )
df_highprob_0 = pd.read_csv( catalog_name )
#df_highprob = df_highprob_0.loc[df_highprob_0['p_val'] >= 0.05]
#df_highprob = df_highprob_0.loc[df_highprob_0['P_B_A'] >= 0.05]
df_highprob = df_highprob_0.loc[df_highprob_0['P_A*P_B_A'] >= 0.0000001]

oc_list = df_highprob.Cluster_name.unique()
oc_list = np.sort(oc_list)
#print oc_list
print 'Number of clusters: ', len(oc_list), '\n'

Matriz = []

problematic = []

i = 1
for f in oc_list :
	
	print f , '\t\t', i, '/', len(oc_list)

	#if f != 'Pozzo_1' and f != 'NGC_6779' and f != 'NGC_6402' and f!= 'NGC_1746' and f!='NGC_6229' and f!= 'ASCC_63' and f!='NGC_6218' and f!= 'ASCC_60' and f!= 'Dolidze_34' and f!='NGC_1893' and f!='Ruprecht_164'and f!='Bochum_7'  and f!='NGC_6293' and f!='Trumpler_3' and f!='NGC_5272' and f!='Tombaugh_1' and f!='NGC_1904' and f!='Teutsch_48':  # problematic clusters
	#if f!= '001_Old' and f!= 'Gaia_3' and f!= 'DES_4' and f!= 'To_1' and f!= 'Palomar_15' and f!= 'DES_5' and f!= 'ESO-SC06' and f!='NGC_6229':
	#if f!= '001_Old' and f!= 'Gaia_3' and f!= 'DES_4' and f!= 'To_1' and f!= 'Palomar_15' and f!= 'DES_5' and f!= 'ESO-SC06' :
	if True and os.path.exists('Clusters/%s/%s_bestIsochrone_CG18.csv'%( f , f )):


		# open best isochrone
		dfiso_0 = pd.read_csv( 'Clusters/%s/%s_bestIsochrone_CG18.csv'%( f , f ) )	
		dfiso = dfiso_0.loc[dfiso_0['label'] == 1]

		G_iso = dfiso['G']
		BP_iso = dfiso['G_BP']
		RP_iso = dfiso['G_RP']

		Z = dfiso['Z'].iloc[0]
		logAge = dfiso['log_age_yr'].iloc[0]
		EBV = dfiso['bestFit_E_B_V'].iloc[0]
		dist = dfiso['bestFit_dist'].iloc[0]

		# check lowest magnitude with the label 1 (main sequence)
		G_turnoff = min(G_iso)

		# check gmag of the cepheid
		subdf_highprob = df_highprob.loc[df_highprob['Cluster_name'] == f]

		#if len(subdf_highprob_prev) == 1: # there is only 1 possible cepheid in that cluster

		G_cep  = subdf_highprob['phot_g_mean_mag'] - 0.85926*3.1*EBV + 5. - 5.*np.log10(dist)
		BP_cep = subdf_highprob['phot_bp_mean_mag'] - 1.06794*3.1*EBV + 5. - 5.*np.log10(dist)
		RP_cep = subdf_highprob['phot_rp_mean_mag'] - 0.65199*3.1*EBV + 5. - 5.*np.log10(dist)
		
		g_cep = subdf_highprob['phot_g_mean_mag'] 
		bp_cep = subdf_highprob['phot_bp_mean_mag'] 
		rp_cep = subdf_highprob['phot_rp_mean_mag'] 

		# if g_cep < g_turnoff, copy the file somewhere else  
		# and if len(cep) > 1

		if len(G_cep) > 1:
			print '\t\t Multiple cepheids form possible combos'
			j = 0
			while j < len(G_cep):

				if G_cep.iloc[j] <= G_turnoff:
					print '\t\t\t Combo   %s - %s    accepted'%(subdf_highprob['Cluster_name'].iloc[j], subdf_highprob['Var_Name'].iloc[j])
					#copy image
					os.system('cp ./Clusters_afterPreselection/%s_radecCMDwISO_CG18.png ./Clusters_afterPreselection/Filtered_by_CMD'%(f))

					#add combo name to list
					datata = subdf_highprob.iloc[j]
					if len(Matriz)<1:
						Matriz = datata
					else:
						Matriz = np.vstack((Matriz, datata))
				else:
					print '\t\t\t Combo   %s - %s    rejected'%(subdf_highprob['Cluster_name'].iloc[j], subdf_highprob['Var_Name'].iloc[j])


				j = j+1

		else:
			#print G_cep.iloc[0], G_turnoff
			print '\t\t Only one cepheid form a possible combos'

			if G_cep.iloc[0] <= G_turnoff:
				print '\t\t\t Combo   %s - %s    accepted'%(subdf_highprob['Cluster_name'].iloc[0], subdf_highprob['Var_Name'].iloc[0])

				os.system('cp ./Clusters_afterPreselection/%s_radecCMDwISO_CG18.png ./Clusters_afterPreselection/Filtered_by_CMD'%(f))

				datata = subdf_highprob.iloc[0]
				if len(Matriz)<1:
					Matriz = datata
				else:
					Matriz = np.vstack((Matriz, datata))

			else:
				print '\t\t\t Combo   %s - %s    rejected'%(subdf_highprob['Cluster_name'].iloc[0], subdf_highprob['Var_Name'].iloc[0])
			


	i = i+1
# ---------------------------------------------------------------------------------- #

# save list for future matching
df = pd.DataFrame( Matriz, columns=['Cluster' , 'Cluster_name' , 'name' , 'RA_deg' , 'DEC_deg' , 'Dist' , 'e_Dist_fromDist' , 'parallax_fromDist' , 'e_parallax_fromDist' , 'E(B-V)' , 'Age_D02' , 'logt_K13' , 'pmRA_1' , 'e_pmRA' , 'pmDE' , 'e_pmDE' , 'r1' , 'r2' , 'Nc' , 'RV' , 'e_RV' , 'o_RV' , 'FeH' , 'e_FeH' , 'o_FeH' , 'From' , 'pmRA_withCG' , 'e_pmRA_withCG' , 'pmDEC_withCG' , 'e_pmDEC_withCG' , 'parallax_cluster_withCG' , 'e_parallax_cluster_withCG' , 'RV_withCarrera19_CARRERA' , 'e_RV_withCarrera19_CARRERA' , 'o_RV_withCarrera19_CARRERA' , 'FeH_withCarrera19_CARRERA' , 'e_FeH_withCarrera19_CARRERA' , 'o_FeH_withCarrera19_CARRERA' , 'RV_FeH_from' , 'source_id_1' , 'designation' , 'ra' , 'dec' , 'l' , 'b' , 'ecl_lon' , 'ecl_lat' , 'parallax' , 'pmra' , 'pmdec' , 'phot_g_mean_mag' , 'phot_bp_mean_mag' , 'phot_rp_mean_mag' , 'ra_error' , 'dec_error' , 'parallax_error' , 'pmra_error' , 'pmdec_error' , 'a_g_percentile_lower' , 'a_g_percentile_upper' , 'a_g_val' , 'astrometric_chi2_al' , 'astrometric_excess_noise' , 'astrometric_excess_noise_sig' , 'astrometric_gof_al' , 'astrometric_matched_observations' , 'astrometric_n_bad_obs_al' , 'astrometric_n_good_obs_al' , 'astrometric_n_obs_ac' , 'astrometric_n_obs_al' , 'astrometric_params_solved' , 'astrometric_primary_flag' , 'astrometric_pseudo_colour' , 'astrometric_pseudo_colour_error' , 'astrometric_sigma5d_max' , 'astrometric_weight_al' , 'bp_g' , 'bp_rp' , 'dec_parallax_corr' , 'dec_pmdec_corr' , 'dec_pmra_corr' , 'duplicated_source' , 'e_bp_min_rp_percentile_lower' , 'e_bp_min_rp_percentile_upper' , 'e_bp_min_rp_val' , 'flame_flags' , 'frame_rotator_object_type' , 'g_rp' , 'lum_percentile_lower' , 'lum_percentile_upper' , 'lum_val' , 'matched_observations' , 'mean_varpi_factor_al' , 'parallax_over_error' , 'parallax_pmdec_corr' , 'parallax_pmra_corr' , 'phot_bp_mean_flux' , 'phot_bp_mean_flux_error' , 'phot_bp_mean_flux_over_error' , 'phot_bp_n_obs' , 'phot_bp_rp_excess_factor' , 'phot_g_mean_flux' , 'phot_g_mean_flux_error' , 'phot_g_mean_flux_over_error' , 'phot_g_n_obs' , 'phot_proc_mode' , 'phot_rp_mean_flux' , 'phot_rp_mean_flux_error' , 'phot_rp_mean_flux_over_error' , 'phot_rp_n_obs' , 'phot_variable_flag' , 'pmra_pmdec_corr' , 'priam_flags' , 'ra_dec_corr' , 'radial_velocity' , 'radial_velocity_error' , 'radius_percentile_lower' , 'radius_percentile_upper' , 'radius_val' , 'random_index' , 'ra_parallax_corr' , 'ra_pmdec_corr' , 'ra_pmra_corr' , 'ref_epoch' , 'rv_nb_transits' , 'rv_template_fe_h' , 'rv_template_logg' , 'rv_template_teff' , 'solution_id_1' , 'teff_percentile_lower' , 'teff_percentile_upper' , 'teff_val' , 'visibility_periods_used' , 'source_id_1a' , 'epoch_bp' , 'epoch_bp_error' , 'epoch_g' , 'epoch_g_error' , 'epoch_rp' , 'epoch_rp_error' , 'g_absorption' , 'g_absorption_error' , 'int_average_bp' , 'int_average_bp_error' , 'int_average_g' , 'int_average_g_error' , 'int_average_rp' , 'int_average_rp_error' , 'metallicity' , 'metallicity_error' , 'mode_best_classification' , 'multi_mode_best_classification' , 'num_clean_epochs_bp' , 'num_clean_epochs_g' , 'num_clean_epochs_rp' , 'p1_o' , 'p1_o_error' , 'p2_o' , 'p2_o_error' , 'p3_o' , 'p3_o_error' , 'peak_to_peak_bp' , 'peak_to_peak_bp_error' , 'peak_to_peak_g' , 'peak_to_peak_g_error' , 'peak_to_peak_rp' , 'peak_to_peak_rp_error' , 'pf' , 'pf_error' , 'phi21_g' , 'phi21_g_error' , 'phi31_g' , 'r21_g' , 'r21_g_error' , 'phi31_g_error' , 'r31_g' , 'r31_g_error' , 'solution_id_1a' , 'type2_best_sub_classification' , 'type_best_classification' , 'RAVE5_Name' , 'RAVE5_HRV' , 'RAVE5_e_HRV' , 'RV_gaiavsrave_GAIA' , 'e_RV_gaiavsrave_GAIA' , 'RV_gaiaravevsmelnik_MELNIK' , 'e_RV_gaiaravevsmelnik_MELNIK' , 'Var_Name' , 'Name_From' , 'VarType' , 'Separation_1' , 'Separation/3600. - 5*r1' , 'x' , 'P_A' , 'Vmag_melnik' , 'HRV_melnik' , 'e_HRV_melnik' , 'ruwe' , 'parallax_diff' , 'vr_diff' , 'pmRA_diff' , 'pmDEC_diff' , 'dfreed' , 'c' , 'p_val' , 'P_B_A' , 'P_A*p_val' , 'P_A*P_B_A'] )
df.to_csv( './Clusters_afterPreselection/afterCMD_combos.csv', index=False )

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #
