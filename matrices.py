# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 16:55:31 2024

@author: avalos-alais.s
"""
#FTRACT-DLPFC&IFG  --- Data preparation
import os
import numpy as np
import re
import tools.stats as ts
#Note : changed to use mask for N >= 50 & Nimp>=3. Same for p and delays (using N significant). Computed delays but not used to filter
#other version is to use get_reliable_atlas_entries_masks
# from tools.get_reliable_atlas_entries_masks import get_reliable_prob_mask as get_p_mask
# from tools.get_reliable_atlas_entries_masks import get_reliable_feat_mask as get_feat_mask
# from tools.get_reliable_atlas_entries_masks import get_N_mask
# from tools import marray
# from tools import index_tool as idx

def read_data(path_file, delays_flag = False):
    #path file should point to inner folder generated with save_collapse method. The folder name contains specific parameters the matrix were generated with. Hardcoded pass.
    #last parameter is the Zth it was generated with, which is part of the name of some feature's folder's
    zth = re.search(r'__?(zth\d+)$', path_file)
    if not zth: raise ValueError(f"Couldnt find zthX suffix in: {path_file}")
    zth = zth.group(1)
    # read probability, N_recordings, N_implantations, and Delays if needed
    p = np.loadtxt(os.path.join(path_file, 'probability.txt'))
    N = np.loadtxt(os.path.join(path_file, 'feature_ampl_' + zth , 'N_total.txt.gz'))
    N_imp = np.loadtxt(os.path.join(path_file, 'implantation_name', 'count_unique_str.txt.gz')) #?

    if delays_flag :
       peak_delay_med = np.loadtxt(os.path.join(path_file, 'feature_peak_delay_' + zth , 'nanquantile_0.5.txt.gz')) #Median Peak delay
       #Save quantiles 0.25 and 0.75
       pd_quant_25 = np.loadtxt(os.path.join(path_file, 'feature_peak_delay_' + zth , 'nanquantile_0.25.txt.gz'))
       pd_quant_75 = np.loadtxt(os.path.join(path_file, 'feature_peak_delay_' + zth, 'nanquantile_0.75.txt.gz'))
       peak_delay_mean = np.loadtxt(os.path.join(path_file, 'feature_peak_delay_' + zth , 'nanmean.txt.gz'))
       #onset_d = np.loadtxt('onset_delay__median.txt.gz') 
    
    assert np.all(np.isnan(p) == (N == 0))
    if not delays_flag : return p, N, N_imp 
    if delays_flag : return p, N, N_imp, peak_delay_med, peak_delay_mean, pd_quant_25,pd_quant_75
    
def save_array(filename, array):
    array = np.asarray(array)
    if array.ndim == 0:
        array = np.array([[array]])  # convertir escalar a 1x1
    np.savetxt(filename, array)
                
#def atlas_mat(path_folder, index_x, index_y,path_folder_output, name = 'xy' , flag_delays = False, min_n_suc = 50, min_n_fail = 50 , min_n_feat = 50, min_n_impl = 3 ,max_ci = 0.1) :# N_thresh = 50 ) : #min_n_suc = 50, min_n_fail = 50 , min_n_feat = 50, min_n_impl = 3 ,max_ci = 0.1) :
def atlas_mat(path_folder, index_x, index_y, path_folder_output, name='xy', flag_delays=False, N_thresh = 50 , N_imp_thresh = 3) : #min_n_suc = 50, min_n_fail = 50 , min_n_feat = 50, min_n_impl = 3 ,max_ci = 0.1) :

    if flag_delays :
        p, N, N_imp, peak_delay_med, peak_delay_mean, pd_quant_25,pd_quant_75  = read_data(path_folder, delays_flag= flag_delays)
    else : 
        p, N, N_imp = read_data(path_folder, delays_flag= flag_delays)

    shape = np.shape(p)
    if p.ndim == 2:
        nX, nY = shape
    elif p.ndim == 1:
        nX, nY = shape[0], 1  # vector
    else:  # escalar
        nX, nY = 1, 1
    #selection
    if index_x is not None and p.ndim >= 2:
        if isinstance(index_x, int):
            index_x = [index_x]
        p = p[index_x, :]
        N = N[index_x, :]
        N_imp = N_imp[index_x, :]

    if index_y is not None and p.ndim >= 2:
        if isinstance(index_y, int):
            index_y = [index_y]
        p = p[:, index_y]
        N = N[:, index_y]
        N_imp = N_imp[:, index_y]
    if flag_delays:
        #Selection
        if index_x is not None :
            peak_delay_med = peak_delay_med[index_x, :]
            peak_delay_mean= peak_delay_mean[index_x, :]
            pd_quant_25 = pd_quant_25[index_x, :]
            pd_quant_75 = pd_quant_75[index_x, :]
        if index_y is not None :
            peak_delay_med = peak_delay_med[:, index_y]
            peak_delay_mean= peak_delay_mean[:, index_y]
            pd_quant_25 = pd_quant_25[:, index_y]
            pd_quant_75 = pd_quant_75[:, index_y]
        #N for delays in merge and mask is Number of Stims above the threshold (N_with_value)
        N_EP = p*N
        # N_EP[np.isnan(p)] = 0
        #mask_d = get_N_mask(N_EP, N_thresh)  #
        #mask_d = get_feat_mask (N_EP, N_imp, min_n_feat, min_n_impl=min_n_impl)
        # mas N recordings
        mask_N = N_EP >= N_thresh
        # MASK ON N_implantations
        mask_imp = N_imp >= N_imp_thresh
        mask_d = mask_N & mask_imp

        peak_delay_med[~mask_d]=np.nan
        peak_delay_mean[~mask_d]=np.nan
        pd_quant_25[~mask_d]=np.nan
        pd_quant_75[~mask_d]=np.nan
        os.chdir(path_folder_output)
        save_array('peak_delay_mean_' + name + '.txt', peak_delay_mean)
        save_array('peak_delay_med_' + name + '.txt', peak_delay_med)
        save_array('pd_quantile_25_' + name + '.txt', pd_quant_25)
        save_array('pd_quantile_75_' + name + '.txt', pd_quant_75)
   
    # mask_p, ci = get_p_mask(p, N, N_imp, alpha=0.05, min_n_suc=min_n_suc, min_n_fail=min_n_fail, min_n_impl=min_n_impl, max_ci=max_ci, debug=False)

    #mas N recordings
    mask_N = N >= N_thresh
    #MASK ON N_implantations
    mask_imp = N_imp >= N_imp_thresh
    mask_p =  mask_N & mask_imp

    ci = ts.get_binom_ci(p, N, alpha = 0.05, ac_corr_flag=True)
    lost_elements = np.sum(~mask_p)
   #  print('Number of elements set to NAN, pmask (min/max n params):    ', lost_elements)
   #  total_N_entries = p.shape[0] * p.shape[1]
   #  print("Testing for min_suc", min_n_suc, "min_fail", min_n_fail, "min_feat", min_n_feat, "min_impl", min_n_impl, "max_ci", max_ci ,":")
   #  print("Mask p / N-entries : " , np.sum(mask_p) / total_N_entries)
    p[~mask_p]=np.nan
    # save merged matrices
    os.chdir(path_folder_output) #ill save just corrected mats
    save_array('p_' + name + '.txt', p)
    save_array('N_' + name + '.txt', N)
    save_array('I_' + name + '.txt', N_imp)
    save_array('CI_' + name + '.txt', ci)
    #Save indexes (relative to original parcel list, to find names)
    with open(f'index_x_{name}.txt', 'w') as f:
        if index_x is not None:
            for i in index_x:
                f.write(f"{int(i)}\n")
        else:
            for i in range(nX):
                f.write(f"{int(i)}\n")

    with open(f'index_y_{name}.txt', 'w') as f:
        if index_y is not None:
            for i in index_y:
                f.write(f"{int(i)}\n")
        else:
            for i in range(nY):
                f.write(f"{int(i)}\n")
    #, pd_quant_25, pd_quant_75
    return (p, N, N_imp, peak_delay_med, peak_delay_mean) if flag_delays else (p, N, N_imp) #N and I returned correspond to p. Number of stims not num of EP signif


