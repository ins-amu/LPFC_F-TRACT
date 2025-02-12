# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 16:55:31 2024

@author: avalos-alais.s
"""

#FTRACT-DLPFC&IFG  --- Data preparation 
import os
import numpy as np
from tools.get_reliable_atlas_entries_masks import get_reliable_prob_mask as get_p_mask
from tools.get_reliable_atlas_entries_masks import get_reliable_feat_mask as get_feat_mask
from tools import marray
from tools import index_tool as idx

def read_data(path_file, delays_flag = False):
    
    atlas_files = np.load(path_file )
    #4th line is the Zth = 5. Discard of unknown
    p = atlas_files['probabilities_early_above_th'][4, 1:,1:] 
    N = atlas_files['N_stims'][1:,1:]                      
    N_imp = atlas_files['N_crfs'][1:,1:]                   
    
    if delays_flag : 
       peak_delay_med = atlas_files['feature__peak_delay__median'][4, 1:,1:]   
       peak_delay_mean = atlas_files['feature__peak_delay__mean'][4, 1:,1:]  
       
       #onset_d = np.loadtxt('onset_delay__median.txt.gz') 
    
    assert np.all(np.isnan(p) == (N == 0)) 
    
    if not delays_flag : return p, N, N_imp 
    if delays_flag : return p, N, N_imp, peak_delay_med, peak_delay_mean 
    

                
def atlas_mat(path_folder, index_x , index_y , path_folder_output, name = 'xy' , flag_delays = False, Eff_flag = True) : 
    # index in shape for merge_in_array : [tuple]. If not merge needed [(idx1,), (indx2,)]
    if flag_delays : 
        p, N, N_imp, peak_delay_med, peak_delay_mean  = read_data(path_folder, delays_flag= flag_delays)
    else : 
        p, N, N_imp = read_data(path_folder, delays_flag= flag_delays)
        
    if Eff_flag : 
        rows_merge = index_x
        cols_merge = index_y
    else :
        rows_merge = index_y
        cols_merge = index_x

    #For no merging pass index in way just to select matrix   
    p_xy, Np_xy, Ip_xy = marray.merge_in_array(p, N, rows_merge, cols_merge, N_imp) 

    
    if flag_delays: 
        #N for delays in merge and mask is Number of Stims above the threshold
        N_EP = p*N
        N_EP[np.isnan(p)] = 0
        D_med_xy, N_xy, I_xy = marray.merge_in_array(peak_delay_med, N_EP, rows_merge, cols_merge, N_imp, debug = False)
        mask_d = get_feat_mask (N_xy, I_xy, 50)
        D_med_xy[~mask_d]=np.nan 
        Np_mergedp = Np_xy * p_xy
        Np_mergedp[np.isnan(p_xy)] = 0
        control = N_xy - Np_mergedp 
        assert np.all(np.abs(control) <= 1/1000)
       
        D_mean_xy, N_xy, I_xy = marray.merge_in_array(peak_delay_mean, N_EP, rows_merge, cols_merge, N_imp, debug = False)
        mask_d = get_feat_mask (N_xy, I_xy, 50)
        D_mean_xy[~mask_d]=np.nan 
        os.chdir(path_folder_output) 
        # print('saving : ' + 'peak_delay_med_' + name + '.txt')
        # np.savetxt('peak_delay_med_' + name + '.txt' , D_med_xy)
        # print('saving : ' + 'peak_delay_mean_' + name + '.txt')
        np.savetxt('peak_delay_mean_' + name + '.txt', D_mean_xy)
        
   
    mask_p, ci = get_p_mask(p_xy, Np_xy, Ip_xy, alpha=0.05, min_n_suc=50, min_n_fail=50, min_n_impl=3, max_ci=0.1, debug=False)
    lost_elements = np.sum(~mask_p)
    print('Number of elements set to NAN, pmask (min/max n params):    ', lost_elements)
    p_xy[~mask_p]=np.nan 
    # save merged matrices
    os.chdir(path_folder_output) #ill save just corrected mats 
    #save indexes related to those matrices 
    if index_x is not None : idx.save_idx_txt( index_x, 'index_x_' + name + '.txt')
    if index_y is not None : idx.save_idx_txt( index_y, 'index_y_'+ name + '.txt')
    np.savetxt('p_' + name +'.txt', p_xy)
    np.savetxt('N_' + name + '.txt', Np_xy)
    # np.savetxt('I_' + name + '.txt' , Ip_xy)
    np.savetxt('CI_' + name + '.txt' , ci)
   
    if flag_delays: return p_xy, Np_xy, Ip_xy, D_med_xy, D_mean_xy #N and I returned correspond to p. Number of stims not num of EP signif
    else: return p_xy, Np_xy, Ip_xy

