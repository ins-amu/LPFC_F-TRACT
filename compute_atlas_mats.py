#Statistic correction and shaping of atlas matrices.
#New version of 'compute_and_plot' independent of plotting.
import os
from numpy.core.defchararray import startswith
import matrices
import numpy as np
import Definitions
from parcellation2.parcellation2 import Parcellation2

# %% DEFINITIONS -  X rows ; Y columns of the matrix
L_s = Definitions.L_s
res33 = Definitions.res33
res125 = Definitions.res125
res250 = Definitions.res250  # FN analysis
res500 = Definitions.res500
# Parcellation labels
labels_33 = Definitions.labels_33
labels_125 = Definitions.labels_125
labels_250 = Definitions.labels_250
labels_500 = Definitions.labels_500
# For plotting
labels_path = Definitions.labels_path
labels_path_33 = Definitions.labels_path_33
labels_path_125 = Definitions.labels_path_125
labels_path_250 = Definitions.labels_path_250
labels_path_500 = Definitions.labels_path_500

# Selection of labels in use
# Labels 33 corrected in Definitions. If problems with no labels found or only subcortical you might be using a version with 'ctx-lh-' ...
labels_rec_33 = [name for name in labels_33 if name.startswith("lh.") or name.startswith("rh.") or name == 'Left-Hippocampus' or name == 'Left-Amygdala' or name == 'Right-Hippocampus' or name == 'Right-Amygdala']
p = Parcellation2(Definitions.parcellation_path, res33)
index_33 = [p.get_index_by_name(i) for i in labels_33]
index_rec_33 = np.transpose([p.get_index_by_name(i) for i in labels_rec_33])
# index all (other parcellations with specific selections or none)
p = Parcellation2(Definitions.parcellation_path, res125)
index_125 = [p.get_index_by_name(i) for i in labels_125]
p = Parcellation2(Definitions.parcellation_path, res250)
index_250 = [p.get_index_by_name(i) for i in labels_250]
p = Parcellation2(Definitions.parcellation_path, res500)
index_500 = [p.get_index_by_name(i) for i in labels_500]

def compute_mats (base_path, output_folder, collapse_conf_name ,N_thresh, N_imp_thresh, gral_tw = [0,100], dir_tw = [0,50], ind_tw = [100,400], min_cb_dir = 20, max_cb_dir= 35, min_cb_ind = 150, max_cb_ind = 200, max_vb_gral = 80, max_avg = 0.175, min_avg = 0.1, max_cb_p = 0.25, max_p_gral = 0.5, max_cb_N = 10000, min_cb_N = 50) :
    #, max_ci = 0.1, n_imp = 3, min_n_suc = 50, min_n_fail = 50, min_n_feat = 50
    # %% AVERAGE CONNECTIVITY :
    # LPFC to AVG BRAIN (FIG 2B)
    # AVG LPFC- AVG BRAIN  (not ploted)
    # Using square matrix of LPFC parcellation (X-X). All LPFC merged as a block toward all the rest of the ipsilateral brain (excluding LPFC).
    # Analysis to prove general efferent-afferent values of the ROI as one. AVG done in compute_atlas_brain_avg
    path_output_avg_all = os.path.join(output_folder, 'AVG_LPFC_AVG_all');os.makedirs(path_output_avg_all, exist_ok=True)

    # Left LPFC
    path_mat_avg_lpfc_avg_all_L_eff = os.path.join(base_path, 'lh_AVG_LPFC_AVG_all_eff', collapse_conf_name[0])
    path_mat_avg_lpfc_avg_all_L_aff = os.path.join(base_path, 'lh_AVG_LPFC_AVG_all_aff', collapse_conf_name[0])
    mat_all_avg_eff, N_all_avg_eff ,  I_all_avg_eff = matrices.atlas_mat(path_mat_avg_lpfc_avg_all_L_eff, [0,], [0,],path_output_avg_all, 'avg_lLPFC_brain_all_eff', flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    mat_all_avg_aff , N_all_avg_aff,  I_all_avg_aff  = matrices.atlas_mat(path_mat_avg_lpfc_avg_all_L_aff, [0,], [0,],path_output_avg_all, 'avg_lLPFC_brain_all_aff', flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    # Right LPFC
    path_mat_avg_lpfc_avg_all_R_eff = os.path.join(base_path, 'rh_AVG_LPFC_AVG_all_eff', collapse_conf_name[0])
    path_mat_avg_lpfc_avg_all_R_aff = os.path.join(base_path, 'rh_AVG_LPFC_AVG_all_aff', collapse_conf_name[0])
    mat_all_avg_eff, N_all_avg_eff ,  I_all_avg_eff = matrices.atlas_mat(path_mat_avg_lpfc_avg_all_R_eff, [0,], [0,],path_output_avg_all, 'avg_rLPFC_brain_all_eff', flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    mat_all_avg_aff , N_all_avg_aff,  I_all_avg_aff  = matrices.atlas_mat(path_mat_avg_lpfc_avg_all_R_aff, [0,], [0,],path_output_avg_all, 'avg_rLPFC_brain_all_aff', flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    # AVG LPFC parcels (Lausanne125) - BRAIN (FIG 2B)
    # ROI (LPFC) parcelled according to Lausanne2008-X towards the rest of the ipsilateral brain (excluding the ROI) merged as one.
    # Use mat xx, with x parcellatioin of roi (for the exclusion of the roi to be possible)
    path_output_avg = os.path.join(output_folder, 'AVG');os.makedirs(path_output_avg, exist_ok=True)

    # Left LPFC - all
    path_mat_avg_lpfc_avg_all_L_eff = os.path.join(base_path, 'lh_AVG_eff', collapse_conf_name[0])
    path_mat_avg_lpfc_avg_all_L_aff = os.path.join(base_path, 'lh_AVG_aff', collapse_conf_name[0])
    l_LPFC = [l for l in Definitions.roi_dlpfc_ifg_125 if l.startswith('lh')]

    mat_all_avg_eff, N_all_avg_eff ,  I_all_avg_eff = matrices.atlas_mat(path_mat_avg_lpfc_avg_all_L_eff, np.arange(0,len(l_LPFC)), [0,],path_output_avg, 'lLPFC_brain_all_eff', flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    mat_all_avg_aff , N_all_avg_aff,  I_all_avg_aff  = matrices.atlas_mat(path_mat_avg_lpfc_avg_all_L_aff, [0,], np.arange(0,len(l_LPFC)), path_output_avg, 'lLPFC_brain_all_aff', flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    # Right LPFC
    path_mat_avg_lpfc_avg_all_R_eff = os.path.join(base_path, 'rh_AVG_eff', collapse_conf_name[0])
    path_mat_avg_lpfc_avg_all_R_aff = os.path.join(base_path, 'rh_AVG_aff', collapse_conf_name[0])
    r_LPFC = [l for l in Definitions.roi_dlpfc_ifg_125 if l.startswith('rh')]
    mat_all_avg_eff, N_all_avg_eff ,  I_all_avg_eff = matrices.atlas_mat(path_mat_avg_lpfc_avg_all_R_eff, np.arange(0,len(r_LPFC)), [0,],path_output_avg, 'rLPFC_brain_all_eff', flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    mat_all_avg_aff , N_all_avg_aff,  I_all_avg_aff  = matrices.atlas_mat(path_mat_avg_lpfc_avg_all_R_aff, [0,],  np.arange(0,len(r_LPFC)), path_output_avg, 'rLPFC_brain_all_aff', flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    #%% PARCELLATION RESOLUTIONS (FIG1D-E-F)
    # Corrected p for square matrices of Lausanne parcellations 33, 125, 500.
    # Stimulation on selected DLPFC equivalent parcels
    path_resolutions = os.path.join(output_folder, 'Resolutions');os.makedirs(path_resolutions, exist_ok=True)

    #Resolution 1
    path_33_33 = os.path.join(base_path, res33 + '_' + res33)
    path_mats_33_33 = os.path.join(path_33_33, collapse_conf_name[0])
    label_stim =  "lh.rostralmiddlefrontal" #Lau33 stim parcel
    p_33, N_33, I_33 = matrices.atlas_mat(path_mats_33_33,  labels_rec_33.index(label_stim) , index_rec_33 , path_resolutions,f'stim33to33_{gral_tw[0]}_{gral_tw[1]}_ms', flag_delays= False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    p_stim_33, N_stim_33, I_stim_33= matrices.atlas_mat(path_mats_33_33,  labels_rec_33.index(label_stim) , labels_rec_33.index(label_stim) , path_resolutions,f'stim33tostim33_{gral_tw[0]}_{gral_tw[1]}_ms', flag_delays= False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    #Resolution 3
    path_125_125 = os.path.join(base_path, res125 + '_' + res125)
    path_mats_125_125 = os.path.join(path_125_125, collapse_conf_name[0])
    label_stim = "lh.rostralmiddlefrontal_1" #Lau125 stim parcel
    p_125, N_125, I_125= matrices.atlas_mat(path_mats_125_125,  labels_125.index(label_stim) , index_125 , path_resolutions,f'stim125to125_{gral_tw[0]}_{gral_tw[1]}_ms', flag_delays= False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    p_stim_125, N_stim_125, I_stim_125 = matrices.atlas_mat(path_mats_125_125,  labels_125.index(label_stim) , labels_125.index(label_stim) , path_resolutions,f'stim125tostim125_{gral_tw[0]}_{gral_tw[1]}_ms', flag_delays= False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    #Resolution 5
    path_500_500 = os.path.join(base_path, res500 + '_' + res500)
    path_mats_500_500 = os.path.join(path_500_500, collapse_conf_name[0])
    label_stim = "lh.rostralmiddlefrontal_22" #Lau500 stim parcel
    p_500, N_500, I_500= matrices.atlas_mat(path_mats_500_500,  labels_500.index(label_stim) , index_500 , path_resolutions,f'stim500to500_{gral_tw[0]}_{gral_tw[1]}_ms', flag_delays= False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    p_stim_500, N_stim_500, I_stim_500= matrices.atlas_mat(path_mats_500_500,  labels_500.index(label_stim) , labels_500.index(label_stim) , path_resolutions,f'stim500tostim500_{gral_tw[0]}_{gral_tw[1]}_ms', flag_delays= False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    #%% EFFERENT CONNECTIVITY (FIG 2E and FIG 2C - first row)  -  CONNECTIVITY FROM X to Y
        # Matrices of combined resolution xy. Fine resolution x over roi and bigger resolution y over all the brain.
        # From the raw original matrix we take the lines of roi stimulated parcels towards columns of interest (in our case all cortical parcels + Amygdala and Hyppocampus)
        # This is in the matrix of parcellations x to y
        # We use xx matrix for stimulations to the roi recorded by the roi.
    path_125_33 = os.path.join(base_path, res125 + '_' + res33)
    path_mats_125_33 = os.path.join(path_125_33, collapse_conf_name[0])
    path_125_125 = os.path.join(base_path, res125 + '_' + res125)
    path_mats_125_125 = os.path.join(path_125_125, collapse_conf_name[0])
    path_folder_output = os.path.join(output_folder, res125 + '_' + res33);os.makedirs(path_folder_output, exist_ok=True)

    #p_xy, N_xy, I_xy, D_med_xy, D_mean_xy, quant25_xy, quant75_xy = matrices.atlas_mat(path_mats_125_33,Definitions.index_dlpfc_ifg,index_rec_33, path_folder_output,'125to33_' + gral_tw,flag_delays=True, min_n_suc=N,min_n_fail=N, min_n_feat=N,min_n_impl=n_imp, max_ci=max_ci)
    #p_xy, N_xy, I_xy= matrices.atlas_mat(path_mats_125_33,Definitions.index_dlpfc_ifg,index_rec_33, path_folder_output,f'125to33_{gral_tw[0]}_{gral_tw[1]}_ms',flag_delays=False, N_thresh=N)
    p_xy, N_xy, I_xy, D_med_xy, D_mean_xy= matrices.atlas_mat(path_mats_125_33, Definitions.index_dlpfc_ifg, index_rec_33,path_folder_output, f'125to33_{gral_tw[0]}_{gral_tw[1]}_ms',flag_delays=True, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    not_nan_xy = np.sum(~np.isnan(p_xy));coverage_xy = not_nan_xy / p_xy.size
    print("Coverage eff xy: ", coverage_xy)
    coverage_file_path = os.path.join(path_folder_output, 'coverage_xy.txt')
    with open(coverage_file_path, 'w') as f: f.write(f"Coverage: {coverage_xy:.6f}\n")

    #p_xx, N_xx, I_xx, D_med_xx, D_mean_xx, quant25_xx, quant75_xx = matrices.atlas_mat(path_mats_125_125,Definitions.index_dlpfc_ifg,Definitions.index_dlpfc_ifg,path_folder_output,'125to125_' + gral_tw,flag_delays=True, min_n_suc=N,min_n_fail=N, min_n_feat=N, min_n_impl=n_imp, max_ci=max_ci)
    #p_xx, N_xx, I_xx = matrices.atlas_mat(path_mats_125_125,Definitions.index_dlpfc_ifg,Definitions.index_dlpfc_ifg,path_folder_output,f'125to125_{gral_tw[0]}_{gral_tw[1]}_ms',  flag_delays=False, N_thresh=N)
    p_xx, N_xx, I_xx,D_med_xx, D_mean_xx, = matrices.atlas_mat(path_mats_125_125, Definitions.index_dlpfc_ifg, Definitions.index_dlpfc_ifg,path_folder_output, f'125to125_{gral_tw[0]}_{gral_tw[1]}_ms',flag_delays=True, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    not_nan_xx = np.sum(~np.isnan(p_xx));coverage_xx = not_nan_xx / p_xx.size
    print("Coverage eff xx: ", coverage_xx)
    coverage_file_path = os.path.join(path_folder_output, 'coverage_xx.txt')
    with open(coverage_file_path, 'w') as f: f.write(f"Coverage: {coverage_xx:.6f}\n")

    # %% AFFERENT CONNECTIVITY (FIG 2E and FIG 2C - second row) -  CONNECTIVITY FROM Y to X
    # From the raw original matrix we take the lines of roi recorded parcels of stimulations done on columns of interest on the rest of the brain
    path_33_125 = os.path.join(base_path, res33 + '_' + res125)
    path_mats_33_125 = os.path.join(path_33_125, collapse_conf_name[0])
    # The square ones (xx) are repeted from eff, but this way its independent
    path_125_125 = os.path.join(base_path, res125 + '_' + res125)
    path_mats_125_125 = os.path.join(path_125_125, collapse_conf_name[0])
    path_folder_output_aff = os.path.join(output_folder, res33 + '_' + res125)
    os.makedirs(path_folder_output_aff, exist_ok=True)
    # probability
    #p_yx, N_yx, I_yx, D_med_yx, D_mean_yx, quant25_yx, quant75_yx = matrices.atlas_mat(path_mats_33_125, index_rec_33,Definitions.index_dlpfc_ifg,path_folder_output_aff,'33to125_' + gral_tw,flag_delays=True, min_n_suc=N, min_n_fail=N, min_n_feat=N,min_n_impl=n_imp, max_ci=max_ci)
    #p_yx, N_yx, I_yx= matrices.atlas_mat(path_mats_33_125, index_rec_33,Definitions.index_dlpfc_ifg,path_folder_output_aff,f'33to125_{gral_tw[0]}_{gral_tw[1]}_ms',flag_delays=False, N_thresh=N)
    p_yx, N_yx, I_yx, D_med_yx, D_mean_yx =  matrices.atlas_mat(path_mats_33_125, index_rec_33, Definitions.index_dlpfc_ifg,path_folder_output_aff, f'33to125_{gral_tw[0]}_{gral_tw[1]}_ms',flag_delays=True, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    not_nan_yx = np.sum(~np.isnan(p_yx));coverage_yx = not_nan_yx / p_yx.size
    print("Coverage aff xy: ", coverage_yx)
    coverage_file_path = os.path.join(path_folder_output_aff, 'coverage_yx.txt')
    with open(coverage_file_path, 'w') as f: f.write(f"Coverage: {coverage_yx:.6f}\n")
    #p_xx, N_xx, I_xx, D_med_xx, D_mean_xx, quant25_yx, quant75_xy = matrices.atlas_mat(path_mats_125_125,Definitions.index_dlpfc_ifg,Definitions.index_dlpfc_ifg,path_folder_output_aff,'125to125_' + gral_tw,flag_delays=True, min_n_suc=N,min_n_fail=N, min_n_feat=N,min_n_impl=n_imp, max_ci=max_ci)
    #p_xx, N_xx, I_xx= matrices.atlas_mat(path_mats_125_125,Definitions.index_dlpfc_ifg,Definitions.index_dlpfc_ifg,path_folder_output_aff,f'125to125_{gral_tw[0]}_{gral_tw[1]}_ms',flag_delays=False, N_thresh=N)
    p_xx, N_xx, I_xx, D_med_xx, D_mean_xx= matrices.atlas_mat(path_mats_125_125, Definitions.index_dlpfc_ifg, Definitions.index_dlpfc_ifg,path_folder_output_aff, f'125to125_{gral_tw[0]}_{gral_tw[1]}_ms',flag_delays=True, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    not_nan_xx = np.sum(~np.isnan(p_xx));coverage_xx = not_nan_xx / p_xx.size
    coverage_file_path = os.path.join(path_folder_output_aff, 'coverage_xx.txt')
    with open(coverage_file_path, 'w') as f: f.write(f"Coverage: {coverage_xx:.6f}\n")
    print("Coverage aff xx: ", coverage_xx)

    # %% MERGED ROIS: DLPFC & IFG. NUMBER of RECORDINGS (FIG 3A) DIRECT (FIG 3B) & INDIRECT (FIG 3C) CONNECTIVITY (probability and delays)
    # DLPFC merged (using Lausanne x) towards the rest of the brain in Lausanne y
    # DIRECT
    path_output_dlpfc = os.path.join(output_folder, 'Mean_ROI_Eff')
    os.makedirs(path_output_dlpfc, exist_ok=True)
    output_paths_dlpfc = os.path.join(path_output_dlpfc, f'Direct_connectivity_{dir_tw[0]}_{dir_tw[1]}_ms')
    os.makedirs(output_paths_dlpfc, exist_ok=True)

    # - DIRECT MATRICES
    #all
    path_dir_roi_33 = os.path.join(base_path, 'DLPFC_IFG_to_Lausanne2008-33')
    stim_parcels_merged = np.loadtxt(os.path.join(path_dir_roi_33, 'stim_parcels_merged.txt'), dtype=str) #stim parcels is roi merged (created when compute atlas : lh.dlpfc, rh.dlpfc, lh.ifg, rh.ifg
    stim_parcels_merged = stim_parcels_merged.tolist()
    index_stims = [stim_parcels_merged.index(i) for i in stim_parcels_merged]
    rec_parcels = np.loadtxt(os.path.join(path_dir_roi_33, 'rec_parcels_Lausanne2008-33.txt'), dtype=str)
    rec_parcels = rec_parcels.tolist()
    index_recs = [rec_parcels.index(i) for i in labels_rec_33] #only labels I need, index from real atlas list. #TODO: maybe compute directly only labels wanted
    # Mats
    path_dir_roi_33 = os.path.join(base_path, 'DLPFC_IFG_to_Lausanne2008-33', 'Direct', collapse_conf_name[1])
    # stim_parcels_merged_fig_names = [f'{l}_to_33' for l in stim_parcels_merged]
    p_xy_ldlpfc, N_xy_ldlpfc, I_xy_ldlpfc, D_med_xy_ldlpfc, D_mean_xy_ldlpfc = matrices.atlas_mat(path_dir_roi_33, index_stims,index_recs,output_paths_dlpfc,'roi_to_33', flag_delays=True, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    #local
    path_dir_roi_roi = os.path.join(base_path, 'DLPFC_IFG_to_DLPFC_IFG')
    stim_parcels_merged = np.loadtxt(os.path.join(path_dir_roi_roi,'stim_parcels_merged.txt'), dtype=str)  # stim parcels is roi merged (created when compute atlas : lh.dlpfc, rh.dlpfc, lh.ifg, rh.ifg
    stim_parcels_merged = stim_parcels_merged.tolist()
    index_stims = [stim_parcels_merged.index(i) for i in stim_parcels_merged]
    rec_parcels = np.loadtxt(os.path.join(path_dir_roi_roi, 'rec_parcels_merged.txt'), dtype=str)
    rec_parcels = rec_parcels.tolist()
    index_recs = [rec_parcels.index(i) for i in rec_parcels]  # just to avoid mistakes, it should be = to stim_parcels_merged
    path_dir_roi_roi = os.path.join(base_path, 'DLPFC_IFG_to_DLPFC_IFG', 'Direct', collapse_conf_name[1])
    # stim_parcels_merged_fig_names = [f'{l}_to_{l}' for l in stim_parcels_merged]
    p_xx_ldlpfc, N_xx_ldlpfc, I_xx_ldlpfc, D_med_xx_ldlpfc, D_mean_xx_ldlpfc = matrices.atlas_mat(path_dir_roi_roi,index_stims, index_stims, output_paths_dlpfc, 'roi_to_roi',flag_delays=True, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    # INDIRECT
    output_paths_dlpfc = os.path.join(path_output_dlpfc, f'Indirect_connectivity_{ind_tw[0]}_{ind_tw[1]}_ms')
    os.makedirs(output_paths_dlpfc, exist_ok=True)
    # all
    path_ind_roi_33 = os.path.join(base_path, 'DLPFC_IFG_to_Lausanne2008-33')
    stim_parcels_merged = np.loadtxt(os.path.join(path_ind_roi_33, 'stim_parcels_merged.txt'),dtype=str)  # stim parcels is roi merged (created when compute atlas : lh.dlpfc, rh.dlpfc, lh.ifg, rh.ifg
    stim_parcels_merged = stim_parcels_merged.tolist()
    index_stims = [stim_parcels_merged.index(i) for i in stim_parcels_merged]
    rec_parcels = np.loadtxt(os.path.join(path_ind_roi_33, 'rec_parcels_Lausanne2008-33.txt'), dtype=str)
    rec_parcels = rec_parcels.tolist()
    index_recs = [rec_parcels.index(i) for i in  labels_rec_33]  # only labels I need, index from real atlas list. #TODO: maybe compute directly only labels wanted
    # Mats
    path_ind_roi_33 = os.path.join(base_path, 'DLPFC_IFG_to_Lausanne2008-33', 'Indirect', collapse_conf_name[2])
    # stim_parcels_merged_fig_names = [f'{l}_to_33' for l in stim_parcels_merged]
    p_xy_ldlpfc, N_xy_ldlpfc, I_xy_ldlpfc, D_med_xy_ldlpfc, D_mean_xy_ldlpfc = matrices.atlas_mat(path_ind_roi_33,  index_stims,index_recs,output_paths_dlpfc, 'roi_to_33', flag_delays=True, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)
    # local
    path_ind_roi_roi = os.path.join(base_path, 'DLPFC_IFG_to_DLPFC_IFG')
    stim_parcels_merged = np.loadtxt(os.path.join(path_ind_roi_roi, 'stim_parcels_merged.txt'),dtype=str)  # stim parcels is roi merged (created when compute atlas : lh.dlpfc, rh.dlpfc, lh.ifg, rh.ifg
    stim_parcels_merged = stim_parcels_merged.tolist()
    index_stims = [stim_parcels_merged.index(i) for i in stim_parcels_merged]
    rec_parcels = np.loadtxt(os.path.join(path_ind_roi_roi, 'rec_parcels_merged.txt'), dtype=str)
    rec_parcels = rec_parcels.tolist()
    index_recs = [rec_parcels.index(i) for i in rec_parcels]  # just to avoid mistakes, it should be = to stim_parcels_merged
    path_ind_roi_roi = os.path.join(base_path, 'DLPFC_IFG_to_DLPFC_IFG', 'Indirect', collapse_conf_name[2])
    p_xx_ldlpfc, N_xx_ldlpfc, I_xx_ldlpfc, D_med_xx_ldlpfc, D_mean_xx_ldlpfc = matrices.atlas_mat(path_ind_roi_roi, index_stims, index_stims,output_paths_dlpfc,'roi_to_roi',flag_delays=True, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    # %% ROI CONNECTIVITY TO FUNCTIONAL NETWORKS (FIG 4C)
    #     Connectivity from roi (Lausanne x) to functional networks as defined by Yeo et al. We compute the functional networks as mergeds of Lausanne 250 parcels.
    path_folder_funct = os.path.join(output_folder, 'Functional_Networks')
    os.makedirs(path_folder_funct, exist_ok=True)

    path_fn_data =  os.path.join(base_path, 'Functional_Networks')
    stim_parcels = np.loadtxt(os.path.join(path_fn_data, 'stim_parcels.txt'),dtype=str)
    stim_parcels = stim_parcels.tolist()

    rec_fn = np.loadtxt(os.path.join(path_fn_data, 'rec_parcels_merged.txt'), dtype=str)
    rec_fn = rec_fn.tolist()

    path_fn_data = os.path.join(path_fn_data, collapse_conf_name[0])
    #Left fn
    index_stims_L = [stim_parcels.index(i) for i in stim_parcels if i.startswith('lh')]
    index_rec_fn_L = [i for i, name in enumerate(rec_fn) if name.startswith('lh.')]
    #_{gral_tw[0]}_{gral_tw[1]}ms - I used to save the tw it s computed with. But always general for fn
    p_net_eff, N_net_eff, I_net_eff = matrices.atlas_mat(path_fn_data, index_stims_L,index_rec_fn_L,path_folder_funct, 'Eff_roi_L_nets' , flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    #Right fn
    index_stims_R = [stim_parcels.index(i) for i in stim_parcels if i.startswith('rh')]
    index_rec_fn_R = [i for i, name in enumerate(rec_fn) if name.startswith('rh.')]
    p_net_eff, N_net_eff, I_net_eff = matrices.atlas_mat(path_fn_data, index_stims_R, index_rec_fn_R, path_folder_funct,'Eff_roi_R_nets', flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

    #Segments
    path_folder_funct_segments = os.path.join(output_folder, 'Functional_Networks', 'Segmented_roi')
    os.makedirs(path_folder_funct_segments, exist_ok=True)

    path_fn_data_segments = os.path.join(base_path, 'Functional_Networks', 'Segmented_roi')
    stim_parcels_merged = np.loadtxt(os.path.join(path_fn_data_segments, 'stim_parcels_merged.txt'),dtype=str)
    stim_parcels = stim_parcels_merged.tolist()

    rec_fn = np.loadtxt(os.path.join(path_fn_data_segments, 'rec_parcels_merged.txt'), dtype=str)
    rec_fn = rec_fn.tolist()
    path_fn_data_segments = os.path.join(path_fn_data_segments, collapse_conf_name[0])
    #left and right, save on mat per segment (adapt to plot bars)
    for s_id in range(0, len(stim_parcels)) :
        stim_segment_tmp = stim_parcels[s_id]
        suffix = stim_segment_tmp.split('_')[-1] if '_' in stim_segment_tmp else stim_segment_tmp.split('.')[-1]
        #intra hemi p
        if stim_segment_tmp.startswith('lh.') :
         index_rec_fn = [i for i, name in enumerate(rec_fn) if name.startswith('lh.')]
         mat_name = 'l_'
        else :
            index_rec_fn = [i for i, name in enumerate(rec_fn) if name.startswith('rh.')]
            mat_name = 'r_'
        #mat_name
        if suffix == 'ifg' : mat_name = mat_name + 'IFG_nets'
        else: mat_name = mat_name + suffix + '_DLPFC_nets'
        p_segment_net, N_segment_net, I_segment_net = matrices.atlas_mat(path_fn_data_segments, [s_id,],index_rec_fn,path_folder_funct_segments, mat_name, flag_delays=False, N_thresh=N_thresh, N_imp_thresh=N_imp_thresh)

def save_params_txt(output_folder,gral_tw,dir_tw,ind_tw,collapse_conf_name, N_thresh, N_imp_thresh, min_cb_dir, max_cb_dir,min_cb_ind,max_cb_ind,max_cb_gral,min_avg,max_avg,max_cb_p,max_p_gral,min_cb_N, max_cb_N):
    os.makedirs(output_folder, exist_ok=True)
    params_text = f"""
    Time Windows:
      General: {gral_tw}
      Direct: {dir_tw}
      Indirect: {ind_tw}
    Collapse Config Names:
      General: {collapse_conf_name[0]}
      Direct: {collapse_conf_name[1]}
      Indirect: {collapse_conf_name[2]}
    N recordings mask: {N_thresh}
    N implantations mask: {N_imp_thresh}
    Delays cb (ms):
      Direct min: {min_cb_dir}
      Direct max: {max_cb_dir}
      Indirect min: {min_cb_ind}
      Indirect max: {max_cb_ind}
      General max: {max_cb_gral}
    P avg cb:
      Min avg: {min_avg}
      Max avg: {max_avg}
    P thresholds:
      Max p for connections: {max_cb_p}
      Max p general: {max_p_gral}
    N cb:
      Min cb N: {min_cb_N}
      Max cb N: {max_cb_N}
    """
    with open(os.path.join(output_folder, 'params.txt'), 'w', encoding='utf-8') as f:  f.write(params_text)

def main():
    Zth = str(5)
    #Data : atlas generation output folders (containing matrices)
    base_path = r'F:\FTRACT\Data_LPFC_FTRACT'
    output_folder = os.path.join(base_path, f'Results_Zth{Zth}')

    gral_tw = [0,100] #'0_100ms'
    dir_tw = [0,50] #'0_50ms'  # time window of 'direct connections'
    ind_tw = [100, 400] #'100_400ms'  # time window of 'indirect connections'
    #collapse_conf_name general, direct, indirect
    collapse_conf_name = [f'max_peak_delay_{gral_tw[1]}__zth{Zth}___min_peak_delay_{gral_tw[0]}__zth{Zth}',f'max_peak_delay_{dir_tw[1]}__zth{Zth}___min_peak_delay_{dir_tw[0]}__zth{Zth}',f'max_peak_delay_{ind_tw[1]}__zth{Zth}___min_peak_delay_{ind_tw[0]}__zth{Zth}'] # Name of folder that saved the collapse, features inside.
    N_thresh=  50
    N_imp_thresh= 3
    # max_ci = 0.1 #not in use

    # delays
    # thresh_d = 20
    min_cb_dir = 20
    max_cb_dir = 35
    min_cb_ind = 150
    max_cb_ind = 200
    max_cb_gral = 80  # min 10
    # p avg
    max_avg = 0.175
    min_avg = 0.1
    # p
    max_cb_p = 0.25
    max_p_gral = 0.5
    # N
    max_cb_N = 10000
    min_cb_N = 50

    compute_mats(base_path, output_folder, collapse_conf_name, N_thresh, N_imp_thresh, gral_tw, dir_tw, ind_tw, min_cb_dir, max_cb_dir, min_cb_ind, max_cb_ind, max_cb_gral, max_avg, min_avg, max_cb_p, max_p_gral, max_cb_N, min_cb_N)

    save_params_txt(output_folder, gral_tw, dir_tw, ind_tw, collapse_conf_name, N_thresh, N_imp_thresh, min_cb_dir, max_cb_dir, min_cb_ind,max_cb_ind, max_cb_gral, min_avg, max_avg, max_cb_p, max_p_gral, min_cb_N, max_cb_N)


if __name__ == "__main__":
        main()
