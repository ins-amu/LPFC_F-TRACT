# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 16:54:03 2024

@author: avalos-alais.s
"""

#Ftract - DLPFC&IFG 
# Function integrating processes to compute the matrices of connctivity. 
# Takes original matrices (stimulated x recorded), do necessary selection of rows/colomns, and merges for each specific analysis. 
# Filters the resulting matrices accordign to statistical significancy. Plot results. 

import os
import matrices
from plotting import plotting
from plotting import plot_bars
from tools.combine_labels import combine_labels
from tools.def_funtional_250 import functional_networks_250
from tools import index_tool as idx
import numpy as np

def main():
    #%% DEFINITIONS -  X rows ; Y columns of the matrix
    L_s =  'Lausanne2008'
    res125 = 'Lausanne2008-125' 
    res33 = 'Lausanne2008-33'
    output_folder = r'C:\Users\avalos-alais.s\FTRACT_LPFC\Results'
    #Data : atlas generation output folders (containing matrices)
    path_general_folder = r'F:\Repo\Data\2024-01-10-10-52_filter_index_0__conf_0_100ms_2mm'
    
    path_gral_matrix_xx = os.path.join(path_general_folder, 'f-tract_atlas__' + res125 + '.npz')
    path_gral_matrix_yy = os.path.join(path_general_folder, 'f-tract_atlas__' + res33 + '.npz')
    path_gral_matrix_xy = os.path.join(path_general_folder, 'f-tract_atlas__' + res125 + '__' + res33 + '.npz')
    path_gral_matrix_yx = os.path.join(path_general_folder, 'f-tract_atlas__' + res33 + '__' + res125 + '.npz')
    
    path_N_contacts = os.path.join(output_folder, 'N_implanted_contacts')#pre-save matrices in this folder
    
    path_direct_folder =  r'F:\Repo\Data\2024-01-10-10-52_filter_index_0__conf_0_50ms_2mm'
   
    
    path_direct_matrix_xx = os.path.join(path_direct_folder, 'f-tract_atlas__' + res125 + '.npz')
    # path_direct_matrix_yy = os.path.join(path_direct_folder, 'f-tract_atlas__' + res33 + '.npz')
    path_direct_matrix_xy = os.path.join(path_direct_folder, 'f-tract_atlas__' + res125 + '__' + res33 + '.npz')
    # path_direct_matrix_yx = os.path.join(path_direct_folder, 'f-tract_atlas__' + res33 + '__' + res125 + '.npz')

    path_indirect_folder = r'F:\Repo\Data\2024-01-10-10-52_filter_index_0__conf_100_400ms_2mm'
    
    path_indirect_matrix_xx = os.path.join(path_indirect_folder,'f-tract_atlas__' +  res125 + '.npz')
    # path_indirect_matrix_yy = os.path.join(path_indirect_folder,'f-tract_atlas__' +  res33 + '.npz')
    path_indirect_matrix_xy = os.path.join(path_indirect_folder,'f-tract_atlas__' +  res125 + '__' + res33 + '.npz')
    # path_indirect_matrix_yx = os.path.join(path_indirect_folder,'f-tract_atlas__' +  res33 + '__' + res125 + '.npz')
    
    gral_tw = '0_100ms'     # time window for general analysis
    dir_tw  = '0_50ms'      # time window of 'direct connections' 
    ind_tw  = '100_400ms'   # time window of 'indirect connections'
    
    #path to object labels for plotting 
    labels_path = r'C:\Users\avalos-alais.s\FTRACT_LPFC\MNE-data\subjects\cvs_avg35_inMNI152' #r'C:\Users\avalos-alais.s\mne_data\MNE-sample-data\subjects\cvs_avg35_inMNI152 '
    labels_pathX =  os.path.join(labels_path + '-Lausanne125' , 'label')  #for plots
    labels_pathY =  os.path.join(labels_path + '-Lausanne33', 'label'   )  #for plots
    labels_path_500 =  os.path.join(labels_path + '-Lausanne500','label') #for resolutions
    labels_path250 =  os.path.join(labels_path + '-Lausanne250','label')  #for networks
    
    meshdirname = labels_path #Default in plotting functions, change parameter direction for implementation 
    
    #For functional network analysis
    res250 = 'Lausanne2008-250'
    path_general_matrix_125_250 = os.path.join(path_general_folder, 'f-tract_atlas__' + res125 + '__' + res250 + '.npz' )
    os.chdir(path_general_folder)
    labels_250 = [line.strip() for line in open(res250 + '.txt')]
    
    #Definitions of parcels of interest (Lausanne 125)
    roi_dlpfc_ifg = [ 
                    "lh.superiorfrontal_6",
                    "lh.rostralmiddlefrontal_1",
                    "lh.rostralmiddlefrontal_2",
                    "lh.rostralmiddlefrontal_3",
                    "lh.caudalmiddlefrontal_1",
                    "lh.caudalmiddlefrontal_2",
                    "lh.caudalmiddlefrontal_3",
                    "lh.parstriangularis_1",
                    "lh.parsopercularis_1",
                    "lh.parsopercularis_2",
                    
                    "rh.superiorfrontal_4",
                    "rh.superiorfrontal_8",
                    "rh.rostralmiddlefrontal_1",
                    "rh.rostralmiddlefrontal_2",
                    "rh.rostralmiddlefrontal_3",
                    "rh.caudalmiddlefrontal_1",
                    "rh.caudalmiddlefrontal_2",
                    "rh.caudalmiddlefrontal_3",
                    "rh.parstriangularis_1",
                    "rh.parstriangularis_2",
                    "rh.parsopercularis_1",
                    "rh.parsopercularis_2",
            ]
    roi_dlpfc =  [  "lh.superiorfrontal_6",
                    "lh.rostralmiddlefrontal_1",
                    "lh.rostralmiddlefrontal_2",
                    "lh.rostralmiddlefrontal_3",
                    "lh.caudalmiddlefrontal_1",
                    "lh.caudalmiddlefrontal_2",
                    "lh.caudalmiddlefrontal_3",
                
                    
                    "rh.superiorfrontal_4",
                    "rh.superiorfrontal_8",
                    "rh.rostralmiddlefrontal_1",
                    "rh.rostralmiddlefrontal_2",
                    "rh.rostralmiddlefrontal_3",
                    "rh.caudalmiddlefrontal_1",
                    "rh.caudalmiddlefrontal_2",
                    "rh.caudalmiddlefrontal_3",
      
        ]
    roi_ifg = [     "lh.parstriangularis_1",
                    "lh.parsopercularis_1",
                    "lh.parsopercularis_2",
                    
                    "rh.parstriangularis_1",
                    "rh.parstriangularis_2",
                    "rh.parsopercularis_1",
                    "rh.parsopercularis_2",
        ]
    
    os.chdir(path_general_folder)
    #labels of Y (33) parcellation
    labels_all_y = [line.strip() for line in open(res33 + '.txt')]
    labels_all_y = [s.replace("ctx-lh-", "lh.").replace("ctx-rh-", "rh.") for s in labels_all_y] #special for Lau33
    with open(res33 + '.txt', 'w') as file:
        file.write('\n'.join(labels_all_y) + '\n')
    labels_all  = [name for name in labels_all_y  if name.startswith("lh.") or name.startswith("rh.") or name== 'Left-Hippocampus' or name=='Left-Amygdala'  or name == 'Right-Hippocampus' or name == 'Right-Amygdala']
    #labels of X (125)  parcellation
    labels_all_x = [line.strip() for line in open(res125 + '.txt')]
    
    index_dlpfc_ifg = idx.get_idx(path_direct_folder, res125 + '.txt', roi_dlpfc_ifg)
    # index_dlpfc = idx.get_idx(path_direct_folder, res125 + '.txt', roi_dlpfc)
    # index_ifg = idx.get_idx(path_direct_folder, res125 + '.txt', roi_ifg)
    index_y = idx.get_idx(path_direct_folder, res33 + '.txt', labels_all)
    index_x = idx.get_idx(path_direct_folder, res125 + '.txt', labels_all_x)
   
    #delays
    # thresh_d = 20
    min_cb_dir = 20
    max_cb_dir = 35
    min_cb_ind = 150
    max_cb_ind = 200
    max_cb_gral = 80 #min 10 
    #p avg 
    max_avg  = 0.175
    min_avg  = 0.1
    #p
    max_cb_p = 0.25
    max_p_gral = 0.5
    #N 
    max_cb_N = 10000
    min_cb_N = 50
    
    #%% DEFINITION OF LPFC (FIG1G and FIG4A) + supplementary for right
    path_lpfc_def = os.path.join(output_folder, 'LPFC_Definition')
    os.makedirs(path_lpfc_def, exist_ok=True)
    plotting.plot_lpfc_definition(path_lpfc_def, labels_path, labels_all )
    #%% PARCELLATION RESOLUTIONS (FIG1D-E-F)
    # Corrected p for square matrices of Lausanne parcellations 33, 125, 500. 
    # Stimulation on selected DLPFC equivalent parcels
    path_resolutions = os.path.join(output_folder, 'Resolutions')
    os.makedirs(path_resolutions, exist_ok=True)
    path_resolutions_figs = os.path.join(path_resolutions, 'Figures')
    os.makedirs(path_resolutions_figs, exist_ok=True)
    
    label_stim =  "lh.rostralmiddlefrontal" #Lau33 stim parcel 
    index_stim_dlpfc = [(labels_all_y.index(label_stim),),]
    index_all_no_merge  = [(i,) for i in index_y]
    p_33, N_33, I_33= matrices.atlas_mat(path_gral_matrix_yy,  index_stim_dlpfc , index_all_no_merge , path_resolutions,'stim33to33_' + gral_tw, flag_delays= False)
    p_stim_33, N_stim_33, I_stim_33= matrices.atlas_mat(path_gral_matrix_yy,  index_stim_dlpfc , index_stim_dlpfc , path_resolutions,'stim33tostim33_' + gral_tw, flag_delays= False)
    plotting.plot_efferent (p_33, p_stim_33, [label_stim,], labels_all, labels_pathY, labels_pathY, path_resolutions_figs, max_p_gral, 0, 'plasma', border_print = True)
    
    label_stim = "lh.rostralmiddlefrontal_1" #Lau125 stim parcel
    index_stim_dlpfc = [(labels_all_x.index(label_stim),),]
    index_all_no_merge  = [(i,) for i in index_x]
    p_125, N_125, I_125= matrices.atlas_mat(path_gral_matrix_xx,  index_stim_dlpfc , index_all_no_merge , path_resolutions,'stim125to125_' + gral_tw, flag_delays= False)
    p_stim_125, N_stim_125, I_stim_125 = matrices.atlas_mat(path_gral_matrix_xx,  index_stim_dlpfc , index_stim_dlpfc , path_resolutions,'stim125tostim125_' + gral_tw, flag_delays= False)
    plotting.plot_efferent (p_125, p_stim_125, [label_stim,], labels_all_x, labels_pathX, labels_pathX, path_resolutions_figs, max_p_gral, 0, 'plasma', border_print = True)
    
    #hard coded 500 because it is the only time we use it 
    path_gral_matrix_500_500 = os.path.join(path_general_folder, 'f-tract_atlas__' + L_s + '-500' + '.npz')
    os.chdir(path_general_folder)
    labels_all_500 = [line.strip() for line in open('Lausanne2008-500.txt')]
    index_500 = idx.get_idx(path_general_folder, 'Lausanne2008-500.txt', labels_all_500)
    label_stim = "lh.rostralmiddlefrontal_22" #Lau500 stim parcel 
    index_stim_dlpfc = [(labels_all_500.index(label_stim),),]
    index_all_no_merge  = [(i,) for i in index_500]
    p_500, N_500, I_500= matrices.atlas_mat(path_gral_matrix_500_500,  index_stim_dlpfc , index_all_no_merge , path_resolutions,'stim500to500_' + gral_tw, flag_delays= False)
    p_stim_500, N_stim_500, I_stim_500= matrices.atlas_mat(path_gral_matrix_500_500,  index_stim_dlpfc , index_stim_dlpfc , path_resolutions,'stim500tostim500_' + gral_tw, flag_delays= False)
    plotting.plot_efferent (p_500, p_stim_500, [label_stim,], labels_all_500, labels_path_500, labels_path_500, path_resolutions_figs, max_p_gral, 0, 'plasma', border_print = True)
    
    #%% AVERAGE CONNECTIVITY : 
    # LPFC to AVG BRAIN (FIG 2B)
    # AVG LPFC- AVG BRAIN  (not ploted) 
    # Using square matrix of LPFC parcellation (X-X). All LPFC merged as a block toward all the rest of the ipsilateral brain (excluding LPFC). 
    # Analysis to prove general efferent-afferent values of the ROI as one. 
    path_output_avg_all = os.path.join(output_folder, 'AVG_all_LPFC')
    os.makedirs(path_output_avg_all, exist_ok=True)
    # Left LPFC
    lh_LPFC =  [l for l in roi_dlpfc_ifg if l.startswith('lh.') ]
    index_lLPFC = idx.get_idx(path_direct_folder, res125 + '.txt', lh_LPFC)
    idx_lLPFC_merge = [tuple(index_lLPFC)]
    #All LEFT brain without LPFC
    labels_all_left   = [l for l in labels_all_x if (( l.startswith('lh.') or  l.startswith ('Left') ) and l not in roi_dlpfc_ifg)]
    index_all_left_merge  = [tuple(idx.get_idx(path_general_folder, res125 + '.txt', labels_all_left))] 
    
    mat_all_avg_eff = matrices.atlas_mat(path_gral_matrix_xx, idx_lLPFC_merge , index_all_left_merge , path_output_avg_all, 'avg_lLPFC_eff', flag_delays=False)
    mat_all_avg_aff = matrices.atlas_mat(path_gral_matrix_xx, idx_lLPFC_merge , index_all_left_merge , path_output_avg_all, 'avg_lLPFC_aff', flag_delays=False, Eff_flag=False)
    
    # Right LPFC
    rh_LPFC =  [l for l in roi_dlpfc_ifg if l.startswith('rh.') ]
    index_rLPFC = idx.get_idx(path_direct_folder, res125 + '.txt', rh_LPFC)
    idx_rLPFC_merge = [tuple(index_rLPFC)]
    #All RIGHT brain without LPFC
    labels_all_right   = [l for l in labels_all_x if (( l.startswith('rh.') or  l.startswith ('Right') ) and l not in roi_dlpfc_ifg)]
    index_all_right_merge  = [tuple(idx.get_idx(path_general_folder, res125 + '.txt', labels_all_right))] 
    
    mat_all_avg_eff = matrices.atlas_mat(path_gral_matrix_xx, idx_rLPFC_merge , index_all_right_merge , path_output_avg_all, 'avg_rLPFC_eff', flag_delays=False)
    mat_all_avg_aff = matrices.atlas_mat(path_gral_matrix_xx, idx_rLPFC_merge , index_all_right_merge , path_output_avg_all, 'avg_rLPFC_aff', flag_delays=False, Eff_flag=False)
    
    # AVG LPFC parcels (Lausanne125) - BRAIN (FIG 2B)
    # ROI (LPFC) parcelled according to Lausanne2008-X towards the rest of the ipsilateral brain (excluding the ROI) merged as one. 
    # Use mat xx, with x parcellatioin of roi (for the exclusion of the roi to be possible)
    path_output_avg = os.path.join(output_folder, 'AVG')
    os.makedirs(path_output_avg, exist_ok=True)
    
    #LEFT
    roi_dlpfc_ifg_l = [l for l in roi_dlpfc_ifg if l.startswith('lh.')]
    index_dlpfc_ifg_l = idx.get_idx(path_general_folder, res125 + '.txt', roi_dlpfc_ifg_l)
    dlpfc_ifg_L = 'lh.DLPFC_IFG' #creo que esto no lo use?? 
    label_dlpfc_ifg_L  = combine_labels(roi_dlpfc_ifg_l, labels_pathX ,dlpfc_ifg_L, labels_pathX)
    index_dlpfc_ifg_l_no_merge = [(i,) for i in index_dlpfc_ifg_l]
    labels_all_left   = [l for l in labels_all_x if (( l.startswith('lh.') or  l.startswith ('Left') ) and l not in roi_dlpfc_ifg)]
    index_all_left_merge  = [tuple(idx.get_idx(path_general_folder, res125 + '.txt', labels_all_left))] 
    p_l_roi_all_eff, N_l_roi_all_eff, I_l_roi_all_eff = matrices.atlas_mat(path_gral_matrix_xx, index_dlpfc_ifg_l_no_merge , index_all_left_merge , path_output_avg, 'avg_L_roi_eff', flag_delays=False)
    p_l_roi_all_aff, N_l_roi_all_aff, I_l_roi_all_aff = matrices.atlas_mat(path_gral_matrix_xx, index_dlpfc_ifg_l_no_merge , index_all_left_merge , path_output_avg, 'avg_L_roi_aff', flag_delays=False, Eff_flag=False)
      
    #RIGHT
    roi_dlpfc_ifg_r = [l for l in roi_dlpfc_ifg if l.startswith('rh.')]
    index_dlpfc_ifg_r = idx.get_idx(path_general_folder, res125 + '.txt', roi_dlpfc_ifg_r)
    dlpfc_ifg_R = 'rh.DLPFC_IFG'
    label_dlpfc_ifg_R  = combine_labels(roi_dlpfc_ifg_r, labels_pathX ,dlpfc_ifg_R, labels_pathX)
    index_dlpfc_ifg_r_no_merge = [ (i,) for i in index_dlpfc_ifg_r]
    labels_all_right  = [l for l in labels_all_x if ((l.startswith('rh.') or  l.startswith ('Right')) and l not in roi_dlpfc_ifg)]
    index_all_right_merge  = [tuple(idx.get_idx(path_general_folder, res125 + '.txt', labels_all_right))]
    p_r_roi_all_eff, N_r_roi_all_eff, I_r_roi_all_eff = matrices.atlas_mat(path_gral_matrix_xx, index_dlpfc_ifg_r_no_merge , index_all_right_merge , path_output_avg, 'avg_R_roi_eff', flag_delays=False)
    p_r_roi_all_aff, N_r_roi_all_aff, I_r_roi_all_aff = matrices.atlas_mat(path_gral_matrix_xx, index_dlpfc_ifg_r_no_merge , index_all_right_merge , path_output_avg, 'avg_R_roi_aff', flag_delays=False, Eff_flag=False)
    
    #AVG PLOT [LEFT, RIGHT]
    avg_eff = np.vstack((p_l_roi_all_eff, p_r_roi_all_eff)) 
    avg_aff = np.hstack((p_l_roi_all_aff, p_r_roi_all_aff)) 
    avg_aff = avg_aff.transpose() 
    
    path_output_avg_figs = os.path.join(path_output_avg, 'Figures')
    os.makedirs(path_output_avg_figs, exist_ok=True)
    #Plot Efferent/Afferent average connectivity of roi in each hemisphere
    plotting.plot_mean_p(avg_eff, roi_dlpfc_ifg, labels_pathX, path_output_avg_figs, 'AVG_EFF.svg', 'Average Efferent Connectivity',  vmax = max_avg, vmin = min_avg , cbar = 'plasma', nan_color = '#c8c8c8')#light gray 
    plotting.plot_mean_p(avg_aff, roi_dlpfc_ifg, labels_pathX, path_output_avg_figs, 'AVG_AFF.svg', 'Average Afferent Connectivity',  vmax = max_avg, vmin = min_avg , cbar = 'plasma' , nan_color = '#c8c8c8')
    

 
    #%% NUMBER OF IMPLANTED CONTACTS (FIG 2A)
    # Vectors extracted separatly from data. Brain plot of number of implanted contacts for Lausanne33 resolution, with Lausanne125 overlapped over the LPFC roi
    os.chdir(path_N_contacts)
    contacts_33 = np.loadtxt(res33 + '.txt') 
    c33= contacts_33.reshape(1,len(contacts_33))
    c33 = c33[:, index_y]
    np.savetxt('contacts_33.txt' , c33)
   
    contacts_125 = np.loadtxt(res125+ '.txt')
    c125= contacts_125.reshape(1,len(contacts_125))
    c125 = c125[:,index_dlpfc_ifg]
    np.savetxt('contacts_125.txt' , c125)
    path_folder_output_N_contacts = os.path.join(path_N_contacts, 'Figures')
    os.makedirs(path_folder_output_N_contacts, exist_ok=True)
    
    roi_dlpfc_ifg_l = [l for l in roi_dlpfc_ifg if l.startswith('lh.') ]
    roi_dlpfc_ifg_r = [l for l in roi_dlpfc_ifg if l.startswith('rh.') ]
    
    plotting.plot_efferent (c33, c125[0:len(roi_dlpfc_ifg_l)], roi_dlpfc_ifg_l, labels_all, labels_pathX, labels_pathY, path_folder_output_N_contacts, 3000, 0, 'coolwarm', border_print = True,  labels_fig_names = ['LH Number of Contacts',] )
    plotting.plot_efferent (c33, c125[:,len(roi_dlpfc_ifg_l)::], roi_dlpfc_ifg_r, labels_all, labels_pathX, labels_pathY, path_folder_output_N_contacts, 3000, 0, 'coolwarm', border_print = True,  labels_fig_names = ['RH Number of Contacts',] )
    
    #%% EFFERENT CONNECTIVITY (FIG 2E and FIG 2C - first row)  -  CONNECTIVITY FROM X to Y 
    # Matrices of combined resolution xy. Fine resolution x over roi and bigger resolution y over all the brain. 
    # From the raw original matrix we take the lines of roi stimulated parcels towards colmuns of interest (in our case all cortical parcels + Amygdala and Hyppocampus)
    # This is in the matrix of parcellations x to y
    # We use xx matrix for stimulations to the roi recorded by the roi.
    
    path_folder_output = os.path.join(output_folder, res125 + '__' + res33)
    os.makedirs(path_folder_output, exist_ok=True)
    index_dlpfc_ifg_no_merge = [(i,) for i in index_dlpfc_ifg]
    index_y_no_merge  = [(i,) for i in index_y]
    p_xy, N_xy, I_xy, D_med_xy, D_mean_xy= matrices.atlas_mat(path_gral_matrix_xy, index_dlpfc_ifg_no_merge , index_y_no_merge , path_folder_output, '125to33_'+ gral_tw, flag_delays= True)
    p_xx, N_xx, I_xx, D_med_xx, D_mean_xx= matrices.atlas_mat(path_gral_matrix_xx,  index_dlpfc_ifg_no_merge , index_dlpfc_ifg_no_merge , path_folder_output,'125to125_'+gral_tw, flag_delays= True)
    
    # EFFERENT PLOT 125-33
    path_folder_output_figs = os.path.join(path_folder_output, 'Figures')
    path_folder_output_p = os.path.join(path_folder_output_figs, 'Probability_' + gral_tw)
    os.makedirs(path_folder_output_p, exist_ok=True)
    path_folder_output_median_pd = os.path.join(path_folder_output_figs, 'Median_Peak_Delays_' + gral_tw)
    os.makedirs(path_folder_output_median_pd, exist_ok=True)
    path_folder_output_mean_pd = os.path.join(path_folder_output_figs, 'Mean_Peak_Delays_' + gral_tw)
    os.makedirs(path_folder_output_mean_pd, exist_ok=True)
    #plot xy matrix and the xx matrix on top of it : probability and delays. 
    plotting.plot_efferent (p_xy, p_xx, roi_dlpfc_ifg, labels_all, labels_pathX, labels_pathY, path_folder_output_p, max_p_gral, 0, 'plasma')
    # plotting.plot_efferent (D_med_xy, D_med_xx, roi_dlpfc_ifg, labels_all, labels_pathX, labels_pathY, path_folder_output_median_pd, max_cb_gral, 10, 'viridis')
    # plotting.plot_efferent (D_med_xy, D_mean_xx, roi_dlpfc_ifg, labels_all, labels_pathX, labels_pathY, path_folder_output_mean_pd, max_cb_gral, 10, 'viridis')
    
    # EFFERENT MATRIX PLOT 
    p_eff_mat = np.hstack((p_xx, p_xy) ) # all rois(125) in rows to : 125 roi columns + 33 all columns 
    plotting.matrix_plot(p_eff_mat, roi_dlpfc_ifg, roi_dlpfc_ifg + labels_all , 'Probability of Efferent Connectivity', 'p_eff_125_125_33.svg', path_folder_output_figs, colorm = 'plasma', vmax = max_p_gral)
    
    #%% AFFERENT CONNECTIVITY (FIG 2E and FIG 2C - second row) -  CONNECTIVITY FROM Y to X 
    # From the raw original matrix we take the lines of roi recorded parcels of stimulations done on colmuns of interest on the rest of the brain 
    path_folder_output_aff = os.path.join(output_folder, res33 + '__' + res125)
    os.makedirs(path_folder_output_aff, exist_ok=True)
    #probability
    index_dlpfc_ifg_no_merge = [(i,) for i in index_dlpfc_ifg]
    index_y_no_merge  = [(i,) for i in index_y]
    p_yx, N_yx, I_yx, D_med_yx, D_mean_yx= matrices.atlas_mat(path_gral_matrix_yx, index_dlpfc_ifg_no_merge , index_y_no_merge , path_folder_output_aff, '33to125_' + gral_tw, flag_delays= True, Eff_flag= False)
    p_xx, N_xx, I_xx, D_med_xx, D_mean_xx= matrices.atlas_mat(path_gral_matrix_xx,  index_dlpfc_ifg_no_merge , index_dlpfc_ifg_no_merge , path_folder_output_aff,'125to125_'+gral_tw, flag_delays= True) 
    #The square ones (xx) are repeted from eff, but this way its independent
    
    #AFFERENT PLOT
    #Afferent plot : plot_efferent with transpose mat_all
    #For efferent the parcel on the title is the stim one, and for aff is the rec one.
    path_folder_output_aff = os.path.join(path_folder_output_aff, 'Figures')
    path_folder_output_aff_p = os.path.join(path_folder_output_aff, 'Probability_' + gral_tw)
    os.makedirs(path_folder_output_aff_p, exist_ok=True)
    path_folder_output_aff_median_pd = os.path.join(path_folder_output_aff, 'Median_Peak_Delays_' + gral_tw)
    os.makedirs(path_folder_output_aff_median_pd, exist_ok=True)
    path_folder_output_aff_mean_pd = os.path.join(path_folder_output_aff, 'Mean_Peak_Delays_' + gral_tw)
    os.makedirs(path_folder_output_aff_mean_pd, exist_ok=True)
    
    plotting.plot_efferent (p_yx.transpose(), p_xx, roi_dlpfc_ifg, labels_all, labels_pathX, labels_pathY, path_folder_output_aff_p,  max_p_gral, 0, 'plasma' )
    # plotting.plot_efferent (D_med_yx.transpose(), D_med_xx, roi_dlpfc_ifg, labels_all, labels_pathX, labels_pathY, path_folder_output_aff_median_pd, max_cb_gral, 10, 'viridis' )
    # plotting.plot_efferent (D_med_yx.transpose(), D_mean_xx, roi_dlpfc_ifg, labels_all, labels_pathX, labels_pathY, path_folder_output_aff_mean_pd, max_cb_gral, 10, 'viridis' )
   
    #MATRIX PLOT 
    p_aff_mat = np.vstack((p_xx, p_yx) ) # all rois(125) in rows to : 125 roi columns + 33 all columns 
    plotting.matrix_plot(p_aff_mat.transpose(), roi_dlpfc_ifg, roi_dlpfc_ifg + labels_all , 'Probability of Afferent Connectivity', 'p_aff_125_125_33.svg', path_folder_output_aff, colorm = 'plasma', vmax = max_p_gral)
    
    
    #%% MERGED ROIS: DLPFC & IFG. NUMBER of RECORDINGS (FIG 3A) DIRECT (FIG 3B) & INDIRECT (FIG 3C) CONNECTIVITY (probability and delays)
    # DLPFC merged (using Lausanne x) towards the rest of the brain in Lausanne y
    # DIRECT 
    path_output_dlpfc = os.path.join(output_folder, 'Mean_DLPFC_Eff')
    os.makedirs(path_output_dlpfc,exist_ok=True)
    output_paths_dlpfc = os.path.join(path_output_dlpfc, 'Diect_connectivity_'+ dir_tw )
    os.makedirs(output_paths_dlpfc,exist_ok=True)
    output_paths_dlpfc_figs = os.path.join(output_paths_dlpfc, 'Figures')
    os.makedirs(output_paths_dlpfc_figs,exist_ok=True)
    output_fig_paths_dlpfc_p = os.path.join(output_paths_dlpfc_figs, 'Probability')
    os.makedirs(output_fig_paths_dlpfc_p,exist_ok=True)
    output_fig_paths_dlpfc_N = os.path.join(output_paths_dlpfc_figs, 'N_stims')
    os.makedirs(output_fig_paths_dlpfc_N,exist_ok=True)
    output_fig_paths_dlpfc_delays = os.path.join(output_paths_dlpfc_figs, 'Peak_delays')
    os.makedirs(output_fig_paths_dlpfc_delays,exist_ok=True)
    # output_fig_paths_dlpfc_median = os.path.join(output_fig_paths_dlpfc_delays, 'Median')
    # os.makedirs(output_fig_paths_dlpfc_median,exist_ok=True)
    output_fig_paths_dlpfc_mean = os.path.join(output_fig_paths_dlpfc_delays, 'Mean')
    os.makedirs(output_fig_paths_dlpfc_mean,exist_ok=True)
    
    #LEFT DLPFC 
    
    roi_dlpfc_l = [l for l in roi_dlpfc if l.startswith('lh.') ]
    index_dlpfc_l = idx.get_idx(path_direct_folder, res125 + '.txt', roi_dlpfc_l)
    Ldlpfc_merge = 'lh.DLPFC'
    label_dlpfc_merged  = combine_labels(roi_dlpfc_l, labels_pathX , Ldlpfc_merge, path_output_dlpfc)
    idx_dlpfc_l_merge = [tuple(index_dlpfc_l)]
    index_y_no_merge  = [(i,) for i in index_y]
    #- DIRECT MATRICES
    p_xy_ldlpfc, N_xy_ldlpfc, I_xy_ldlpfc, D_med_xy_ldlpfc, D_mean_xy_ldlpfc = matrices.atlas_mat(path_direct_matrix_xy, idx_dlpfc_l_merge , index_y_no_merge , output_paths_dlpfc, 'lDLPFCto33', flag_delays= True)
    p_xx_ldlpfc, N_xx_ldlpfc, I_xx_ldlpfc, D_med_xx_ldlpfc, D_mean_xx_ldlpfc = matrices.atlas_mat(path_direct_matrix_xx,  idx_dlpfc_l_merge , idx_dlpfc_l_merge , output_paths_dlpfc,'lDLPFCtolDLPFC', flag_delays= True)
    #- PLOT
    #   N recordings
    plotting.plot_efferent (N_xy_ldlpfc, N_xx_ldlpfc, [Ldlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_N, max_cb_N, min_cb_N, 'coolwarm', border_print = True, labels_fig_names = ['LH N Recordings',])
    #   p
    plotting.plot_efferent (p_xy_ldlpfc, p_xx_ldlpfc, [Ldlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_p, max_cb_p, 0, 'plasma', labels_fig_names = ['LH Probability of Direct Connectivity',])
    #   Delays
    # plotting.plot_efferent (D_med_xy_ldlpfc, D_med_xx_ldlpfc, [Ldlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_median, max_cb_dir,min_cb_dir, 'viridis',  labels_fig_names = ['LH Median Peak Delays',])
    plotting.plot_efferent (D_mean_xy_ldlpfc, D_mean_xx_ldlpfc, [Ldlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_mean, max_cb_dir, min_cb_dir, 'viridis', labels_fig_names = ['LH Mean Peak Delays',])
       
    #RIGHT DLPFC 
    
    roi_dlpfc_r = [l for l in roi_dlpfc if l.startswith('rh.') ]
    index_dlpfc_r = idx.get_idx(path_direct_folder, res125 + '.txt', roi_dlpfc_r)
    Rdlpfc_merge = 'rh.DLPFC'
    label_dlpfc_merged  = combine_labels(roi_dlpfc_r, labels_pathX , Rdlpfc_merge, path_output_dlpfc)
    idx_dlpfc_r_merge = [tuple(index_dlpfc_r)]
    index_y_no_merge  = [(i,) for i in index_y]
    #- DIRECT MATRICES
    p_xy_rdlpfc, N_xy_rdlpfc, I_xy_rdlpfc, D_med_xy_rdlpfc, D_mean_xy_rdlpfc = matrices.atlas_mat(path_direct_matrix_xy, idx_dlpfc_r_merge , index_y_no_merge , output_paths_dlpfc, 'rDLPFCto33', flag_delays= True)
    p_xx_rdlpfc, N_xx_rdlpfc, I_xx_rdlpfc, D_med_xx_rdlpfc, D_mean_xx_rdlpfc = matrices.atlas_mat(path_direct_matrix_xx,  idx_dlpfc_r_merge , idx_dlpfc_r_merge , output_paths_dlpfc,'rDLPFCtorDLPFC', flag_delays= True)
    #- PLOT
    #   N recordings
    plotting.plot_efferent (N_xy_rdlpfc, N_xx_rdlpfc, [Rdlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_N, max_cb_N, min_cb_N, 'coolwarm', border_print = True ,labels_fig_names = ['RH N Recordings',])
    #   p 
    plotting.plot_efferent (p_xy_rdlpfc, p_xx_rdlpfc, [Rdlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_p, max_cb_p, 0, 'plasma', labels_fig_names = ['RH Probability of Direct Connectivity',])
    #   Delays
    # plotting.plot_efferent (D_med_xy_rdlpfc, D_med_xx_rdlpfc, [Rdlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_median, max_cb_dir,min_cb_dir, 'viridis',  labels_fig_names = ['RH Median Peak Delays',])
    plotting.plot_efferent (D_mean_xy_rdlpfc, D_mean_xx_rdlpfc, [Rdlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_mean, max_cb_dir, min_cb_dir, 'viridis', labels_fig_names = ['RH Mean Peak Delays',])
    
    #INDIRECT
    output_paths_dlpfc = os.path.join(path_output_dlpfc, 'Indirect_connectivity_' + ind_tw)
    os.makedirs(output_paths_dlpfc,exist_ok=True)
    output_paths_dlpfc_figs = os.path.join(output_paths_dlpfc, 'Figures')
    os.makedirs(output_paths_dlpfc_figs,exist_ok=True)
    output_fig_paths_dlpfc_p = os.path.join(output_paths_dlpfc_figs, 'Probability')
    os.makedirs(output_fig_paths_dlpfc_p,exist_ok=True)
    output_fig_paths_dlpfc_N = os.path.join(output_paths_dlpfc_figs, 'N_stims')
    os.makedirs(output_fig_paths_dlpfc_N,exist_ok=True)
    output_fig_paths_dlpfc_delays = os.path.join(output_paths_dlpfc_figs, 'Peak_delays')
    os.makedirs(output_fig_paths_dlpfc_delays,exist_ok=True)
    # output_fig_paths_dlpfc_median = os.path.join(output_fig_paths_dlpfc_delays, 'Median')
    # os.makedirs(output_fig_paths_dlpfc_median,exist_ok=True)
    output_fig_paths_dlpfc_mean = os.path.join(output_fig_paths_dlpfc_delays, 'Mean')
    os.makedirs(output_fig_paths_dlpfc_mean,exist_ok=True)
    
    #LEFT DLPFC - INDIRECT MATRICES 
    p_xy_ldlpfc, N_xy_ldlpfc, I_xy_ldlpfc, D_med_xy_ldlpfc, D_mean_xy_ldlpfc = matrices.atlas_mat(path_indirect_matrix_xy, idx_dlpfc_l_merge , index_y_no_merge , output_paths_dlpfc, 'lDLPFCto33', flag_delays= True)
    p_xx_ldlpfc, N_xx_ldlpfc, I_xx_ldlpfc, D_med_xx_ldlpfc, D_mean_xx_ldlpfc = matrices.atlas_mat(path_indirect_matrix_xx,  idx_dlpfc_l_merge , idx_dlpfc_l_merge , output_paths_dlpfc, 'lDLPFCtolDLPFC', flag_delays= True)
    
    #-PLOT
    #   N recordings - indirect N recs = direct N recs 
    # plotting.plot_efferent (N_xy_ldlpfc, N_xx_ldlpfc, [Ldlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_N, max_cb_N, 50, 'coolwarm',  border_print = True, labels_fig_names = ['LH N Recordings',])
    #   p
    plotting.plot_efferent (p_xy_ldlpfc, p_xx_ldlpfc, [Ldlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_p, max_cb_p , 0, 'plasma',  labels_fig_names = ['LH Probability of Indirect Connectivity',] )
    #   Delays
    # plotting.plot_efferent (D_med_xy_ldlpfc, D_med_xx_ldlpfc, [Ldlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_median, max_cb_ind, min_cb_ind, 'viridis',  labels_fig_names = ['LH Median Peak Delays',])
    plotting.plot_efferent (D_mean_xy_ldlpfc, D_mean_xx_ldlpfc, [Ldlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_mean, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names = ['LH Mean Peak Delays',])
    
    #%RIGHT DLPFC - INDIRECT MATRICES 
    p_xy_rdlpfc, N_xy_rdlpfc, I_xy_rdlpfc, D_med_xy_rdlpfc, D_mean_xy_rdlpfc = matrices.atlas_mat(path_indirect_matrix_xy, idx_dlpfc_r_merge , index_y_no_merge , output_paths_dlpfc, 'rDLPFCto33', flag_delays= True)
    p_xx_rdlpfc, N_xx_rdlpfc, I_xx_rdlpfc, D_med_xx_rdlpfc, D_mean_xx_rdlpfc = matrices.atlas_mat(path_indirect_matrix_xx,  idx_dlpfc_r_merge , idx_dlpfc_r_merge , output_paths_dlpfc, 'rDLPFCtorDLPFC', flag_delays= True)
    #-PLOT
    #   N recordings
    # plotting.plot_efferent (N_xy_rdlpfc, N_xx_rdlpfc, [Rdlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_N, max_cb_N,min_cb_N, 'coolwarm', border_print = True, labels_fig_names = ['RH N Recordings',])
    #   p
    plotting.plot_efferent (p_xy_rdlpfc, p_xx_rdlpfc, [Rdlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_p, max_cb_p , 0, 'plasma',  labels_fig_names = ['RH Probability of Indirect Connectivity',] )
    #   Delays
    # plotting.plot_efferent (D_med_xy_rdlpfc, D_med_xx_rdlpfc, [Rdlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_median, max_cb_ind, min_cb_ind, 'viridis',  labels_fig_names = ['RH Median Peak Delays',])
    plotting.plot_efferent (D_mean_xy_rdlpfc, D_mean_xx_rdlpfc, [Rdlpfc_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_dlpfc_mean, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names = ['RH Mean Peak Delays',])
   
    #IFG merged
    #DIRECT
    path_output_ifg = os.path.join(output_folder, 'Mean_IFG_Eff')
    os.makedirs(path_output_ifg,exist_ok=True)
    output_paths_ifg= os.path.join(path_output_ifg, 'Diect_connectivity_'+ dir_tw)
    os.makedirs(output_paths_ifg,exist_ok=True)
    output_paths_ifg_figs = os.path.join(output_paths_ifg, 'Figures')
    os.makedirs(output_paths_ifg_figs,exist_ok=True)
    output_paths_ifg_figs_p = os.path.join(output_paths_ifg_figs, 'Probability')
    os.makedirs(output_paths_ifg_figs_p,exist_ok=True)
    output_fig_paths_ifg_N = os.path.join(output_paths_ifg_figs, 'N_stims')
    os.makedirs(output_fig_paths_ifg_N,exist_ok=True)
    output_fig_paths_ifg_delays = os.path.join(output_paths_ifg_figs, 'Peak_delays')
    os.makedirs(output_fig_paths_ifg_delays,exist_ok=True)
    # output_fig_paths_ifg_median = os.path.join(output_fig_paths_ifg_delays, 'Median')
    # os.makedirs(output_fig_paths_ifg_median,exist_ok=True)
    output_fig_paths_ifg_mean = os.path.join(output_fig_paths_ifg_delays, 'Mean')
    os.makedirs(output_fig_paths_ifg_mean,exist_ok=True)
    
    #LEFT IFG
    roi_ifg_l = [l for l in roi_ifg if l.startswith('lh.') ]
    index_ifg_l = idx.get_idx(path_direct_folder, res125 + '.txt', roi_ifg_l)
    lifg_merge = 'lh.IFG'
    label_ifg_merged  = combine_labels(roi_ifg_l, labels_pathX, lifg_merge, path_output_ifg)
    idx_ifg_l_merge = [tuple(index_ifg_l)]
    index_y_no_merge  = [(i,) for i in index_y]
    #-DIRECT MATRICES 
    p_xy_lifg, N_xy_lifg, I_xy_lifg, D_med_xy_lifg, D_mean_xy_lifg = matrices.atlas_mat(path_direct_matrix_xy, idx_ifg_l_merge , index_y_no_merge , output_paths_ifg, 'lIFGto33', flag_delays=True)
    p_xx_lifg, N_xx_lifg, I_xx_lifg, D_med_xx_lifg, D_mean_xx_lifg = matrices.atlas_mat(path_direct_matrix_xx,  idx_ifg_l_merge , idx_ifg_l_merge , output_paths_ifg,'lIFGtolIFG', flag_delays=True)
    #-PLOT
    plotting.plot_efferent (N_xy_lifg, N_xx_lifg, [lifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_N, max_cb_N, min_cb_N, 'coolwarm', labels_fig_names = ['LH N Recordings',], border_print= True)
    plotting.plot_efferent(p_xy_lifg, p_xx_lifg, [lifg_merge,], labels_all, labels_pathX, labels_pathY, output_paths_ifg_figs_p, max_cb_p, 0, 'plasma', labels_fig_names=['LH Probability of Direct Connectivity',])
    # plotting.plot_efferent (D_med_xy_lifg, D_med_xx_lifg, [lifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_median, max_cb_dir, min_cb_dir, 'viridis',  labels_fig_names = ['LH Median Peak Delays',])
    plotting.plot_efferent (D_mean_xy_lifg, D_mean_xx_lifg, [lifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_mean, max_cb_dir, min_cb_dir, 'viridis', labels_fig_names = ['LH Mean Peak Delays',])
    
    #RIGHT IFG
    roi_ifg_r = [l for l in roi_ifg if l.startswith('rh.') ]
    index_ifg_r = idx.get_idx(path_direct_folder, res125 + '.txt', roi_ifg_r)
    rifg_merge = 'rh.IFG'
    label_ifg_merged  = combine_labels(roi_ifg_r, labels_pathX, rifg_merge, path_output_ifg)
    idx_ifg_r_merge = [tuple(index_ifg_r)]
    index_y_no_merge  = [(i,) for i in index_y]
    #-DIRECT MATRICES 
    p_xy_rifg, N_xy_rifg, I_xy_rifg, D_med_xy_rifg, D_mean_xy_rifg = matrices.atlas_mat(path_direct_matrix_xy, idx_ifg_r_merge , index_y_no_merge , output_paths_ifg, 'rIFGto33', flag_delays=True)
    p_xx_rifg, N_xx_rifg, I_xx_rifg, D_med_xx_rifg, D_mean_xx_rifg = matrices.atlas_mat(path_direct_matrix_xx,  idx_ifg_r_merge , idx_ifg_r_merge , output_paths_ifg,'rIFGtorIFG', flag_delays=True)
    #-PLOT
    plotting.plot_efferent (N_xy_rifg, N_xx_rifg, [rifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_N, max_cb_N, min_cb_N, 'coolwarm', labels_fig_names = ['RH N Recordings',], border_print= True)
    plotting.plot_efferent(p_xy_rifg, p_xx_rifg, [rifg_merge,], labels_all, labels_pathX, labels_pathY, output_paths_ifg_figs_p, max_cb_p, 0, 'plasma', labels_fig_names=['RH Probability of Direct Connectivity',])
    # plotting.plot_efferent (D_med_xy_rifg, D_med_xx_rifg, [rifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_median, max_cb_dir, min_cb_dir, 'viridis',  labels_fig_names = ['RH Median Peak Delays',])
    plotting.plot_efferent (D_mean_xy_rifg, D_mean_xx_rifg, [rifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_mean, max_cb_dir, min_cb_dir, 'viridis', labels_fig_names = ['RH Mean Peak Delays',])
   
    #INDIRECT
    output_paths_ifg= os.path.join(path_output_ifg, 'Indirect_connectivity_'+ ind_tw)
    os.makedirs(output_paths_ifg,exist_ok=True)
    output_paths_ifg_figs = os.path.join(output_paths_ifg, 'Figures')
    os.makedirs(output_paths_ifg_figs,exist_ok=True)
    output_paths_ifg_figs_p = os.path.join(output_paths_ifg_figs, 'Probability')
    os.makedirs(output_paths_ifg_figs_p,exist_ok=True)
    output_fig_paths_ifg_N = os.path.join(output_paths_ifg_figs, 'N_stims')
    os.makedirs(output_fig_paths_ifg_N,exist_ok=True)
    output_fig_paths_ifg_delays = os.path.join(output_paths_ifg_figs, 'Peak_delays')
    os.makedirs(output_fig_paths_ifg_delays,exist_ok=True)
    # output_fig_paths_ifg_median = os.path.join(output_fig_paths_ifg_delays, 'Median')
    # os.makedirs(output_fig_paths_ifg_median,exist_ok=True)
    output_fig_paths_ifg_mean = os.path.join(output_fig_paths_ifg_delays, 'Mean')
    os.makedirs(output_fig_paths_ifg_mean,exist_ok=True)
    
    #LEFT IFG - INDIRECT MATRICES 
    p_xy_lifg, N_xy_lifg, I_xy_lifg, D_med_xy_lifg, D_mean_xy_lifg = matrices.atlas_mat(path_indirect_matrix_xy, idx_ifg_l_merge , index_y_no_merge , output_paths_ifg, 'lIFGto33', flag_delays=True)
    p_xx_lifg, N_xx_lifg, I_xx_lifg, D_med_xx_lifg, D_mean_xx_lifg = matrices.atlas_mat(path_indirect_matrix_xx, idx_ifg_l_merge , idx_ifg_l_merge , output_paths_ifg, 'lIFGtolIFG', flag_delays=True)
    #-PLOT
    plotting.plot_efferent (N_xy_lifg, N_xx_lifg, [lifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_N, max_cb_N, min_cb_N, 'coolwarm', labels_fig_names = ['LH N Recordings',], border_print= True)
    plotting.plot_efferent(p_xy_lifg, p_xx_lifg, [lifg_merge,], labels_all, labels_pathX, labels_pathY, output_paths_ifg_figs_p, max_cb_p, 0, 'plasma', labels_fig_names=['LH Probability of Indirect Connectivity',])
    # plotting.plot_efferent(D_med_xy_lifg, D_med_xx_lifg, [lifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_median, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names=['LH Median Peak Delays',])
    plotting.plot_efferent(D_mean_xy_lifg, D_mean_xx_lifg, [lifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_mean, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names=['LH Mean Peak Delays',])
    
    #RIGHT IFG - INDIRECT MATRICES
    p_xy_rifg, N_xy_rifg, I_xy_rifg, D_med_xy_rifg, D_mean_xy_rifg = matrices.atlas_mat(path_indirect_matrix_xy, idx_ifg_r_merge , index_y_no_merge , output_paths_ifg, 'rIFGto33', flag_delays=True)
    p_xx_rifg, N_xx_rifg, I_xx_rifg, D_med_xx_rifg, D_mean_xx_rifg = matrices.atlas_mat(path_indirect_matrix_xx, idx_ifg_r_merge , idx_ifg_r_merge , output_paths_ifg, 'rIFGtorIFG', flag_delays=True)
    #-PLOT
    plotting.plot_efferent (N_xy_rifg, N_xx_rifg, [rifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_N, max_cb_N, min_cb_N, 'coolwarm', labels_fig_names = ['RH N Recordings',])
    plotting.plot_efferent(p_xy_rifg, p_xx_rifg, [rifg_merge,], labels_all, labels_pathX, labels_pathY, output_paths_ifg_figs_p, max_cb_p, 0, 'plasma', labels_fig_names=['RH Probability of Indirect Connectivity',])
    # plotting.plot_efferent(D_med_xy_rifg, D_med_xx_rifg, [rifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_median, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names=['RH Median Peak Delays',])
    plotting.plot_efferent(D_mean_xy_rifg, D_mean_xx_rifg, [rifg_merge,], labels_all, labels_pathX, labels_pathY, output_fig_paths_ifg_mean, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names=['RH Mean Peak Delays',])
   
    #%% ROI CONNECTIVITY TO FUNCTIONAL NETWORKS (FIG 4C) 
    # Connectivity from roi (Lausanne x) to functional networks as defined by Yeo et al. We compute the functional networks as mergeds of Lausanne 250 parcels. 
    path_folder_funct = os.path.join(output_folder, 'Functional_nets')
    os.makedirs(path_folder_funct, exist_ok=True)
    path_folder_funct_fig = os.path.join(path_folder_funct, 'Figures')
    os.makedirs(path_folder_funct_fig, exist_ok = True)
    path_folder_funct_fig_brains = os.path.join(path_folder_funct_fig, 'Brain_Plots')
    os.makedirs(path_folder_funct_fig_brains, exist_ok = True)
    path_folder_funct_fig_bars = os.path.join(path_folder_funct_fig, 'Bar_Plots')
    os.makedirs(path_folder_funct_fig_bars, exist_ok = True)
    #folder for object labels for plotting networks 
    path_folder_labels_nets = os.path.join(path_folder_funct, 'Labels_nets')
    os.makedirs(path_folder_labels_nets, exist_ok = True)
    
    [names_labels_regions, index_net, index_roi] = functional_networks_250(output_folder) #(labels_250)
    labels_nets = [[labels_250[idx] for idx in tup] for tup in index_net]
    colors  = ['#89259a', '#64aff3','#2ecd00', '#f566d6', '#ffe496', '#ff972b', '#ff1010'] #purple,blue, green, pink, yellow, orange, red
    roi_color = colors + colors
    #create label objects for networks (left and right)
    for l in range(0,len(names_labels_regions)) :
        combine_labels(labels_nets[l], labels_path250, new_label_name = names_labels_regions[l],  labels_output_path = path_folder_labels_nets )
    #Plot definition of networks (FIG 4A second block)
    plotting.plot_functional_net (output_folder, labels_all_x, labels_pathX, path_folder_labels_nets, roi_color, path_folder_funct_fig)
    
    #LEFT ROI - (FIG 4C first row)
    roi_dlpfc_ifg_l = [l for l in roi_dlpfc_ifg if l.startswith('lh.')]
    index_dlpfc_ifg_l = idx.get_idx(path_general_folder, res125 + '.txt', roi_dlpfc_ifg_l)
    index_roi_L_no_merge = [(i,) for i in index_dlpfc_ifg_l]
    p_net_eff, N_net_eff, I_net_eff = matrices.atlas_mat(path_general_matrix_125_250, index_roi_L_no_merge , index_net[0:int(len(names_labels_regions)/2)] , path_folder_funct, 'Eff_roi_L_nets_' + gral_tw, flag_delays=False)
    plotting.plot_eff(path_folder_funct, p_net_eff, names_labels_regions[0:int(len(names_labels_regions)/2)],roi_dlpfc_ifg_l,path_folder_labels_nets, labels_pathX, res125 + '.txt', path_folder_funct_fig_brains,  max_cb_p, roi_color,  'plasma',  'cvs_avg35_inMNI152 -Lausanne250' )
    
    #RIGHT ROI (supplementary)
    roi_dlpfc_ifg_r = [l for l in roi_dlpfc_ifg if l.startswith('rh.')]
    index_dlpfc_ifg_r = idx.get_idx(path_general_folder, res125 + '.txt', roi_dlpfc_ifg_r)
    index_roi_R_no_merge = [(i,) for i in index_dlpfc_ifg_r]
    p_net_eff, N_net_eff, I_net_eff = matrices.atlas_mat(path_general_matrix_125_250, index_roi_R_no_merge , index_net[int(len(names_labels_regions)/2):len(names_labels_regions)] , path_folder_funct, 'Eff_roi_R_nets_' + gral_tw, flag_delays=False)
    plotting.plot_eff(path_folder_funct, p_net_eff, names_labels_regions[int(len(names_labels_regions)/2):len(names_labels_regions)],roi_dlpfc_ifg_r, path_folder_labels_nets, labels_pathX, res125 + '.txt', path_folder_funct_fig_brains,  max_cb_p, roi_color,  'plasma',  'cvs_avg35_inMNI152 -Lausanne250' )
    
    #bar plots (FIG 4C second row)
    os.chdir(path_folder_funct_fig_bars)
    plot_bars.plot_fine_data(output_folder)
    
    #bar plots rh (supplementary)
    os.chdir(path_folder_funct_fig_bars)
    plot_bars.plot_fine_data_rh(output_folder)
    
    #%% CONNECTIVITY FROM ROI SEGMENTATIONS TO FUNCTIONAL NETWORKS (FIG 4A & 4B):  SUP_DLPFC INF_DLPFC ANT_DLPFC POST_DLPFC IFG
    path_folder_funct_m_roi = os.path.join(path_folder_funct, 'Segmented_Roi')
    os.makedirs(path_folder_funct_m_roi, exist_ok=True)
    path_folder_labels_segments = os.path.join(path_folder_funct_m_roi, 'labels_segments')
    os.makedirs(path_folder_labels_segments, exist_ok=True)
    path_folder_ant_post = os.path.join(path_folder_funct_m_roi, 'ant_post_DLPFC_IFG')
    os.makedirs(path_folder_ant_post, exist_ok=True)
    path_folder_sup_inf = os.path.join(path_folder_funct_m_roi, 'sup_inf_DLPFC_IFG')
    os.makedirs(path_folder_sup_inf, exist_ok=True)
    
    path_folder_ant_post_figs = os.path.join(path_folder_ant_post, 'Figures')
    os.makedirs(path_folder_ant_post_figs, exist_ok=True)
    path_folder_sup_inf_figs = os.path.join(path_folder_sup_inf, 'Figures')
    os.makedirs(path_folder_sup_inf_figs, exist_ok=True)

    path_folder_funct_fig_bars = os.path.join(path_folder_funct_m_roi, 'Bar_Plots')
    os.makedirs(path_folder_funct_fig_bars, exist_ok = True)
    #LEFT
    roi_ant_dlpfc = [
        "lh.rostralmiddlefrontal_1",
        "lh.rostralmiddlefrontal_2",
        "lh.rostralmiddlefrontal_3",]
    roi_post_dlpfc = [
        "lh.caudalmiddlefrontal_1", 
        "lh.caudalmiddlefrontal_2",
        "lh.caudalmiddlefrontal_3",
        "lh.superiorfrontal_6",]
    roi_sup_dlpfc = [
        "lh.superiorfrontal_6",
        "lh.caudalmiddlefrontal_1",
        "lh.caudalmiddlefrontal_3",
        "lh.rostralmiddlefrontal_2",
        ]
    roi_inf_dlpfc = [
        "lh.rostralmiddlefrontal_3",
        "lh.rostralmiddlefrontal_1",
        "lh.caudalmiddlefrontal_2",  
        ]
    roi_ifg_L = [  
        "lh.parsopercularis_1",
        "lh.parsopercularis_2",
        "lh.parstriangularis_1",
                 ]
    idx_ant_dlpfc = idx.get_idx(path_general_folder, res125 + '.txt', roi_ant_dlpfc)
    idx_post_dlpfc = idx.get_idx(path_general_folder, res125 + '.txt', roi_post_dlpfc)
    idx_sup_dlpfc = idx.get_idx(path_general_folder, res125 + '.txt', roi_sup_dlpfc)
    idx_inf_dlpfc = idx.get_idx(path_general_folder, res125 + '.txt', roi_inf_dlpfc)
    idx_ifg_L = idx.get_idx(path_general_folder, res125 + '.txt', roi_ifg_L)
    
    ant_dlpfc = 'lh.ant_DLPFC'
    label_ant_dlpfc  = combine_labels(roi_ant_dlpfc, labels_pathX, ant_dlpfc, path_folder_labels_segments)
    post_dlpfc = 'lh.post_DLPFC'
    label_post_dlpfc  = combine_labels(roi_post_dlpfc, labels_pathX, post_dlpfc, path_folder_labels_segments)
    ifg_l = 'lh.IFG'
    label_l_ifg  = combine_labels(roi_ifg_L, labels_pathX, ifg_l, path_folder_labels_segments)
   
    #Anterior-Posterior + IFG
    p_lifg_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_ifg_L)] , index_net[0:int(len(names_labels_regions)/2)] , path_folder_ant_post, 'l_IFG_nets', flag_delays=False)
    p_ant_dlpfc_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_ant_dlpfc )], index_net[0:int(len(names_labels_regions)/2)] , path_folder_ant_post, 'l_antDLPFC_nets', flag_delays=False)
    p_post_dlpfc_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_post_dlpfc)] , index_net[0:int(len(names_labels_regions)/2)] , path_folder_ant_post, 'l_postDLPFC_nets', flag_delays=False)
    
    p_ant_post = np.vstack((p_ant_dlpfc_net, p_post_dlpfc_net , p_lifg_net)) 
    roi_ant_post = [ant_dlpfc,post_dlpfc,ifg_l]
    plotting.plot_eff(path_folder_funct, p_ant_post, names_labels_regions[0:int(len(names_labels_regions)/2)], roi_ant_post, path_folder_labels_nets, path_folder_labels_segments, res125 + '.txt', path_folder_ant_post_figs,  max_cb_p , roi_color,  'plasma',  'cvs_avg35_inMNI152 -Lausanne250' )
    
    #Superior-Inferior +IFG
    p_sup_dlpfc_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_sup_dlpfc)] , index_net[0:int(len(names_labels_regions)/2)] , path_folder_sup_inf, 'l_supDLPFC_nets', flag_delays=False)
    p_inf_dlpfc_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_inf_dlpfc)] , index_net[0:int(len(names_labels_regions)/2)] , path_folder_sup_inf, 'l_infDLPFC_nets', flag_delays=False)
    p_lifg_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_ifg_L)] , index_net[0:int(len(names_labels_regions)/2)] , path_folder_sup_inf, 'l_IFG_nets', flag_delays=False)
    
    sup_dlpfc = 'lh.sup_DLPFC'
    label_sup_dlpfc  = combine_labels(roi_sup_dlpfc, labels_pathX, sup_dlpfc, path_folder_labels_segments)
    inf_dlpfc = 'lh.inf_DLPFC'
    label_inf_dlpfc  = combine_labels(roi_inf_dlpfc, labels_pathX, inf_dlpfc, path_folder_labels_segments)
    # ifg_l = 'lh.IFG'
    # label_l_ifg  = combine_labels(roi_ifg_L, labels_pathX, ifg_l, path_folder_labels_segments)
    p_sup_inf = np.vstack((p_sup_dlpfc_net, p_inf_dlpfc_net , p_lifg_net)) 
    roi_sup_inf = [sup_dlpfc,inf_dlpfc,ifg_l]
    plotting.plot_eff(path_folder_funct, p_sup_inf, names_labels_regions[0:int(len(names_labels_regions)/2)], roi_sup_inf,  path_folder_labels_nets, path_folder_labels_segments, res125 + '.txt', path_folder_sup_inf_figs,  max_cb_p , roi_color,  'plasma',  'cvs_avg35_inMNI152 -Lausanne250' )
    
    #RIGHT
    roi_ant_dlpfc = [
        "rh.superiorfrontal_4",
        "rh.rostralmiddlefrontal_1",
        "rh.rostralmiddlefrontal_3",]
    roi_post_dlpfc = [
        "rh.caudalmiddlefrontal_1", 
        "rh.caudalmiddlefrontal_2",
        "rh.caudalmiddlefrontal_3",
        "rh.superiorfrontal_8",
        "rh.rostralmiddlefrontal_2",] #check this rostral here
    roi_sup_dlpfc = [
        "rh.superiorfrontal_4",
        "rh.superiorfrontal_8",
        "rh.caudalmiddlefrontal_3",
        "rh.rostralmiddlefrontal_1",
        ]
    roi_inf_dlpfc = [
        "rh.rostralmiddlefrontal_3",
        "rh.rostralmiddlefrontal_2",
        "rh.caudalmiddlefrontal_1",  
        "rh.caudalmiddlefrontal_2",  
        ]
    roi_ifg_R = [  
        "rh.parsopercularis_1",
        "rh.parsopercularis_2",
        "rh.parstriangularis_1",
        "rh.parstriangularis_2",
                 ]
    idx_ant_dlpfc = idx.get_idx(path_general_folder, res125 + '.txt', roi_ant_dlpfc)
    idx_post_dlpfc = idx.get_idx(path_general_folder, res125 + '.txt', roi_post_dlpfc)
    idx_sup_dlpfc = idx.get_idx(path_general_folder, res125 + '.txt', roi_sup_dlpfc)
    idx_inf_dlpfc = idx.get_idx(path_general_folder, res125 + '.txt', roi_inf_dlpfc)
    idx_ifg_L = idx.get_idx(path_general_folder, res125 + '.txt', roi_ifg_R)
    
    ant_dlpfc = 'rh.ant_DLPFC'
    label_ant_dlpfc  = combine_labels(roi_ant_dlpfc, labels_pathX, ant_dlpfc, path_folder_labels_segments)
    post_dlpfc = 'rh.post_DLPFC'
    label_post_dlpfc  = combine_labels(roi_post_dlpfc, labels_pathX, post_dlpfc, path_folder_labels_segments)
    ifg_r = 'rh.IFG'
    label_r_ifg  = combine_labels(roi_ifg_R, labels_pathX, ifg_r, path_folder_labels_segments)
   
    #Anterior-Posterior + IFG
    p_rifg_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_ifg_L)] , index_net[int(len(names_labels_regions)/2):len(names_labels_regions)] , path_folder_ant_post, 'r_IFG_nets', flag_delays=False)
    p_ant_dlpfc_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_ant_dlpfc )], index_net[int(len(names_labels_regions)/2):len(names_labels_regions)] , path_folder_ant_post, 'r_antDLPFC_nets', flag_delays=False)
    p_post_dlpfc_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_post_dlpfc)] , index_net[int(len(names_labels_regions)/2):len(names_labels_regions)] , path_folder_ant_post, 'r_postDLPFC_nets', flag_delays=False)
    
    p_ant_post = np.vstack((p_ant_dlpfc_net, p_post_dlpfc_net , p_rifg_net)) 
    roi_ant_post = [ant_dlpfc,post_dlpfc,ifg_r]
    plotting.plot_efferent_fn(path_folder_funct, p_ant_post, names_labels_regions[int(len(names_labels_regions)/2):len(names_labels_regions)], roi_ant_post, path_folder_labels_nets, path_folder_labels_segments, res125 + '.txt', path_folder_ant_post_figs,  max_cb_p , roi_color,  'plasma',  'cvs_avg35_inMNI152 -Lausanne250' )
    
    #Superior-Inferior +IFG
    p_sup_dlpfc_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_sup_dlpfc)] , index_net[int(len(names_labels_regions)/2):len(names_labels_regions)] , path_folder_sup_inf, 'r_supDLPFC_nets', flag_delays=False)
    p_inf_dlpfc_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_inf_dlpfc)] , index_net[int(len(names_labels_regions)/2):len(names_labels_regions)] , path_folder_sup_inf, 'r_infDLPFC_nets', flag_delays=False)
    p_lifg_net, x,x = matrices.atlas_mat(path_general_matrix_125_250, [tuple(idx_ifg_L)] , index_net[int(len(names_labels_regions)/2):len(names_labels_regions)] , path_folder_sup_inf, 'r_IFG_nets', flag_delays=False)
    
    #im repeting ifg to have it in both folders
    sup_dlpfc = 'rh.sup_DLPFC'
    label_sup_dlpfc  = combine_labels(roi_sup_dlpfc, labels_pathX, sup_dlpfc, path_folder_labels_segments)
    inf_dlpfc = 'rh.inf_DLPFC'
    label_inf_dlpfc  = combine_labels(roi_inf_dlpfc, labels_pathX, inf_dlpfc, path_folder_labels_segments)
    # ifg_r = 'rh.IFG'
    # label_r_ifg  = combine_labels(roi_ifg_R, labels_pathX, ifg_r, path_folder_sup_inf)
    p_sup_inf = np.vstack((p_sup_dlpfc_net, p_inf_dlpfc_net , p_rifg_net)) 
    roi_sup_inf = [sup_dlpfc,inf_dlpfc,ifg_r]
    plotting.plot_efferent_fn(path_folder_funct, p_sup_inf, names_labels_regions[int(len(names_labels_regions)/2):len(names_labels_regions)], roi_sup_inf, path_folder_labels_nets, path_folder_labels_segments, res125 + '.txt', path_folder_sup_inf_figs,  max_cb_p , roi_color,  'plasma',  'cvs_avg35_inMNI152 -Lausanne250' )
    
    #bar plots (FIG4B): ant/post DLPFC + IFG and inf/sup DLPFC + IFG
    os.chdir(path_folder_funct_fig_bars)
    dv = plot_bars.get_coarse_data(output_folder, L_s)
    # print(dv.keys())
    for hemi in ('l', 'r'): #Maybe sep left and right 
        plot_bars.plot_coarse_data(dv, hemi, 'p', 'CI', 0.25) #maybe automat this
        plot_bars.plot_coarse_data(dv, hemi, 'N', None, 20000)
        
    
        
if __name__ == "__main__":
        main()
        
        