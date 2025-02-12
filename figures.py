# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 17:19:52 2025

@author: avalos-alais.s
"""

# Generation of all figures in Avalos-Alais & Jedynak et al 2025.
        
import os
import numpy as np
import matplotlib.pyplot as plt
from tools.combine_labels import combine_labels
from tools.def_funtional_250 import functional_networks_250
from tools import index_tool as idx
from tools import common_lpfc as lpfc
from plotting.scatter_p_delay import analyze_connectivity
from plotting import imshow_brace 
from plotting import plot_scatters as scatter_eff_aff
from plotting import plot_scatters as scatter_p_delay
from plotting import plotting
from plotting import plot_bars

def main():
    #%% DEFINITIONS -  X rows ; Y columns of the matrix
    matrix_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Results') # with saved Lausanne parcellations
    output_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Figures_paper') 
    L_s =  'Lausanne2008'
    res125 = 'Lausanne2008-125'  
    res33  = 'Lausanne2008-33'
    res250 = 'Lausanne2008-250'
    res500 = 'Lausanne2008-500'

    #Data : atlas generation output folders (containing matrices)
    
    gral_tw = '0_100ms'     # time window for general analysis
    dir_tw  = '0_50ms'      # time window of 'direct connections' 
    ind_tw  = '100_400ms'   # time window of 'indirect connections'

    #Path to folder of mne data
    labels_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'MNE-data\subjects\cvs_avg35_inMNI152')
    labels_pathX =  os.path.join(labels_path + '-Lausanne125' , 'label') #for plots
    labels_pathY =  os.path.join(labels_path + '-Lausanne33', 'label'   )  #for plots
    labels_path_500 =  os.path.join(labels_path + '-Lausanne500','label')  #for resolutions
    labels_path250 =  os.path.join(labels_path + '-Lausanne250','label')   #for networks
    
    meshdirname = labels_path #Default in plotting functions, change parameter direction for implementation 
    
    #For functional network analysis
    os.chdir(matrix_folder)
    labels_250 = [line.strip() for line in open( res250 + '.txt' )]
    #def
    labels_500 = [line.strip() for line in open(res500 + '.txt' )]
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
  
    os.chdir(matrix_folder)
    #labels of Y (33) parcellation
    labels_all_y = [line.strip() for line in open(res33 + '.txt')]
    labels_all_y = [s.replace("ctx-lh-", "lh.").replace("ctx-rh-", "rh.") for s in labels_all_y] #special for Lau33
    with open(res33 + '.txt', 'w') as file:
        file.write('\n'.join(labels_all_y) + '\n')
    labels_all  = [name for name in labels_all_y  if name.startswith("lh.") or name.startswith("rh.") or name== 'Left-Hippocampus' or name=='Left-Amygdala'  or name == 'Right-Hippocampus' or name == 'Right-Amygdala']
    #labels of X (125)  parcellation
    labels_all_x = [line.strip() for line in open(res125 + '.txt')]

   
    #delays
    # thresh_d = 20
    min_cb_dir = 20
    max_cb_dir = 35
    min_cb_ind = 150
    max_cb_ind = 200
    #p avg 
    max_avg  = 0.175
    min_avg  = 0.1
    #p
    max_cb_p = 0.25
    max_p_gral = 0.5
    #N 
    max_cb_N = 10000
    min_cb_N = 50
    #%% OUTPUT FOLDERS - general
    path_figs_1 = os.path.join(output_folder, 'Fig1')
    os.makedirs(path_figs_1, exist_ok=True)    
    path_figs_2 = os.path.join(output_folder, 'Fig2')
    os.makedirs(path_figs_2, exist_ok=True)        
    path_folder_output_figs3 = os.path.join(output_folder, 'Fig3')
    os.makedirs(path_folder_output_figs3, exist_ok=True)
    path_folder_output_figs_4 = os.path.join(output_folder, 'Fig4')
    os.makedirs(path_folder_output_figs_4, exist_ok=True)
    path_folder_output_Supplementary = os.path.join(output_folder, 'Supplementary')
    os.makedirs(path_folder_output_Supplementary, exist_ok=True)
    #%% LPFC DEFINITION 
    path_figs_1G = os.path.join(path_figs_1, 'G')
    os.makedirs(path_figs_1G, exist_ok=True)
    plotting.plot_lpfc_definition(path_figs_1G, labels_path, labels_all, meshdirname = meshdirname ) #Fig1G = Fig4A

    #%% PARCELLATION RESOLUTIONS (FIG1D-E-F)
    # Corrected p for square matrices of Lausanne parcellations 33, 125, 500. 
    # Stimulation on selected DLPFC equivalent parcels
    path_resolutions = os.path.join(matrix_folder, 'Resolutions')
    
    path_figs_1_D = os.path.join(path_figs_1, 'D')
    os.makedirs(path_figs_1_D, exist_ok=True)
    path_figs_1_E = os.path.join(path_figs_1, 'E')
    os.makedirs(path_figs_1_E, exist_ok=True)
    path_figs_1_F = os.path.join(path_figs_1, 'F')
    os.makedirs(path_figs_1_F, exist_ok=True)
    
    label_stim =  "lh.rostralmiddlefrontal" #Lau33 stim parcel 
    p_33 = np.loadtxt(path_resolutions + '\\p_' + 'stim33to33_' + gral_tw +'.txt').reshape(1, -1)
    p_stim_33 =  np.loadtxt(path_resolutions + '\\p_' + 'stim33tostim33_' + gral_tw +'.txt').reshape(1, -1)
    p_stim_33 = p_stim_33.reshape(1, -1)
    plotting.plot_efferent (np.array(p_33), np.array(p_stim_33), [label_stim,], labels_all, labels_pathY, labels_pathY, path_figs_1_D, max_p_gral, 0, 'plasma', border_print = True, meshdirname = meshdirname)
    
    label_stim = "lh.rostralmiddlefrontal_1" #Lau125 stim parcel
    p_125 = np.loadtxt(path_resolutions + '\\p_' + 'stim125to125_' + gral_tw +'.txt').reshape(1, -1)
    p_stim_125 = np.loadtxt(path_resolutions + '\\p_' + 'stim125tostim125_' + gral_tw +'.txt').reshape(1, -1)
    plotting.plot_efferent (p_125, p_stim_125, [label_stim,], labels_all_x, labels_pathX, labels_pathX, path_figs_1_E, max_p_gral, 0, 'plasma', border_print = True, meshdirname = meshdirname)
    

    label_stim = "lh.rostralmiddlefrontal_22" #Lau500 stim parcel 
    p_500 = np.loadtxt(path_resolutions + '\\p_' + 'stim500to500_' + gral_tw +'.txt').reshape(1, -1)
    p_stim_500 = np.loadtxt(path_resolutions + '\\p_' + 'stim500tostim500_' + gral_tw +'.txt').reshape(1, -1)
    plotting.plot_efferent (p_500, p_stim_500, [label_stim,], labels_500, labels_path_500, labels_path_500, path_figs_1_F, max_p_gral, 0, 'plasma', border_print = True, meshdirname = meshdirname)
    #%% NUMBER OF IMPLANTED CONTACTS (FIG 2A)
    # Vectors extracted separatly from data. Brain plot of number of implanted contacts for Lausanne33 resolution, with Lausanne125 overlapped over the LPFC roi
    path_N_contacts = os.path.join(matrix_folder, 'N_implanted_contacts')
    os.chdir(path_N_contacts)
    c33 = np.loadtxt('contacts_33.txt').reshape(1, -1)
    c125 = np.loadtxt('contacts_125.txt').reshape(1, -1)

    path_figs_2_A = os.path.join(path_figs_2, 'A')
    os.makedirs(path_figs_2_A, exist_ok=True)
    roi_dlpfc_ifg_l = [l for l in roi_dlpfc_ifg if l.startswith('lh.') ]
    roi_dlpfc_ifg_r = [l for l in roi_dlpfc_ifg if l.startswith('rh.') ]

    plotting.plot_efferent (c33, c125[0:len(roi_dlpfc_ifg_l)], roi_dlpfc_ifg_l, labels_all, labels_pathX, labels_pathY, path_figs_2_A, 3000, 0, 'coolwarm', border_print = True,  labels_fig_names = ['LH Number of Contacts',] )
    plotting.plot_efferent (c33, c125[:,len(roi_dlpfc_ifg_l)::], roi_dlpfc_ifg_r, labels_all, labels_pathX, labels_pathY, path_figs_2_A, 3000, 0, 'coolwarm', border_print = True,  labels_fig_names = ['RH Number of Contacts',] )
    
    #%% AVERAGE CONNECTIVITY (FIG 2B)
    # Using square matrix of LPFC parcellation (X-X). All LPFC merged as a block toward all the rest of the ipsilateral brain (excluding LPFC). 
    # Analysis to prove general efferent-afferent values of the ROI as one. 
    # AVG LPFC parcels (Lausanne125) - Brain (FIG 2B)
    # ROI (LPFC) parcelled according to Lausanne2008-125 towards the rest of the ipsilateral brain (excluding the ROI) merged as one. 
    # Use mat xx, with x parcellatioin of roi (for the exclusion of the roi to be possible)
    path_output_avg = os.path.join(matrix_folder, 'AVG')
    #AVG PLOT [LEFT, RIGHT]
    os.chdir(path_output_avg)
    p_l_roi_all_eff = np.loadtxt('p_avg_l_roi_eff.txt').reshape(-1, 1)
    p_r_roi_all_eff = np.loadtxt('p_avg_r_roi_eff.txt').reshape(-1, 1)
    p_l_roi_all_aff = np.loadtxt('p_avg_l_roi_aff.txt').reshape(1, -1)
    p_r_roi_all_aff = np.loadtxt('p_avg_r_roi_aff.txt').reshape(1, -1)
    avg_eff = np.vstack((p_l_roi_all_eff, p_r_roi_all_eff)) 
    avg_aff = np.hstack((p_l_roi_all_aff, p_r_roi_all_aff)) 
    avg_aff = avg_aff.transpose() 
    
    path_figs_2_B = os.path.join(path_figs_2, 'B')
    os.makedirs(path_figs_2_B, exist_ok=True)
    #Plot Efferent/Afferent average connectivity of roi in each hemisphere
    plotting.plot_mean_p(avg_eff, roi_dlpfc_ifg, labels_pathX, path_figs_2_B, 'AVG_EFF.svg', 'Average Efferent Connectivity',  vmax = max_avg, vmin = min_avg , cbar = 'plasma', nan_color = '#c8c8c8', meshdirname = meshdirname)#light gray 
    plotting.plot_mean_p(avg_aff, roi_dlpfc_ifg, labels_pathX, path_figs_2_B, 'AVG_AFF.svg', 'Average Afferent Connectivity',  vmax = max_avg, vmin = min_avg , cbar = 'plasma' , nan_color = '#c8c8c8', meshdirname = meshdirname)
    #%% EFFERENT CONNECTIVITY (FIG 2C and FIG 2E - first row)  -  CONNECTIVITY FROM X to Y 
    # Including right hemisphere stimulations - supplementary
    # Matrices of combined resolution xy. Fine resolution x over roi and bigger resolution y over all the brain. 
    # From the raw original matrix we take the lines of roi stimulated parcels towards colmuns of interest (in our case all cortical parcels + Amygdala and Hyppocampus)
    # This is in the matrix of parcellations x to y
    # We use xx matrix for stimulations to the roi recorded by the roi.
    
    path_folder_output = os.path.join(matrix_folder, res125 + '__' + res33)
    os.chdir(path_folder_output)
    p_xy = np.loadtxt('p_125to33_0_100ms.txt')
    p_xx = np.loadtxt('p_125to125_0_100ms.txt')
    path_figs_2_C = os.path.join(path_figs_2, 'C')
    os.makedirs(path_figs_2_C, exist_ok=True)
    path_figs_2_C_I= os.path.join(path_figs_2_C, 'I')
    os.makedirs(path_figs_2_C_I, exist_ok=True)
    #plot xy matrix and the xx matrix on top of it : probability
    plotting.plot_efferent (p_xy, p_xx, roi_dlpfc_ifg, labels_all, labels_pathX, labels_pathY, path_figs_2_C_I, max_p_gral, 0, 'plasma', meshdirname = meshdirname)
    
    #%% AFFERENT CONNECTIVITY (FIG 2E and FIG 2C - second row) -  CONNECTIVITY FROM Y to X 
    # Including right hemisphere stimulations - supplementary
    # From the raw original matrix we take the lines of roi recorded parcels of stimulations done on colmuns of interest on the rest of the brain 
    path_folder_output_aff = os.path.join(matrix_folder, res33 + '__' + res125)
    os.chdir(path_folder_output_aff)
    p_yx = np.loadtxt('p_33to125_0_100ms.txt')
    p_xx = np.loadtxt('p_125to125_0_100ms.txt')
    path_figs_2_C_II= os.path.join(path_figs_2_C, 'II')
    os.makedirs(path_figs_2_C_II, exist_ok=True)
    plotting.plot_efferent (p_yx.transpose(), p_xx, roi_dlpfc_ifg, labels_all, labels_pathX, labels_pathY, path_figs_2_C_II,  max_p_gral, 0, 'plasma', meshdirname = meshdirname)
    
    #%% MATRIX PLOT (FIG2E)
    path_figs_2_E = os.path.join(path_figs_2, 'E')
    os.makedirs(path_figs_2_E, exist_ok=True)
    imshow_brace.plot_aff_eff_with_braces(matrix_folder,path_figs_2_E, 'p', 'plasma', None, cb_max = 0.5)
    
    #%% SCATTER PLOT / SYMMETRY OF CONNECTIVITY DIRECTIONALITY 
    path_figs_2_G = os.path.join(path_figs_2, 'G') 
    os.makedirs(path_figs_2_G, exist_ok=True)
    
    diff = 0.1
    min_ = -0.01
    alpha = 0.05
    # read data
    names_125, names_33, p_125_33, N_125_33, CI_125_33, p_33_125, N_33_125, CI_33_125 = lpfc.read_data()
    # eff vs aff 
    asym_mask, p1_p2 = scatter_eff_aff._get_asymmetric_connections(p_125_33, p_33_125.T, N_125_33, N_33_125.T, diff=diff, alpha=alpha, debug=True)
    p1_p2[asym_mask] = np.nan
   
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot()
    regression_q_l, regression_ols_l = scatter_eff_aff._plot_scatters(p_125_33, p_33_125.T, lpfc.colors_lpfc, diff, ax, asym_mask, min_)
    R2_l, a_l, a_stderr_l = [], [], []
    for regression_ols in regression_ols_l:
        a_l.append(regression_ols.params['x'])
        R2_l.append(regression_ols.rsquared)
        a_stderr_l.append(regression_ols.bse['x'])
    print('R2, a, a_stderr for OLS:')
    for name_125, R_2, a, a_stderr in zip(names_125, R2_l, a_l, a_stderr_l):
        print(name_125, R_2, a, a_stderr)
    R2_mean = sum(R2_l)/len(R2_l)
    print('Average R2', R2_mean, 'Std R2', np.std(np.array(R2_l)))
    R2_min = min(R2_l)
    print('Min R2', (R2_min), names_125[R2_l.index(R2_min)])
    fontsize = 26 
    fontname = 'Arial'
    ax.text(0.46, 0.32, r'average $R^2$ = {}'.format(round(R2_mean, 2)), fontsize=fontsize + 4, fontname=fontname)
    ax.spines[['right', 'top']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(4)
    ax.set_xlim([min_, 0.7])
    ax.set_ylim([min_, 0.7])
    ax.set_xlabel('Efferent prob. connectivity', fontsize=fontsize, fontname=fontname)
    ax.set_ylabel('Afferent prob. connectivity', fontsize=fontsize, fontname=fontname)
    plt.tight_layout()
    
    file_path = os.path.join(path_figs_2_G, 'scatter_eff_aff.svg')
    fig.savefig(file_path)
    plt.show()
    #%% MERGED ROIS: DLPFC & IFG. NUMBER of RECORDINGS (FIG 3A) DIRECT (FIG 3B) & INDIRECT (FIG 3C) CONNECTIVITY (probability and delays)
    # Including right hemisphere stimulations - supplementary
    # DLPFC merged (using Lausanne x) towards the rest of the brain in Lausanne y
    # DIRECT 
    path_output_dlpfc = os.path.join(matrix_folder, 'Mean_DLPFC_Eff')
    output_paths_dlpfc = os.path.join(path_output_dlpfc, 'Diect_connectivity_'+ dir_tw )

    figs3_DLPFC = os.path.join(path_folder_output_figs3, 'DLPFC')
    os.makedirs(figs3_DLPFC, exist_ok=True)
    figs3_DLPFC_A = os.path.join(figs3_DLPFC, 'A') # N recordings
    os.makedirs(figs3_DLPFC_A, exist_ok=True)
    figs3_DLPFC_B = os.path.join(figs3_DLPFC, 'B') # Direct p and delays
    os.makedirs(figs3_DLPFC_B, exist_ok=True)
    
    #LEFT DLPFC 
    Ldlpfc_merge = 'lh.DLPFC'
    #- DIRECT MATRICES
    os.chdir(output_paths_dlpfc)
    name_xy = 'lDLPFCto33'
    name_xx = 'lDLPFCtolDLPFC'
    N_xy_ldlpfc = np.loadtxt('N_' + name_xy + '.txt').reshape(1, -1)
    N_xx_ldlpfc = np.loadtxt('N_' + name_xx + '.txt').reshape(1, -1)
    
    p_xy_dir_ldlpfc = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    p_xx_dir_ldlpfc = np.loadtxt('p_' + name_xx + '.txt').reshape(1, -1)
    
    D_mean_xy_dir_ldlpfc = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
    D_mean_xx_dir_ldlpfc = np.loadtxt('peak_delay_mean_' + name_xx + '.txt').reshape(1, -1)

    #   N recordings
    plotting.plot_efferent (N_xy_ldlpfc, N_xx_ldlpfc, [Ldlpfc_merge,], labels_all, path_output_dlpfc, labels_pathY, figs3_DLPFC_A, max_cb_N, min_cb_N, 'coolwarm', border_print = True, labels_fig_names = ['LH N Recordings',], meshdirname = meshdirname)
    #   p
    plotting.plot_efferent (p_xy_dir_ldlpfc, p_xx_dir_ldlpfc, [Ldlpfc_merge,], labels_all, path_output_dlpfc, labels_pathY, figs3_DLPFC_B, max_cb_p, 0, 'plasma', labels_fig_names = ['LH Probability of Direct Connectivity',], meshdirname = meshdirname)
    #   Delays
    plotting.plot_efferent (D_mean_xy_dir_ldlpfc, D_mean_xx_dir_ldlpfc, [Ldlpfc_merge,], labels_all, path_output_dlpfc, labels_pathY, figs3_DLPFC_B, max_cb_dir, min_cb_dir, 'viridis', labels_fig_names = ['LH Mean Peak Delays',], meshdirname = meshdirname)
       
    #RIGHT DLPFC 
    Rdlpfc_merge = 'rh.DLPFC'
    #- DIRECT MATRICES
    os.chdir(output_paths_dlpfc)
    name_xy = 'rDLPFCto33'
    name_xx = 'rDLPFCtorDLPFC'
    N_xy_rdlpfc = np.loadtxt('N_' + name_xy + '.txt').reshape(1, -1)
    N_xx_rdlpfc = np.loadtxt('N_' + name_xx + '.txt').reshape(1, -1)
    
    p_xy_rdlpfc = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    p_xx_rdlpfc = np.loadtxt('p_' + name_xx + '.txt').reshape(1, -1)
    
    D_mean_xy_rdlpfc = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
    D_mean_xx_rdlpfc = np.loadtxt('peak_delay_mean_' + name_xx + '.txt').reshape(1, -1)
    
    #   N recordings
    plotting.plot_efferent (N_xy_rdlpfc, N_xx_rdlpfc, [Rdlpfc_merge,], labels_all, path_output_dlpfc, labels_pathY, figs3_DLPFC_A, max_cb_N, min_cb_N, 'coolwarm', border_print = True ,labels_fig_names = ['RH N Recordings',], meshdirname = meshdirname)
    #   p 
    plotting.plot_efferent (p_xy_rdlpfc, p_xx_rdlpfc, [Rdlpfc_merge,], labels_all, path_output_dlpfc, labels_pathY, figs3_DLPFC_B, max_cb_p, 0, 'plasma', labels_fig_names = ['RH Probability of Direct Connectivity',], meshdirname = meshdirname)
    #   Delays
    plotting.plot_efferent (D_mean_xy_rdlpfc, D_mean_xx_rdlpfc, [Rdlpfc_merge,], labels_all, path_output_dlpfc, labels_pathY, figs3_DLPFC_B, max_cb_dir, min_cb_dir, 'viridis', labels_fig_names = ['RH Mean Peak Delays',], meshdirname = meshdirname)
    
    #INDIRECT
    output_paths_dlpfc = os.path.join(path_output_dlpfc, 'Indirect_connectivity_' + ind_tw)
    figs3_DLPFC_C = os.path.join(figs3_DLPFC, 'C') # Indirect p and delays
    os.makedirs(figs3_DLPFC_C, exist_ok=True)

    #LEFT DLPFC - INDIRECT MATRICES 
    os.chdir(output_paths_dlpfc)
    name_xy = 'lDLPFCto33'
    name_xx = 'lDLPFCtolDLPFC'
    
    p_xy_ind_ldlpfc = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    p_xx_ind_ldlpfc = np.loadtxt('p_' + name_xx + '.txt').reshape(1, -1)
    
    D_mean_xy_ind_ldlpfc = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
    D_mean_xx_ind_ldlpfc = np.loadtxt('peak_delay_mean_' + name_xx + '.txt').reshape(1, -1)

    #-PLOT
    #   p
    plotting.plot_efferent (p_xy_ind_ldlpfc, p_xx_ind_ldlpfc, [Ldlpfc_merge,], labels_all, path_output_dlpfc, labels_pathY, figs3_DLPFC_C, max_cb_p , 0, 'plasma',  labels_fig_names = ['LH Probability of Indirect Connectivity',] , meshdirname = meshdirname)
    #   Delays
    plotting.plot_efferent (D_mean_xy_ind_ldlpfc, D_mean_xx_ind_ldlpfc, [Ldlpfc_merge,], labels_all, path_output_dlpfc, labels_pathY, figs3_DLPFC_C, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names = ['LH Mean Peak Delays',], meshdirname = meshdirname)
    
    #RIGHT DLPFC - INDIRECT MATRICES 
    os.chdir(output_paths_dlpfc)
    name_xy = 'rDLPFCto33'
    name_xx = 'rDLPFCtorDLPFC'
    p_xy_rdlpfc = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    p_xx_rdlpfc = np.loadtxt('p_' + name_xx + '.txt').reshape(1, -1)
    
    D_mean_xy_rdlpfc = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
    D_mean_xx_rdlpfc = np.loadtxt('peak_delay_mean_' + name_xx + '.txt').reshape(1, -1)
    #   p
    plotting.plot_efferent (p_xy_rdlpfc, p_xx_rdlpfc, [Rdlpfc_merge,], labels_all, path_output_dlpfc, labels_pathY, figs3_DLPFC_C, max_cb_p , 0, 'plasma',  labels_fig_names = ['RH Probability of Indirect Connectivity',], meshdirname = meshdirname)
    #   Delays
    plotting.plot_efferent (D_mean_xy_rdlpfc, D_mean_xx_rdlpfc, [Rdlpfc_merge,], labels_all, path_output_dlpfc, labels_pathY, figs3_DLPFC_C, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names = ['RH Mean Peak Delays',], meshdirname = meshdirname)
    
    #IFG merged
    #DIRECT
    path_output_ifg = os.path.join(matrix_folder, 'Mean_IFG_Eff')
    output_paths_ifg= os.path.join(path_output_ifg, 'Diect_connectivity_'+ dir_tw)
    
    figs3_IFG = os.path.join(path_folder_output_figs3, 'IFG')
    os.makedirs(figs3_IFG, exist_ok=True)
    figs3_IFG_A = os.path.join(figs3_IFG, 'A') # N recordings
    os.makedirs(figs3_IFG_A, exist_ok=True)
    figs3_IFG_B = os.path.join(figs3_IFG, 'B') # Direct p and delays
    os.makedirs(figs3_IFG_B, exist_ok=True)
    
    #LEFT IFG
    lifg_merge = 'lh.IFG'

    #-DIRECT MATRICES 
    os.chdir(output_paths_ifg)
    name_xy = 'lIFGto33'
    name_xx = 'lIFGtolIFG'
    N_xy_lifg = np.loadtxt('N_' + name_xy + '.txt').reshape(1, -1)
    N_xx_lifg = np.loadtxt('N_' + name_xx + '.txt').reshape(1, -1)
    
    p_xy_dir_lifg = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    p_xx_dir_lifg = np.loadtxt('p_' + name_xx + '.txt').reshape(1, -1)
    
    D_mean_xy_dir_lifg = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
    D_mean_xx_dir_lifg = np.loadtxt('peak_delay_mean_' + name_xx + '.txt').reshape(1, -1)
    

    plotting.plot_efferent (N_xy_lifg, N_xx_lifg, [lifg_merge,], labels_all, path_output_ifg, labels_pathY, figs3_IFG_A, max_cb_N, min_cb_N, 'coolwarm', labels_fig_names = ['LH N Recordings',], border_print= True, meshdirname = meshdirname)
    plotting.plot_efferent(p_xy_dir_lifg, p_xx_dir_lifg, [lifg_merge,], labels_all, path_output_ifg, labels_pathY, figs3_IFG_B, max_cb_p, 0, 'plasma', labels_fig_names=['LH Probability of Direct Connectivity',], meshdirname = meshdirname)
    plotting.plot_efferent (D_mean_xy_dir_lifg, D_mean_xx_dir_lifg, [lifg_merge,], labels_all, path_output_ifg, labels_pathY, figs3_IFG_B, max_cb_dir, min_cb_dir, 'viridis', labels_fig_names = ['LH Mean Peak Delays',], meshdirname = meshdirname)
    
    #RIGHT IFG
    rifg_merge = 'rh.IFG'
    #-DIRECT MATRICES 
    os.chdir(output_paths_ifg)
    name_xy = 'rIFGto33'
    name_xx = 'rIFGtorIFG'
    N_xy_rifg = np.loadtxt('N_' + name_xy + '.txt').reshape(1, -1)
    N_xx_rifg = np.loadtxt('N_' + name_xx + '.txt').reshape(1, -1)
    
    p_xy_rifg = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    p_xx_rifg = np.loadtxt('p_' + name_xx + '.txt').reshape(1, -1)
    
    D_mean_xy_rifg = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
    D_mean_xx_rifg = np.loadtxt('peak_delay_mean_' + name_xx + '.txt').reshape(1, -1)
    
    plotting.plot_efferent (N_xy_rifg, N_xx_rifg, [rifg_merge,], labels_all, path_output_ifg, labels_pathY, figs3_IFG_A, max_cb_N, min_cb_N, 'coolwarm', labels_fig_names = ['RH N Recordings',], border_print= True, meshdirname = meshdirname)
    plotting.plot_efferent(p_xy_rifg, p_xx_rifg, [rifg_merge,], labels_all, path_output_ifg, labels_pathY, figs3_IFG_B, max_cb_p, 0, 'plasma', labels_fig_names=['RH Probability of Direct Connectivity',], meshdirname = meshdirname)
    plotting.plot_efferent (D_mean_xy_rifg, D_mean_xx_rifg, [rifg_merge,], labels_all, path_output_ifg, labels_pathY, figs3_IFG_B, max_cb_dir, min_cb_dir, 'viridis', labels_fig_names = ['RH Mean Peak Delays',], meshdirname = meshdirname)
    
    #INDIRECT
    output_paths_ifg= os.path.join(path_output_ifg, 'Indirect_connectivity_'+ ind_tw)
    figs3_IFG_C = os.path.join(figs3_IFG, 'C') # Indirect p and delays
    os.makedirs(figs3_IFG_C, exist_ok=True)
    
    #LEFT IFG - INDIRECT MATRICES 
    os.chdir(output_paths_ifg)
    name_xy = 'lIFGto33'
    name_xx = 'lIFGtolIFG'
    
    p_xy_ind_lifg = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    p_xx_ind_lifg = np.loadtxt('p_' + name_xx + '.txt').reshape(1, -1)
    
    D_mean_xy_ind_lifg = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
    D_mean_xx_ind_lifg = np.loadtxt('peak_delay_mean_' + name_xx + '.txt').reshape(1, -1)

    plotting.plot_efferent(p_xy_ind_lifg, p_xx_ind_lifg, [lifg_merge,], labels_all, path_output_ifg, labels_pathY, figs3_IFG_C, max_cb_p, 0, 'plasma', labels_fig_names=['LH Probability of Indirect Connectivity',], meshdirname = meshdirname)
    plotting.plot_efferent(D_mean_xy_ind_lifg, D_mean_xx_ind_lifg, [lifg_merge,], labels_all, path_output_ifg, labels_pathY, figs3_IFG_C, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names=['LH Mean Peak Delays',], meshdirname = meshdirname)
    
    #RIGHT IFG - INDIRECT MATRICES
    os.chdir(output_paths_ifg)
    name_xy = 'rIFGto33'
    name_xx = 'rIFGtorIFG'
    
    p_xy_rifg = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    p_xx_rifg = np.loadtxt('p_' + name_xx + '.txt').reshape(1, -1)
    
    D_mean_xy_rifg = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
    D_mean_xx_rifg = np.loadtxt('peak_delay_mean_' + name_xx + '.txt').reshape(1, -1)
    
    plotting.plot_efferent(p_xy_rifg, p_xx_rifg, [rifg_merge,], labels_all, path_output_ifg, labels_pathY, figs3_IFG_C, max_cb_p, 0, 'plasma', labels_fig_names=['RH Probability of Indirect Connectivity',], meshdirname = meshdirname)
    plotting.plot_efferent(D_mean_xy_rifg, D_mean_xx_rifg, [rifg_merge,], labels_all, path_output_ifg, labels_pathY, figs3_IFG_C, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names=['RH Mean Peak Delays',], meshdirname = meshdirname)  
    
    # Supplementary FIG4 : Probability of connectivity vs mean peak delay
    #Analysis of DLPFC and IFG direct and indirect connections to the rest of the brain 
    analyze_connectivity(p_xy_dir_ldlpfc, D_mean_xy_dir_ldlpfc, path_folder_output_Supplementary, 'Direct connectivity, p vs mean peak delay', "dir_DLPFC.svg", 'r', (20,40), (0,0.3), 'DLPFC')
    analyze_connectivity(p_xy_ind_ldlpfc, D_mean_xy_ind_ldlpfc, path_folder_output_Supplementary, 'Indirect connectivity, p vs mean peak delay', "ind_DLPFC.svg", 'g', (135,225), (0,0.3), 'DLPFC')
    analyze_connectivity(p_xy_dir_lifg, D_mean_xy_dir_lifg, path_folder_output_Supplementary, 'Direct connectivity, p vs mean peak delay', "dir_IFG.svg", 'b', (20,40), (0,0.3), 'IFG')
    analyze_connectivity(p_xy_ind_lifg, D_mean_xy_ind_lifg, path_folder_output_Supplementary, 'Indirect connectivity, p vs mean peak delay', "ind_IFG.svg", 'm', (135,225), (0,0.3), 'IFG')

    #%% ROI CONNECTIVITY TO FUNCTIONAL NETWORKS (FIG 4ABC) 
    # Connectivity from roi (Lausanne x) to functional networks as defined by Yeo et al. We compute the functional networks as mergeds of Lausanne 250 parcels. 
    path_folder_funct = os.path.join(matrix_folder, 'Functional_nets')
    #folder for object labels for plotting networks 
    path_folder_labels_nets = os.path.join(path_folder_funct, 'Labels_nets')

    path_figs_4A = os.path.join(path_folder_output_figs_4, 'A')
    os.makedirs(path_figs_4A, exist_ok=True)
    path_figs_4B = os.path.join(path_folder_output_figs_4, 'B')
    os.makedirs(path_figs_4B, exist_ok=True)
    path_figs_4C = os.path.join(path_folder_output_figs_4, 'C')
    os.makedirs(path_figs_4C, exist_ok=True) 
    path_figs_4C_I = os.path.join(path_figs_4C, 'I')
    os.makedirs(path_figs_4C_I, exist_ok=True)
    path_figs_4C_II = os.path.join(path_figs_4C, 'II')
    os.makedirs(path_figs_4C_II, exist_ok=True)
    
    plotting.plot_lpfc_definition(path_figs_4A, labels_path, labels_all ) #Fig4A = Fig1G
    
    [names_labels_regions, index_net, index_roi] = functional_networks_250(matrix_folder) #labels_250)
    labels_nets = [[labels_250[idx] for idx in tup] for tup in index_net]
    colors  = ['#89259a', '#64aff3','#2ecd00', '#f566d6', '#ffe496', '#ff972b', '#ff1010'] #purple,blue, green, pink, yellow, orange, red
    roi_color = colors + colors
    #create label objects for networks (left and right)
    for l in range(0,len(names_labels_regions)) :
        combine_labels(labels_nets[l], labels_path250, new_label_name = names_labels_regions[l],  labels_output_path = path_folder_labels_nets )
    #Plot definition of networks (FIG 4A second block)
    plotting.plot_functional_net (matrix_folder, labels_all_x, labels_pathX, path_folder_labels_nets, roi_color, path_figs_4A)
    
    #LEFT ROI - (FIG 4C first row)
    roi_dlpfc_ifg_l = [l for l in roi_dlpfc_ifg if l.startswith('lh.')]
    os.chdir(path_folder_funct)
    p_net_eff = np.loadtxt('p_Eff_roi_L_nets_0_100ms.txt')
    plotting.plot_eff(path_folder_funct, p_net_eff, names_labels_regions[0:int(len(names_labels_regions)/2)],roi_dlpfc_ifg_l, path_folder_labels_nets, labels_pathX, res125 + '.txt', path_figs_4C_I,  max_cb_p, roi_color,  'plasma',  'cvs_avg35_inMNI152 -Lausanne250', meshdirname = meshdirname)
    
    #RIGHT ROI (supplementary)
    roi_dlpfc_ifg_r = [l for l in roi_dlpfc_ifg if l.startswith('rh.')]
    os.chdir(path_folder_funct)
    p_net_eff = np.loadtxt('p_Eff_roi_R_nets_0_100ms.txt')
    plotting.plot_eff(path_folder_funct, p_net_eff, names_labels_regions[int(len(names_labels_regions)/2):len(names_labels_regions)],roi_dlpfc_ifg_r, path_folder_labels_nets, labels_pathX, res125 + '.txt', path_figs_4C_I,  max_cb_p, roi_color,  'plasma',  'cvs_avg35_inMNI152 -Lausanne250', meshdirname = meshdirname)
    
    #bar plots (FIG 4C second row)
    os.chdir(path_figs_4C_II)
    plot_bars.plot_fine_data(matrix_folder)
    
    #bar plots rh (supplementary)
    os.chdir(path_figs_4C_II)
    plot_bars.plot_fine_data_rh(matrix_folder)
    #TODO: plot bars need generalization to one ft for both hemis, now dependent on nm of parcels an hardcoded path 
    #bar plots (FIG4B): ant/post DLPFC + IFG and inf/sup DLPFC + IFG
    os.chdir(path_figs_4B)
    dv = plot_bars.get_coarse_data(matrix_folder, L_s)
    # print(dv.keys())
    for hemi in ('l', 'r'): #Maybe sep left and right 
        plot_bars.plot_coarse_data(dv, hemi, 'p', 'CI', 0.25) #maybe automat this
        plot_bars.plot_coarse_data(dv, hemi, 'N', None, 20000)    
        

#%%
if __name__ == "__main__":
        main()