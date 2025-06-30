# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 17:19:52 2025
@author: avalos-alais.s
"""

# Generation of all figures in Avalos-Alais & Jedynak et al 2025.
import os
import numpy as np
import matplotlib.pyplot as plt
import Definitions
from tools.combine_labels import combine_labels
from tools import common_lpfc as lpfc
from plotting.scatter_p_delay import analyze_connectivity
from plotting import imshow_brace 
from plotting import plot_scatters as scatter_eff_aff
from plotting import plotting
from plotting import plot_bars



#TODO : Segment into functions - sep defs
def main():

    #%% DEFINITIONS -  X rows ; Y columns of the matrix
    matrix_folder = r'F:\FTRACT\Data_LPFC_FTRACT\Results_Zth5'
    output_folder = os.path.join(matrix_folder, 'Figures')

    L_s = Definitions.L_s
    res33  = Definitions.res33
    res125 = Definitions.res125
    res250 = Definitions.res250  # FN analysis
    res500 = Definitions.res500  # For resolution analysis

    #Parcellation labels
    labels_33   = Definitions.labels_33
    labels_125  = Definitions.labels_125
    labels_250  = Definitions.labels_250
    labels_500  = Definitions.labels_500
    labels_rec_33 = [name for name in labels_33 if name.startswith("lh.") or name.startswith("rh.") or name == 'Left-Hippocampus' or name == 'Left-Amygdala' or name == 'Right-Hippocampus' or name == 'Right-Amygdala']

    #For plotting
    labels_path = Definitions.labels_path
    labels_path_33  = Definitions.labels_path_33
    labels_path_125 = Definitions.labels_path_125
    labels_path_250  = Definitions.labels_path_250
    labels_path_500 = Definitions.labels_path_500

    #Data : atlas generation output folders (containing matrices)
    gral_tw = '0_100_ms'     # time window for general analysis
    dir_tw  = '0_50_ms'      # time window of 'direct connections'
    ind_tw  = '100_400_ms'   # time window of 'indirect connections'

    #delays
    # thresh_d = 20
    min_cb_dir = 20
    max_cb_dir = 37
    min_cb_ind = 130
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

    # 1
    path_figs_1 = os.path.join(output_folder, 'Fig1')
    os.makedirs(path_figs_1, exist_ok=True)

    #%% LPFC DEFINITION
    path_figs_1I = os.path.join(path_figs_1, 'I')
    os.makedirs(path_figs_1I, exist_ok=True)
    plotting.plot_lpfc_definition(path_figs_1I, Definitions.labels_path, labels_33, meshdirname = Definitions.meshdirname ) #Fig1I = Fig4A

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
    plotting.plot_efferent (np.array(p_33), np.array(p_stim_33), [label_stim,], labels_rec_33, Definitions.labels_path_33, Definitions.labels_path_33, path_figs_1_D, max_p_gral, 0, 'plasma', border_print = True, meshdirname = Definitions.meshdirname)

    label_stim = "lh.rostralmiddlefrontal_1" #Lau125 stim parcel
    p_125 = np.loadtxt(path_resolutions + '\\p_' + 'stim125to125_' + gral_tw +'.txt').reshape(1, -1)
    p_stim_125 = np.loadtxt(path_resolutions + '\\p_' + 'stim125tostim125_' + gral_tw +'.txt').reshape(1, -1)
    plotting.plot_efferent (p_125, p_stim_125, [label_stim,], labels_125, Definitions.labels_path_125, Definitions.labels_path_125, path_figs_1_E, max_p_gral, 0, 'plasma', border_print = True, meshdirname = Definitions.meshdirname)


    label_stim = "lh.rostralmiddlefrontal_22" #Lau500 stim parcel
    p_500 = np.loadtxt(path_resolutions + '\\p_' + 'stim500to500_' + gral_tw +'.txt').reshape(1, -1)
    p_stim_500 = np.loadtxt(path_resolutions + '\\p_' + 'stim500tostim500_' + gral_tw +'.txt').reshape(1, -1)
    plotting.plot_efferent (p_500, p_stim_500, [label_stim,], labels_500, Definitions.labels_path_500, Definitions.labels_path_500, path_figs_1_F, max_p_gral, 0, 'plasma', border_print = True, meshdirname = Definitions.meshdirname)

    # 2
    path_figs_2 = os.path.join(output_folder, 'Fig2')
    os.makedirs(path_figs_2, exist_ok=True)
    #%% NUMBER OF IMPLANTED CONTACTS (FIG 2A)
    # Vectors extracted separatly from data. Brain plot of number of implanted contacts for Lausanne33 resolution, with Lausanne125 overlapped over the LPFC roi
    path_N_contacts = os.path.join(matrix_folder, 'N_implanted_contacts')
    os.chdir(path_N_contacts)
    c33 = np.loadtxt('contacts_33.txt').reshape(1, -1)
    c125 = np.loadtxt('contacts_125.txt').reshape(1, -1)

    path_figs_2_A = os.path.join(path_figs_2, 'A')
    os.makedirs(path_figs_2_A, exist_ok=True)
    roi_dlpfc_ifg_l = [l for l in Definitions.roi_dlpfc_ifg_125 if l.startswith('lh.')]
    roi_dlpfc_ifg_r = [l for l in Definitions.roi_dlpfc_ifg_125 if l.startswith('rh.')]

    plotting.plot_efferent (c33, c125[0:len(roi_dlpfc_ifg_l)], roi_dlpfc_ifg_l, labels_rec_33, Definitions.labels_path_125, Definitions.labels_path_33, path_figs_2_A, 3000, 0, 'coolwarm', border_print = True,  labels_fig_names = ['LH Number of Contacts',] )
    plotting.plot_efferent (c33, c125[:,len(roi_dlpfc_ifg_l)::], roi_dlpfc_ifg_r, labels_rec_33, Definitions.labels_path_125, Definitions.labels_path_33, path_figs_2_A, 3000, 0, 'coolwarm', border_print = True,  labels_fig_names = ['RH Number of Contacts',] )

    #%% AVERAGE CONNECTIVITY (FIG 2B)
        # Using square matrix of LPFC parcellation (X-X). All LPFC merged as a block toward all the rest of the ipsilateral brain (excluding LPFC).
        # Analysis to prove general efferent-afferent values of the ROI as one.
        # AVG LPFC parcels (Lausanne125) - Brain (FIG 2B)
        # ROI (LPFC) parcelled according to Lausanne2008-125 towards the rest of the ipsilateral brain (excluding the ROI) merged as one.
        # Use mat xx, with x parcellatioin of roi (for the exclusion of the roi to be possible)
    path_output_avg= os.path.join(matrix_folder, 'AVG')
    #AVG PLOT [LEFT, RIGHT]
    os.chdir(path_output_avg)
    p_l_roi_all_eff = np.loadtxt('p_lLPFC_brain_all_eff.txt').reshape(-1, 1)
    p_r_roi_all_eff = np.loadtxt('p_rLPFC_brain_all_eff.txt').reshape(-1, 1)
    p_l_roi_all_aff = np.loadtxt('p_lLPFC_brain_all_aff.txt').reshape(1, -1)
    p_r_roi_all_aff = np.loadtxt('p_rLPFC_brain_all_aff.txt').reshape(1, -1)
    avg_eff = np.vstack((p_l_roi_all_eff, p_r_roi_all_eff))
    avg_aff = np.hstack((p_l_roi_all_aff, p_r_roi_all_aff))
    avg_aff = avg_aff.transpose()
    path_figs_2_B = os.path.join(path_figs_2, 'B')
    os.makedirs(path_figs_2_B, exist_ok=True)
    #Plot Efferent/Afferent average connectivity of roi in each hemisphere
    plotting.plot_mean_p(avg_eff, Definitions.roi_dlpfc_ifg_125, Definitions.labels_path_125, path_figs_2_B, 'AVG_EFF.svg', 'Average Efferent Connectivity', vmax = max_avg, vmin = min_avg, cbar ='plasma', nan_color ='#c8c8c8', meshdirname = Definitions.meshdirname)#light gray
    plotting.plot_mean_p(avg_aff, Definitions.roi_dlpfc_ifg_125, Definitions.labels_path_125, path_figs_2_B, 'AVG_AFF.svg', 'Average Afferent Connectivity', vmax = max_avg, vmin = min_avg, cbar ='plasma', nan_color ='#c8c8c8', meshdirname = Definitions.meshdirname)

    # %% EFFERENT CONNECTIVITY (FIG 2C and FIG 2E - first row)  -  CONNECTIVITY FROM X to Y
    #     Including right hemisphere stimulations - supplementary
    #     Matrices of combined resolution xy. Fine resolution x over roi and bigger resolution y over all the brain.
    #     From the raw original matrix we take the lines of roi stimulated parcels towards colmuns of interest (in our case all cortical parcels + Amygdala and Hyppocampus)
    #     This is in the matrix of parcellations x to y
    #     We use xx matrix for stimulations to the roi recorded by the roi.

    path_folder_output = os.path.join(matrix_folder, res125 + '_' + res33)
    os.chdir(path_folder_output)
    # #TODO : automat this names with parameters.
    p_xy = np.loadtxt('p_125to33_0_100_ms.txt')
    p_xx = np.loadtxt('p_125to125_0_100_ms.txt')
    index_y = np.loadtxt('index_y_125to33_0_100_ms.txt', dtype= int)
    path_figs_2_C = os.path.join(path_figs_2, 'C')
    os.makedirs(path_figs_2_C, exist_ok=True)
    path_figs_2_C_I= os.path.join(path_figs_2_C, 'I')
    os.makedirs(path_figs_2_C_I, exist_ok=True)
    # plot xy matrix and the xx matrix on top of it : probability
    plotting.plot_efferent (p_xy, p_xx, Definitions.roi_dlpfc_ifg_125,  [labels_33[idx] for idx in index_y], Definitions.labels_path_125, Definitions.labels_path_33, path_figs_2_C_I, max_p_gral, 0, 'plasma', meshdirname = Definitions.meshdirname)

    # %% AFFERENT CONNECTIVITY (FIG 2E and FIG 2C - second row) -  CONNECTIVITY FROM Y to X
    # Including right hemisphere stimulations - supplementary
    # From the raw original matrix we take the lines of roi recorded parcels of stimulations done on colmuns of interest on the rest of the brain
    path_folder_output_aff = os.path.join(matrix_folder, res33 + '_' + res125)
    os.chdir(path_folder_output_aff)
    p_yx = np.loadtxt('p_33to125_0_100_ms.txt')
    p_xx = np.loadtxt('p_125to125_0_100_ms.txt')
    index_x = np.loadtxt('index_x_33to125_0_100_ms.txt', dtype= int)
    path_figs_2_C_II= os.path.join(path_figs_2_C, 'II')
    os.makedirs(path_figs_2_C_II, exist_ok=True)
    plotting.plot_efferent (p_yx.transpose(), p_xx, Definitions.roi_dlpfc_ifg_125, [labels_33[idx] for idx in index_x], Definitions.labels_path_125, Definitions.labels_path_33, path_figs_2_C_II, max_p_gral, 0, 'plasma', meshdirname = Definitions.meshdirname)

    # %% MATRIX PLOT (FIG2E)
    path_figs_2_E = os.path.join(path_figs_2, 'E')
    os.makedirs(path_figs_2_E, exist_ok=True)
    imshow_brace.plot_aff_eff_with_braces(matrix_folder,path_figs_2_E, 'p', 'plasma', None, cb_max = max_p_gral)

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
    #
    file_path = os.path.join(path_figs_2_G, 'scatter_eff_aff.svg')
    fig.savefig(file_path)
    plt.show()

    # 3
    path_folder_output_figs3 = os.path.join(output_folder, 'Fig3')
    os.makedirs(path_folder_output_figs3, exist_ok=True)
    #%% MERGED ROIS: DLPFC & IFG. NUMBER of RECORDINGS (FIG 3A) DIRECT (FIG 3B) & INDIRECT (FIG 3C) CONNECTIVITY (probability and delays)
    # Including right hemisphere stimulations - supplementary
    path_output_roi_eff = os.path.join(matrix_folder, 'Mean_ROI_Eff')
    # Create labels for merged rois
    #left
    Ldlpfc_merge = 'lh.DLPFC'
    labels_dlpfc_l = [labels_125[i] for i in Definitions.index_dlpfc_L]
    label_dlpfc_merged = combine_labels(labels_dlpfc_l, Definitions.labels_path_125, Ldlpfc_merge, path_output_roi_eff)
    Lifg_merge = 'lh.IFG'
    labels_ifg_l = [labels_125[i] for i in Definitions.index_ifg_L]
    label_ifg_merged = combine_labels(labels_ifg_l, Definitions.labels_path_125, Lifg_merge, path_output_roi_eff)
    # right - this goes to supplementary
    Rdlpfc_merge = 'rh.DLPFC'
    labels_dlpfc_r = [labels_125[i] for i in Definitions.index_dlpfc_R]
    label_dlpfc_merged = combine_labels(labels_dlpfc_r, Definitions.labels_path_125, Rdlpfc_merge, path_output_roi_eff)
    Rifg_merge = 'rh.IFG'
    labels_ifg_r = [labels_125[i] for i in Definitions.index_ifg_R]
    label_ifg_merged = combine_labels(labels_ifg_r, Definitions.labels_path_125, Rifg_merge, path_output_roi_eff)

    figs3_A = os.path.join(path_folder_output_figs3, 'A') # N recordings
    os.makedirs(figs3_A, exist_ok=True)
    figs3_B = os.path.join(path_folder_output_figs3, 'B') # Direct p and delays
    os.makedirs(figs3_B, exist_ok=True)
    figs3_C = os.path.join(path_folder_output_figs3, 'C')  # Indirect p and delays
    os.makedirs(figs3_C, exist_ok=True)

    # #- DIRECT MATRICES
    output_paths_roi_eff_dir = os.path.join(path_output_roi_eff, 'Direct_connectivity_'+ dir_tw )
    os.chdir(output_paths_roi_eff_dir)
    name_xy = 'roi_to_33'
    name_xx = 'roi_to_roi'
    N_xy = np.loadtxt('N_' + name_xy + '.txt')
    p_xy_dir = np.loadtxt('p_' + name_xy + '.txt')
    D_mean_xy_dir = np.loadtxt('peak_delay_mean_' + name_xy + '.txt')
    D_median_xy_dir = np.loadtxt('peak_delay_med_' + name_xy + '.txt')

    #local (merged-merged)
    N_xx = np.loadtxt('N_' + name_xx + '.txt')
    p_xx_dir = np.loadtxt('p_' + name_xx + '.txt')
    D_mean_xx_dir = np.loadtxt('peak_delay_mean_' + name_xx + '.txt')
    D_median_xx_dir = np.loadtxt('peak_delay_med_' + name_xx + '.txt')

    #Diagonal to have just the stim one recording
    N_xx_diag = np.diag(N_xx).reshape(-1, 1)
    p_xx_dir_diag = np.diag(p_xx_dir).reshape(-1, 1)
    D_mean_dir_diag = np.diag(D_mean_xx_dir).reshape(-1, 1)
    D_median_dir_diag = np.diag(D_median_xx_dir).reshape(-1, 1)

    index_rec = np.loadtxt('index_y_roi_to_33.txt', dtype=int)
    labels_roi = [Ldlpfc_merge, Rdlpfc_merge, Lifg_merge, Rifg_merge] #this should be automatized, I know it is in this order but to avoid error it shouldnt depend on that

    figs_names_N  = ['LH DLPFC N Recordings','RH DLPFC N Recordings','LH IFG N Recordings','RH IFG N Recordings']
    figs_p_names = ['LH DLPFC Probability of Direct Connectivity', 'RH DLPFC Probability of Direct Connectivity', 'LH IFG Probability of Direct Connectivity', 'RH IFG Probability of Direct Connectivity', ]
    figs_mean_delays_names = ['LH DLPFC Mean Peak Delays Direct', 'RH DLPFC Mean Peak Delays Direct', 'LH IFG Mean Peak Delays Direct', 'RH IFG Mean Peak Delays Direct']
    figs_median_delays_names = ['LH DLPFC Median Peak Delays Direct', 'RH DLPFC Median Peak Delays Direct','LH IFG Median Peak Delays Direct', 'RH IFG Median Peak Delays Direct']

    for i in range(len(figs_names_N)):
        #   N recordings
        plotting.plot_efferent (N_xy[i,:].reshape(1, -1), N_xx_diag[i].reshape(1, 1), [labels_roi[i],], [labels_33[idx] for idx in index_rec], path_output_roi_eff, Definitions.labels_path_33, figs3_A, max_cb_N, min_cb_N, 'coolwarm', border_print = True, labels_fig_names = [figs_names_N[i],], meshdirname = Definitions.meshdirname)
        #   p
        plotting.plot_efferent (p_xy_dir[i,:].reshape(1, -1), p_xx_dir_diag[i].reshape(1, 1), [labels_roi[i],], [labels_33[idx] for idx in index_rec], path_output_roi_eff, Definitions.labels_path_33, figs3_B, max_cb_p, 0, 'plasma', labels_fig_names = [figs_p_names[i],], meshdirname = Definitions.meshdirname)
        #   Delays
        plotting.plot_efferent (D_mean_xy_dir[i,:].reshape(1, -1), D_mean_dir_diag[i].reshape(1, 1), [labels_roi[i],], [labels_33[idx] for idx in index_rec], path_output_roi_eff, Definitions.labels_path_33, figs3_B, max_cb_dir, min_cb_dir, 'viridis', labels_fig_names = [figs_mean_delays_names[i]], meshdirname = Definitions.meshdirname)
        plotting.plot_efferent(D_median_xy_dir[i,:].reshape(1, -1), D_median_dir_diag[i].reshape(1, 1), [labels_roi[i],], [labels_33[idx] for idx in index_rec], path_output_roi_eff, Definitions.labels_path_33, figs3_B, max_cb_dir, min_cb_dir, 'viridis',labels_fig_names= [figs_median_delays_names[i]], meshdirname=Definitions.meshdirname)

    # - INDIRECT MATRICES
    output_paths_roi_eff_ind = os.path.join(path_output_roi_eff, 'Indirect_connectivity_' + ind_tw)
    os.chdir(output_paths_roi_eff_ind)
    name_xy = 'roi_to_33'
    name_xx = 'roi_to_roi'
    p_xy_ind = np.loadtxt('p_' + name_xy + '.txt')
    D_mean_xy_ind = np.loadtxt('peak_delay_mean_' + name_xy + '.txt')
    D_median_xy_ind = np.loadtxt('peak_delay_med_' + name_xy + '.txt')

    # local (merged-merged)
    p_xx_ind = np.loadtxt('p_' + name_xx + '.txt')
    D_mean_xx_ind = np.loadtxt('peak_delay_mean_' + name_xx + '.txt')
    D_median_xx_ind = np.loadtxt('peak_delay_med_' + name_xx + '.txt')

    # Diagonal to have just the stim one recording
    p_xx_ind_diag = np.diag(p_xx_ind).reshape(-1, 1)
    D_mean_ind_diag = np.diag(D_mean_xx_ind).reshape(-1, 1)
    D_median_ind_diag = np.diag(D_median_xx_ind).reshape(-1, 1)

    index_rec = np.loadtxt('index_y_roi_to_33.txt', dtype=int)
    labels_roi = [Ldlpfc_merge, Rdlpfc_merge, Lifg_merge, Rifg_merge]  #this should be automatized, I know it is in this order but to avoid error it shouldnt depend on that

    figs_p_names = ['LH DLPFC Probability of Indirect Connectivity', 'RH DLPFC Probability of Indirect Connectivity',
                    'LH IFG Probability of Indirect Connectivity', 'RH IFG Probability of Indirect Connectivity', ]
    figs_mean_delays_names = ['LH DLPFC Mean Peak Delays Indirect', 'RH DLPFC Mean Peak Delays Indirect', 'LH IFG Mean Peak Delays Indirect','RH IFG Mean Peak Delays Indirect']
    figs_median_delays_names = ['LH DLPFC Median Peak Delays Indirect', 'RH DLPFC Median Peak Delays Indirect', 'LH IFG Median Peak Delays Indirect', 'RH IFG Median Peak Delays Indirect']

    for i in range(len(figs_p_names)):
        #   p
        plotting.plot_efferent(p_xy_ind[i, :].reshape(1, -1), p_xx_ind_diag[i].reshape(1, 1), [labels_roi[i], ],[labels_33[idx] for idx in index_rec], path_output_roi_eff, Definitions.labels_path_33,  figs3_B, max_cb_p, 0, 'plasma', labels_fig_names=[figs_p_names[i], ], meshdirname=Definitions.meshdirname)
        #   Delays
        plotting.plot_efferent(D_mean_xy_ind[i, :].reshape(1, -1), D_mean_ind_diag[i].reshape(1, 1), [labels_roi[i], ],[labels_33[idx] for idx in index_rec], path_output_roi_eff, Definitions.labels_path_33,figs3_B, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names=[figs_mean_delays_names[i]],meshdirname=Definitions.meshdirname)
        plotting.plot_efferent(D_median_xy_ind[i, :].reshape(1, -1), D_median_ind_diag[i].reshape(1, 1),[labels_roi[i], ], [labels_33[idx] for idx in index_rec], path_output_roi_eff,Definitions.labels_path_33, figs3_B, max_cb_ind, min_cb_ind, 'viridis', labels_fig_names=[figs_median_delays_names[i]], meshdirname=Definitions.meshdirname)

    #Supplementary
    path_folder_output_Supplementary = os.path.join(output_folder, 'Supplementary')
        # Supplementary FIG4 : Probability of connectivity vs mean peak delay
        # Analysis of DLPFC and IFG direct and indirect connections to the rest of the brain
    path_folder_output_Supplementary_s4 = os.path.join(path_folder_output_Supplementary, 'Sup_4')
    fig_p_mean_pd_dir_names  = ["dir_lDLPFC_mean_pd.svg", "dir_rDLPFC_mean_pd.svg", "dir_lIFG_mean_pd.svg", "dir_rIFG_mean_pd.svg"]
    fig_p_mean_pd_ind_names  = ["ind_lDLPFC_mean_pd.svg", "ind_rDLPFC_mean_pd.svg", "ind_lIFG_mean_pd.svg", "ind_rIFG_mean_pd.svg"]
    fig_p_median_pd_dir_names = ["dir_lDLPFC_median_pd.svg", "dir_rDLPFC_median_pd.svg", "dir_lIFG_median_pd.svg", "dir_rIFG_median_pd.svg"]
    fig_p_median_pd_ind_names = ["ind_lDLPFC_median_pd.svg", "ind_rDLPFC_median_pd.svg", "ind_lIFG_median_pd.svg", "ind_rIFG_median_pd.svg"]
    label = ['lDLPFC', 'rDLPFC', 'lIFG', 'rIFG']
    colors_dir = ['r', 'r', 'g', 'g']
    colors_ind = ['b', 'b', 'm', 'm']
    for i in range(len(p_xy_dir)):
         print(label[i])
         analyze_connectivity(p_xy_dir[i, :].reshape(1, -1), D_mean_xy_dir[i, :].reshape(1, -1), path_folder_output_Supplementary_s4,
                         'Direct connectivity, p vs mean peak delay', fig_p_mean_pd_dir_names[i], colors_dir[i], (15, 45), (0, 0.3), label[i])
         analyze_connectivity(p_xy_ind[i, :].reshape(1, -1), D_mean_xy_ind[i, :].reshape(1, -1), path_folder_output_Supplementary_s4,
                         'Indirect connectivity, p vs mean peak delay', fig_p_mean_pd_ind_names[i], colors_ind[i], (120, 250), (0, 0.3),label[i])

         analyze_connectivity(p_xy_dir[i, :].reshape(1, -1), D_median_xy_dir[i, :].reshape(1, -1), path_folder_output_Supplementary_s4,
                         'Direct connectivity, p vs median peak delay', fig_p_median_pd_dir_names[i], colors_dir[i], (15, 45), (0, 0.3), label[i])
         analyze_connectivity(p_xy_ind[i, :].reshape(1, -1), D_median_xy_ind[i, :].reshape(1, -1), path_folder_output_Supplementary_s4,
                         'Indirect connectivity, p vs median peak delay', fig_p_median_pd_ind_names[i], colors_ind[i], (120, 250), (0, 0.3),label[i])

    # 4
    #%% ROI CONNECTIVITY TO FUNCTIONAL NETWORKS (FIG 4ABC)
    path_folder_output_figs_4 = os.path.join(output_folder, 'Fig4')
    os.makedirs(path_folder_output_figs_4, exist_ok=True)
    # Connectivity from roi (Lausanne x) to functional networks as defined by Yeo et al. We compute the functional networks as merged parcels of Lausanne 250.
    path_folder_funct = os.path.join(matrix_folder, 'Functional_Networks')
    #folder for object labels for plotting networks
    path_folder_labels_nets = os.path.join(path_folder_funct, 'Labels_nets')

    path_figs_4A = os.path.join(path_folder_output_figs_4, 'A')
    os.makedirs(path_figs_4A, exist_ok=True)
    path_figs_4C = os.path.join(path_folder_output_figs_4, 'C')
    os.makedirs(path_figs_4C, exist_ok=True)
    path_figs_4C_I = os.path.join(path_figs_4C, 'I')
    os.makedirs(path_figs_4C_I, exist_ok=True)
    path_figs_4C_II = os.path.join(path_figs_4C, 'II')
    os.makedirs(path_figs_4C_II, exist_ok=True)

    plotting.plot_lpfc_definition(path_figs_4A, Definitions.labels_path, Definitions.labels_33 ) #Fig4A = Fig1G

    #create label objects for networks (left and right)
    L_fn_labels = Definitions.fn_L_labels
    L_fn_index_sub_parcels = Definitions.idx_fn_L
    R_fn_labels = Definitions.fn_R_labels
    R_fn_index_sub_parcels = Definitions.idx_fn_R
    roi_color = Definitions.colors_fn + Definitions.colors_fn
    for l in range(0,len(L_fn_labels)) :
        L_sub_index = L_fn_index_sub_parcels[l]  # tuple of index net l
        L_labels_net = [labels_250[i] for i in L_sub_index]  #labels parcels lau250
        combine_labels(L_labels_net, labels_path_250, new_label_name = L_fn_labels[l],  labels_output_path = path_folder_labels_nets )
        R_sub_index = R_fn_index_sub_parcels[l]  # tuple of index net l
        R_labels_net = [labels_250[i] for i in R_sub_index]  # labels parcels lau250
        combine_labels(R_labels_net, labels_path_250, new_label_name= R_fn_labels[l],labels_output_path=path_folder_labels_nets)

    # #Plot definition of networks (FIG 4A second block)
    plotting.plot_segmentation_def(L_fn_labels+R_fn_labels, path_folder_labels_nets, Definitions.roi_dlpfc_ifg_125,labels_path_125, roi_color, path_figs_4A, 'Functional_Networks' )
    #LEFT ROI - (FIG 4C first row)
    roi_dlpfc_ifg_l = [l for l in Definitions.roi_dlpfc_ifg_125 if l.startswith('lh.')]
    os.chdir(path_folder_funct)
    p_net_eff = np.loadtxt('p_Eff_roi_L_nets.txt')
    plotting.plot_efferent_fn(p_net_eff, L_fn_labels, roi_dlpfc_ifg_l, path_folder_labels_nets, labels_path_125, path_figs_4C_I,  0.3, roi_color,  'plasma',  'cvs_avg35_inMNI152 -Lausanne250', meshdirname = Definitions.meshdirname)
    #
    # #RIGHT ROI (supplementary)
    roi_dlpfc_ifg_r = [l for l in Definitions.roi_dlpfc_ifg_125 if l.startswith('rh.')]
    os.chdir(path_folder_funct)
    p_net_eff = np.loadtxt('p_Eff_roi_R_nets.txt')
    plotting.plot_efferent_fn( p_net_eff,R_fn_labels,roi_dlpfc_ifg_r, path_folder_labels_nets, labels_path_125, path_figs_4C_I,  0.3, roi_color,  'plasma',  'cvs_avg35_inMNI152 -Lausanne250', meshdirname = Definitions.meshdirname)
    #bar plots (FIG 4C second row)
    os.chdir(path_figs_4C_II)
    plot_bars.plot_fine_data(matrix_folder)
    #bar plots rh (supplementary)
    os.chdir(path_figs_4C_II)
    plot_bars.plot_fine_data_rh(matrix_folder)
    # #%% plot bars - efferent connectivity roi segments to nets
    #bar plots (FIG4B): ant/post DLPFC + IFG and inf/sup DLPFC + IFG
    path_figs_4B = os.path.join(path_folder_output_figs_4, 'B')
    os.makedirs(path_figs_4B, exist_ok=True)
    os.chdir(path_figs_4B)
    dv = plot_bars.get_coarse_data(matrix_folder, L_s)
    print(dv.keys())
    for hemi in ('l', 'r'): #Maybe sep left and right
        plot_bars.plot_coarse_data(dv, hemi, 'p', 'CI', 0.3) #maybe automat this
        plot_bars.plot_coarse_data(dv, hemi, 'N', None, 20000)


if __name__ == "__main__":
        main()