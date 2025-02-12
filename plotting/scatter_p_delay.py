# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 17:16:17 2025

@author: avalos-alais.s
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import os

def analyze_connectivity(p_matrix, delay_matrix, output_path, title, filename, color, xlim, ylim, label):
    """
    Analyzes and plots de relationsip between probability and delay of connectivity data, saving scatter plots and computing correlations.
    
    Parameters:
    p_matrix: numpy array (p-values or probability matrix)
    delay_matrix: numpy array (corresponding delay values)
    output_path: str, directory to save plots
    title: str, title of the plot
    filename: str, filename to save the plot
    color: str, color of the scatter points
    xlim: tuple, x-axis limits
    ylim: tuple, y-axis limits
    label: str, legend label
    """
    os.makedirs(output_path, exist_ok=True)
    
    def filter_vectors_without_nan(a, b):
        mask = ~(np.isnan(a) | np.isnan(b))
        return a[mask].reshape(1, -1), b[mask].reshape(1, -1)
    
    def plot_scatter(x, y, title, filename, color, xlim, ylim, label):
        plt.figure(figsize=(8, 6))
        plt.scatter(x, y, alpha=0.7, color=color, s=50)
        plt.ylim(ylim)
        plt.xlim(xlim)
        plt.tick_params(axis='both', which='major', labelsize=16, width=3)
        plt.title(title)
        plt.legend([label])
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_linewidth(3)
        ax.spines['left'].set_linewidth(3)
        plt.show()
        plt.savefig(os.path.join(output_path, filename))
        plt.close()
    
    # Filter out NaN values
    p_matrix, delay_matrix = filter_vectors_without_nan(p_matrix, delay_matrix)
    
    # Plot scatter plot
    plot_scatter(delay_matrix, p_matrix, title, filename, color, xlim, ylim, label)
    
    # Perform linear regression and print results
    def compute_regression(x, y, label):
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x[0], y[0])
        print(f"{label}: R = {r_value:.4f}, p-value = {p_value:.4g}")
    
    compute_regression(delay_matrix, p_matrix, title)

def main():
    #TEST
    rootpath = r"F:\Test_Results"
    dir_tw = '0_50ms'
    ind_tw = '100_400ms'
    
    path_output_dlpfc = os.path.join(rootpath, 'Mean_DLPFC_Eff')
    output_paths_dlpfc = os.path.join(path_output_dlpfc, 'Diect_connectivity_'+ dir_tw )
    
#DLPFC
    os.chdir(output_paths_dlpfc)
    name_xy = 'lDLPFCto33'
    p_DLPFC_direct = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    pd_DLPFC_direct = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
    output_paths_dlpfc = os.path.join(path_output_dlpfc, 'Indirect_connectivity_'+ ind_tw )
    os.chdir(output_paths_dlpfc)
    p_DLPFC_indirect = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    pd_DLPFC_indirect = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
   
#IFG    
    path_output_ifg = os.path.join(rootpath, 'Mean_IFG_Eff')
    output_paths_ifg = os.path.join(path_output_ifg, 'Diect_connectivity_'+ dir_tw )
    os.chdir(output_paths_ifg)
    name_xy = 'lIFGto33'
    p_IFG_direct = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    pd_IFG_direct = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
    output_paths_ifg = os.path.join(path_output_ifg, 'Indirect_connectivity_'+ ind_tw )
    os.chdir(output_paths_ifg)
    p_IFG_indirect = np.loadtxt('p_' + name_xy + '.txt').reshape(1, -1)
    pd_IFG_indirect = np.loadtxt('peak_delay_mean_' + name_xy + '.txt').reshape(1, -1)
   


#Plots 
    output_path = os.path.join(rootpath,'p_delays_scatters')
    os.makedirs(output_path, exist_ok=True)
    analyze_connectivity(p_DLPFC_direct, pd_DLPFC_direct, output_path, 'Direct connectivity, p vs mean peak delay', "dir_DLPFC.svg", 'r', (20,40), (0,0.3), 'DLPFC')
    analyze_connectivity(p_DLPFC_indirect, pd_DLPFC_indirect, output_path, 'Indirect connectivity, p vs mean peak delay', "ind_DLPFC.svg", 'g', (135,225), (0,0.3), 'DLPFC')
    analyze_connectivity(p_IFG_direct, pd_IFG_direct, output_path, 'Direct connectivity, p vs mean peak delay', "dir_IFG.svg", 'b', (20,40), (0,0.3), 'IFG')
    analyze_connectivity(p_IFG_indirect, pd_IFG_indirect, output_path, 'Indirect connectivity, p vs mean peak delay', "ind_IFG.svg", 'm', (135,225), (0,0.3), 'IFG')


if __name__ == "__main__":
    main()
