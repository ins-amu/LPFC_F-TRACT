# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 17:03:59 2024

@author: avalos-alais.s
"""

#    Author: Maciej Jedynak, <maciej.jedynak@protonmail.com>
#
#    Developed for the F-TRACT project
#    Copyright: 2024-, Maciej Jedynak, Aix-Marseille Univeriste
#

import matplotlib.pyplot as plt
import numpy as np
import os

# import tools.stats as ts
# import tools.marray as tm
from tools import common_lpfc as cd


network_labels = ['Visual', 'Somatomotor', 'Dorsal Attention', 'Ventral Attention', 'Limbic', 'Fronto-parietal', 'Default']
network_labels = ['VN', 'SN', 'DAN', 'VAN', 'LN', 'FPN', 'DN']
# network_colors = ['#89259a', '#64aff3','#2ecd00', '#f566d6', '#ffe496', '#ff972b', '#ff1010']
network_colors = ['#89259a', '#64aff3','#2ecd00', '#f566d6', 'gold', '#ff972b', '#ff1010']

   
def _plot_grouped_bars( ax, 
                        vs, 
                        colors, 
                        hatch=None, 
                        yerr=None, 
                        xticklabels=None, 
                        xtickcolors=None,
                        yticks=None,
                        yticklabels=None, 
                        labels=None, 
                        legend_fontsize=None, 
                        legend_loc=None,
                        ticklabels_fontsize=None, 
                        xlabel=None,
                        ylabel=None,
                        label_fontsize=None,
                        ylim=None,
                        rotation=None,
                        two_axes_flag=False,
                        vertical_flag=True,
                        debug=False):
    """
    vs is a vector of values
    If colors is iterable, the len of colors is the len of a group
    If colors is a scalar, the len of hatch is the len of a group    
    first entries of vs correrspond to the first group
    """
    if isinstance(colors, list) and isinstance(hatch, list):
        assert len(colors) == len(hatch)
    if isinstance(colors, list):
        group_len = len(colors)
    else:
        group_len = len(hatch)
        colors = [colors for _ in range(len(hatch))]
    if hatch is None:
        hatch = [None] * len(colors)
    if labels is None:
        labels = [None] * len(colors)
    group_n = len(vs) / group_len
    assert (group_n % 1) == 0, group_n
    group_n = int(group_n)
    x = np.array((range(group_n)))
    bar_width = 1./(group_len + 1)
    if debug:
        print('colors', colors)
        print('group_len', group_len)
        print('group_n', group_n)
        print('height', vs[0::group_n])
        print('bar_width', bar_width)
        print('hatch', hatch)
        print('x', x)
    for i in range(group_len):
        # each iteration of this loop plots one bar from each group
        kwargs = {'x':x+i*bar_width, 'height':vs[i::group_len], 'width':bar_width, 'color':colors[i], 'hatch':hatch[i], 'edgecolor':'k', 'error_kw':{'elinewidth':2, 'capsize':5}, 'label':labels[i]}
        if yerr is not None:
            kwargs['yerr'] = yerr[i::group_len]
        # ax.barh(y=x+i*bar_width, width=vs[i::group_len], yerr=yerr[i::group_len], height=bar_width, color=colors[i], hatch=hatch[i], edgecolor='k', error_kw={'elinewidth':2, 'capsize':5}, label=labels[i])
        ax.bar(**kwargs)
    # ax.set_xticks([])
    # if just x below, label will be below the first bar
    # to move to the second bar, you need to add bar_width
    # below 0.5 to be in the middle of a group, -1 because of the gap between groups
    if vertical_flag:
        ax.set_xticks(x + 0.5 * (group_len - 1) * bar_width)
    else:
        ax.set_yticks(x + 0.5 * (group_len - 1) * bar_width)
    ax.tick_params(axis='both', which='major', labelsize=ticklabels_fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=ticklabels_fontsize)
    if xticklabels:
        if vertical_flag:
            ax.set_xticklabels(xticklabels, fontsize=ticklabels_fontsize)
        else:
            ax.set_yticklabels(xticklabels, fontsize=ticklabels_fontsize)
    else:
        if vertical_flag:
            ax.set_xticklabels([])
        else:
            ax.set_yticklabels([])
        # ax.set_xticks([])
    if xtickcolors:
        if vertical_flag:
            for color, t in zip(xtickcolors, ax.xaxis.get_ticklabels()):
                t.set_color(color) 
        else:
            for color, t in zip(xtickcolors, ax.yaxis.get_ticklabels()):
                t.set_color(color) 
    if yticklabels:
        if vertical_flag:
            ax.set_yticklabels(yticklabels, fontsize=ticklabels_fontsize)
        else:
            ax.set_xticklabels(yticklabels, fontsize=ticklabels_fontsize)
    else:
        if vertical_flag:
            ax.set_yticklabels([])
        else:
            ax.set_xticklabels([])
    if yticks:
        if vertical_flag:
            ax.set_yticks(yticks)
        else:
            ax.set_xticks(yticks)
    if ylabel:
        if vertical_flag:
            ax.set_ylabel(ylabel, fontsize=label_fontsize)
        else:
            ax.set_xlabel(ylabel, fontsize=label_fontsize)
    if xlabel:
        if vertical_flag:
            ax.set_xlabel(xlabel, fontsize=label_fontsize)
        else:
            ax.set_ylabel(xlabel, fontsize=label_fontsize)
    if labels and legend_loc:
        ax.legend(prop={'size':legend_fontsize}, loc=legend_loc)
    if ylim:
        if vertical_flag:    
            ax.set_ylim(ylim)
        else:
            ax.set_xlim(ylim)
    if rotation:
        if vertical_flag:
            plt.xticks(rotation=rotation)
        else:
            plt.yticks(rotation=rotation)
    if two_axes_flag:
        ax.spines[['right', 'top']].set_visible(False)
        ax.spines[['left', 'bottom']].set_linewidth(4)
        
    
    

def test_color():
    fig = plt.figure()
    ax = fig.add_subplot()
    vs = np.arange(0, 12)
    yerr = np.random.rand(len(vs))
    colors = ['r', 'g', 'b']
    labels = ['first', 'second', 'third', 'fourth']
    _plot_grouped_bars(ax, vs, colors, None, yerr, labels=labels)
    plt.show()

def test_hatch_one_color():
    fig = plt.figure()
    ax = fig.add_subplot()
    vs = np.arange(0, 12)
    yerr = np.random.rand(len(vs))
    colors = 'r'
    labels = ['first', 'second', 'third', 'fourth']
    hatch = ['//', '\\', '|']
    _plot_grouped_bars(ax, vs, colors, hatch, yerr, labels=labels)
    plt.show()

def test_hatch_colors():
    fig = plt.figure()
    ax = fig.add_subplot()
    vs = np.arange(0, 12)
    yerr = np.random.rand(len(vs))
    colors = ['r', 'g', 'b']
    labels = ['first', 'second', 'third', 'fourth']
    hatch = ['//', '\\', '|']
    _plot_grouped_bars(ax, vs, colors, hatch, yerr, labels=labels)
    plt.show()
    # main()

def plot_coarse_data(dv, hemi, var_name, yerr_name=None, ymax=None):
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot()
    vs = dv[var_name][hemi] 
    if yerr_name:
        yerr = dv[yerr_name][hemi]
    else:
        yerr = None
    print(vs.shape)
    # colors = ['r', 'g', 'b', 'yellow', 'magenta',]
    colors = 'lightgrey'
    hatch = ['\\\\]', '///', '..', 'xxx', '*']
    labels = ['Anterior', 'Posterior', 'Superior', 'Inferior', 'IFG']
    fontsize = 28
    _plot_grouped_bars( ax, 
                        vs, 
                        colors, 
                        hatch, 
                        yerr, 
                        xticklabels=network_labels, 
                        xtickcolors=network_colors,
                        labels=labels, 
                        legend_fontsize=fontsize, 
                        ticklabels_fontsize=fontsize, 
                        legend_loc=2, 
                        ylabel=None,#'Probabilistic connectivity',
                        label_fontsize=fontsize,
                        ylim=(0, ymax),
                        yticks=([0, ymax]),
                        yticklabels=([0, ymax]),
                        #rotation=15,
                        vertical_flag=True,
                        debug=True)
    plt.tight_layout()
    fig.savefig('coarse_data_{}_{}.svg'.format(hemi, var_name), bbox_inches='tight')
    plt.show()

###########################
#  Data reading functions #
###########################
# 
# These functions should return arrays as will be used for barplots, i.e. 
# one dimensional vectors where the first N entries correspond to the first
# group of len N.
#
 
vars_ = ('p', 'CI', 'N')
hemis = ('l', 'r')

def get_coarse_data(data_rootpath,L_s):
    """
    Reads in data needed to plot anterior/posterior, superior/inferior, IFG
    connectivity to all functional networks.
    The returned data should be a vector of p and a vector of CI, each of
    len = 5 (ant, post, inf, sup, IFG) x 14 = 70. In fact there will be 4
    such vectors: always stimulation in left, but response in L and R, thus
    p_l->l, CI_l->l, p_l->r, CI_l->r.
    
    """
    # Load arrays
    data_rootpath_local = os.path.join(data_rootpath, 'Functional_nets', 'Segmented_Roi')
    universal_s = '{}_{}_{}_nets.txt'
    sitess = (('ant', 'post'), ('sup', 'inf'))
    def _extend_d(suffix):
        for h in hemis:
            # v = _load(var, h, site, suffix)
            v = np.loadtxt(os.path.join(data_rootpath_local, subdir, universal_s.format(var, h, site + suffix)))
            d[var][h] = np.hstack((d[var][h], v))
    d = {}
    for var in vars_:
        d[var] = {}
        for h in hemis:
            d[var][h] = np.array([])
        for sites in sitess:
            subdir = '{}_{}_DLPFC_IFG'.format(*sites)
            for site in sites:
                _extend_d('DLPFC')
        site = 'IFG'
        _extend_d('')
    # reshape arrays
    # by now each vector is of len 35, begin 5 (coarse regions) x 7 (networks)
    # the order is ant->net1, ant->net2, ..., IFG->net7
    # dv is dictionary view
    dv = {}
    for var in vars_:
        dv[var] = {}
        for h in hemis:
            dv[var][h] = d[var][h].reshape(5, 7).T.reshape(35)
    # dv order is ant->net1, post->net1, ..., IFG->net7
    return dv


def plot_fine_data(data_rootpath, plot_in_one_flag=False):
    #TODO : automat for l and right 
    p = np.loadtxt(os.path.join(data_rootpath, 'Functional_nets', 'p_Eff_roi_L_nets_0_100ms.txt'))
    # take left to left
    n_parcels = 10
    n_networks = 7
    p_ll = p[:n_parcels, :n_networks].T.flatten() #the index arent necessary 
    CI = np.loadtxt(os.path.join(data_rootpath, 'Functional_nets', 'CI_Eff_roi_L_nets_0_100ms.txt'))
    CI_ll = CI[:n_parcels, :n_networks].T.flatten()
    colors = list(cd.colors_dlpfc_dict.values())[:n_parcels]
    labels = list(cd.colors_dlpfc_dict.keys())[:n_parcels]
    fontsize = 50 
    figsize = (9, 10)
    ymax = 0.25
    if plot_in_one_flag:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        _plot_grouped_bars( ax, 
                            p_ll, 
                            colors, 
                            hatch=None, 
                            yerr=CI_ll, 
                            # xticklabels=xticklabels, #ERROR
                            labels=labels, 
                            legend_fontsize=fontsize, 
                            ticklabels_fontsize=fontsize, 
                            legend_loc=2, 
                            ylabel=None,#'Probabilistic connectivity',
                            label_fontsize=fontsize,
                            ylim=(0, ymax),
                            debug=True)
        plt.tight_layout()
        fig.savefig('fine_data_ll.svg', bbox_inches='tight')
    else:
        for i in range(n_networks):
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot()
            _plot_grouped_bars( ax, 
                                p_ll[n_parcels*i:n_parcels*(i+1)], 
                                colors, 
                                hatch=None, 
                                yerr=CI_ll[n_parcels*i:n_parcels*(i+1)], 
                                xticklabels=None, 
                                yticks=[0, ymax/2, ymax],
                                yticklabels=['0', ymax/2, ymax],
                                labels=None, 
                                legend_fontsize=fontsize, 
                                ticklabels_fontsize=fontsize, 
                                legend_loc=None, 
                                xlabel=None,#xticklabels[i],
                                ylabel=None,#'Probabilistic connectivity',
                                label_fontsize=fontsize,
                                ylim=(0, ymax),
                                two_axes_flag=True,
                                debug=True)
            plt.tight_layout()
            fig.savefig('fine_data_ll_' + str(i) + '.svg', bbox_inches='tight')
    plt.show()
def plot_fine_data_rh(data_rootpath, plot_in_one_flag=False):
    #TODO : automat for l and right 
    p = np.loadtxt(os.path.join(data_rootpath, 'Functional_nets', 'p_Eff_roi_R_nets_0_100ms.txt'))
    # take left to left
    n_parcels_l = 10 
    n_parcels = 12 #10 left
    n_networks = 7
    p_rr = p.T.flatten()
    CI = np.loadtxt(os.path.join(data_rootpath, 'Functional_nets', 'CI_Eff_roi_R_nets_0_100ms.txt'))
    CI_rr = CI.T.flatten()
    colors = list(cd.colors_dlpfc_dict.values())[n_parcels_l:]
    labels = list(cd.colors_dlpfc_dict.keys())[n_parcels_l:]
    fontsize = 50 
    figsize = (9, 10)
    ymax = 0.25
    if plot_in_one_flag:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        _plot_grouped_bars( ax, 
                            p_rr, 
                            colors, 
                            hatch=None, 
                            yerr=CI_rr, 
                            # xticklabels=xticklabels, #ERROR
                            labels=labels, 
                            legend_fontsize=fontsize, 
                            ticklabels_fontsize=fontsize, 
                            legend_loc=2, 
                            ylabel=None,#'Probabilistic connectivity',
                            label_fontsize=fontsize,
                            ylim=(0, ymax),
                            debug=True)
        plt.tight_layout()
        fig.savefig('fine_data_rr.svg', bbox_inches='tight')
    else:
        for i in range(n_networks):
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot()
            _plot_grouped_bars( ax, 
                                p_rr[n_parcels*i:n_parcels*(i+1)], 
                                colors, 
                                hatch=None, 
                                yerr=CI_rr[n_parcels*i:n_parcels*(i+1)], 
                                xticklabels=None, 
                                yticks=[0, ymax/2, ymax],
                                yticklabels=['0', ymax/2, ymax],
                                labels=None, 
                                legend_fontsize=fontsize, 
                                ticklabels_fontsize=fontsize, 
                                legend_loc=None, 
                                xlabel=None,#xticklabels[i],
                                ylabel=None,#'Probabilistic connectivity',
                                label_fontsize=fontsize,
                                ylim=(0, ymax),
                                two_axes_flag=True,
                                debug=True)
            plt.tight_layout()
            fig.savefig('fine_data_rr_' + str(i) + '.svg', bbox_inches='tight')
    plt.show()

def main():
    
    output_folder = r'F:\Results\Test_Results'
    path_folder_funct = os.path.join(output_folder, 'Functional_nets')
    os.makedirs(path_folder_funct, exist_ok=True)
    path_folder_funct_fig = os.path.join(path_folder_funct, 'Figures')
    os.makedirs(path_folder_funct_fig, exist_ok = True)
    
    path_folder_funct_fig_bars = os.path.join(path_folder_funct_fig, 'Bar_Plots')
    os.makedirs(path_folder_funct_fig_bars, exist_ok = True)
    
    os.chdir(path_folder_funct_fig_bars)
    plot_fine_data_rh(output_folder)
    # coarse
    # data_rootpath = r'F:\Results\Res_30_12'
    # L_s = 'Lausanne2008'
    # dv = get_coarse_data(data_rootpath, L_s)
    # print(dv.keys())
    # os.chdir(r'F:\Results\Res_30_12\Functional_comb\Merged_Roi')
    # for hemi in ('l', 'r'):
    #     plot_coarse_data(dv, hemi, 'p', 'CI', 0.25)
    #     plot_coarse_data(dv, hemi, 'N', None, 20000)
    # # fine
    # plot_fine_data(data_rootpath)
    # test_color()
    # test_hatch_one_color()

if __name__ == '__main__':
    main()

