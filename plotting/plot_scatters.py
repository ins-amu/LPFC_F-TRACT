#    Author: Maciej Jedynak, <maciej.jedynak@protonmail.com>
#
#    Developed for the F-TRACT project
#    Copyright: 2024-, Maciej Jedynak, Aix-Marseille Univeriste
#

import matplotlib.pyplot as plt
import scipy.stats as ss
from statsmodels.formula.api import ols
from statsmodels.formula.api import quantreg
import numpy as np
import os

import tools.stats as ts
import tools.marray as tm
import tools.common_lpfc as cd


# Not used, will wait for 2018 with LR asymmetry
pairs_lr = (
    ('lh.rostralmiddlefrontal_2',        'rh.rostralmiddlefrontal_1'),
    ('lh.rostralmiddlefrontal_3',        'rh.rostralmiddlefrontal_3'),
    ('lh.rostralmiddlefrontal_1',        'rh.rostralmiddlefrontal_2'),
    ('lh.caudalmiddlefrontal_3',         'rh.caudalmiddlefrontal_1'),
    ('lh.caudalmiddlefrontal_1',         'rh.caudalmiddlefrontal_3'),
    ('lh.caudalmiddlefrontal_2',         'rh.caudalmiddlefrontal_2'),
    ('lh.parstriangularis_1',            'rh.parstriangularis_2'), 
    ('lh.parsopercularis_1',             'rh.parsopercularis_2'),
    ('lh.parsopercularis_2',             'rh.parsopercularis_1'),
    ('lh.superiorfrontal_6',             'rh.superiorfrontal_8'))

def _plot_scatters(a1, a2, colors, diff, ax, star_mask=None, min_=0):
    """
    ...
    """
    assert a1.ndim == a2.ndim == 2
    assert a1.shape == a2.shape
    assert len(a1) == len(colors)
    if star_mask is not None:
        assert star_mask.shape == a1.shape
    regression_q_l, regression_ols_l, = [], []
    x = np.linspace(min_, 1, 2)
    # shaded area
    for i in range(len(a1)):
        # stats computation and lines plotting
        # -1 forces it to go through zero
        # regression_q = quantreg("data ~ x -1", data=dict(data=a2[i], x=a1[i])).fit(q=0.5)
        # regression_q_l.append(regression_q)
        # ax.plot(x, regression_q.params['x'] * x, color=cd.colors_lpfc[i], linestyle=':')
        regression_ols = ols("data ~ x -1", data=dict(data=a2[i], x=a1[i])).fit()
        regression_ols_l.append(regression_ols)
        ax.plot(x, regression_ols.params['x'] * x, color=cd.colors_lpfc[i], linestyle='--', alpha=0.5)
        # plot scatter
        if star_mask is None:
            ax.scatter(a1[i], a2[i], c=colors[i], label=cd.names_lpfc[i])
        else:
            ax.scatter(a1[i][~star_mask[i]], a2[i][~star_mask[i]], c=colors[i], label=cd.names_lpfc[i])
            ax.scatter(a1[i][star_mask[i]], a2[i][star_mask[i]], c=colors[i], marker='*', sizes=(400,))
    ax.plot(x, x, color='k', linewidth=2)
    ax.fill_between(x, x - diff, x + diff, alpha=0.3, color='gray')
    # ax.legend(ncol=2)
    ax.tick_params(axis='both', which='major', labelsize=20)
    return regression_q_l, regression_ols_l

# TODO: should not be in the plot module
def _get_asymmetric_connections(p1, p2, n1, n2, diff, alpha, debug=False):
    """
    ...
    """
    ts._check_np(p1, n1)
    ts._check_np(p2, n2)
    assert p1.shape == p2.shape
    print('_get_asymmetric_connections(): matrix1 avg {} +/- {}, matrix2 avg {} +/- {}'.format(np.nanmean(p1), np.nanstd(p1), np.nanmean(p2), np.nanstd(p2)))
    pvs = ts.pv_binom_test(p1, p2, n1, n2, two_sided_flag=True)
    pv_mask = pvs < alpha
    pvs_corr = ts.nan_false_discovery_control(pvs)
    pv_corr_mask = pvs_corr < alpha
    p1_p2 = p1 - p2
    diff_mask = np.abs(p1_p2) > diff
    final_mask = pv_corr_mask & diff_mask
    if debug:
        print('_get_asymmetric_connections(): the number of entries', p1.shape[0] * p1.shape[1])
        print('_get_asymmetric_connections(): pv mask was passed by', np.nansum(pv_mask))
        print('_get_asymmetric_connections(): pv mask corr was passed by', np.nansum(pv_corr_mask))
        print('_get_asymmetric_connections(): diff mask corr was passed by', np.nansum(diff_mask))
        print('_get_asymmetric_connections(): final mask corr was passed by', np.nansum(final_mask))
    return final_mask, p1_p2
    
def test2():
    k = 2
    a = np.array([[1, 3], [2, 4]])
    k_indices = np.argpartition(-a, k-1, axis=-1)[:, :k]
    # adjust indices to apply in flat array
    adjuster = np.arange(a.shape[0]) * a.shape[1]
    adjuster = np.broadcast_to(adjuster[:, None], k_indices.shape)
    k_indices_flat = k_indices + adjuster
    k_values = a.flatten()[k_indices_flat]
    print(k_values)

def main():
    # const
    diff = 0
    min_ = -0.01
    alpha = 0.05
    # read data
    names_125, names_33, p_125_33, N_125_33, CI_125_33, p_33_125, N_33_125, CI_33_125 = cd.read_data()
    if True:
        print(names_125.shape)
        print(names_125.tolist().index('lh.caudalmiddlefrontal_1'))
        print(names_33.shape)
        print(names_33.tolist().index('ctx-rh-parstriangularis'))
    ##############
    # eff vs aff 
    ##############
    # The convention below is efferent vs. afferent
    asym_mask, p1_p2 = _get_asymmetric_connections(p_125_33, p_33_125.T, N_125_33, N_33_125.T, diff=diff, alpha=alpha, debug=True)
    p1_p2[asym_mask] = np.nan
    print('sum of asyms', np.nansum(p1_p2))
    print('sums of rows of eff - aff difference:')
    for i, (diff_local, name) in enumerate(zip(np.nansum(p1_p2, axis=1), names_125)):
        print(name.ljust(30), round(diff_local, 3))
        # if i in (0, 3, 4, 5, 9, 14, 17):
        #     print(p1_p2[i], p_125_33[i], p_33_125.T[i])
    print('Strongest efferent difference:')
    largest_idx = tm.get_extreme_indeces(p1_p2, 100)
    for idx in largest_idx:
        if True:#not ((idx[0] < 10 and idx[1] < 33) or (idx[0] > 10 and idx[1] > 33)):
            print(p1_p2[*idx], 'is for', names_125[idx[0]], names_33[idx[1]])
    # Below afferent, because it is negative, with our convention eff - aff
    print('Strongest afferent difference:')
    smallest_idx = tm.get_extreme_indeces(p1_p2, -100)
    for idx in smallest_idx:
        if True:#not ((idx[0] < 10 and idx[1] < 33) or (idx[0] > 10 and idx[1] > 33)):
            print(p1_p2[*idx], 'is for', names_125[idx[0]], names_33[idx[1]])
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot()
    regression_q_l, regression_ols_l = _plot_scatters(p_125_33, p_33_125.T, cd.colors_lpfc, diff, ax, asym_mask, min_)
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
    
    path_output = r"F:\Results\Res_30oct\Directionality_symmetry"
    
    
    file_path = os.path.join(path_output, 'scatter_eff_aff.svg')
    os.makedirs(path_output, exist_ok=True)
    fig.savefig(file_path)
    plt.show()
    return
    ##############
    # LR indeces
    ##############
    # keys_list = list(colores_lpfc_dict.keys())
    # indeces_lr = np.zeros((len(pairs_lr), 2)).astype(int)
    # for i in range(len(pairs_lr)):
    #     for j in (0, 1):
    #         indeces_lr[i, j] = keys_list.index(pairs_lr[i][j])    
    # ##############
    # # eff, L vs R 
    # ##############
    # # NOTE I don't do lpfc vs lpfc, because it is too asymmetric.
    # # Only stim in lpfc, response elsewhere
    # half_idx = int(p_125_33.shape[1] / 2)
    # fig = plt.figure(figsize=(12, 8))
    # ax = fig.add_subplot()
    # for idx in indeces_lr:
    #     print(p_125_33[idx[0], :half_idx], p_125_33[idx[1], half_idx:])
    #     _plot_scatters(a1, a2, colors, ax) 
    # plt.show()
    
def test():
    a1 = np.array([
                    [1, 2, 3, 4],
                    [3, 4, 5, 6],
                    [5, 6, 7, 8],
                ])
    a2 = np.array([
                    [1.5, 1.5, 3.3, 4.4],
                    [2.5, 3.5, 4, 10],
                    [6  , 7  , 9, 1],
                ])
    colors = ('blue', 'red', 'green')
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot()
    regression_l, outliers_l = _plot_scatters(a1, a2, colors, ax)
    print(outliers_l)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot()
    # _plot_bars([regression.rsquared for regression in regression_l], colors, ax)
    plt.show()


if __name__ == '__main__':
    main()
    # test2()

