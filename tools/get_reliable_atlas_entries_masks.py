import numpy as np

import tools.stats as ts
import tools.marray as tm


def _data_check(p, N_stim, N_impl, debug):
    p = tm.arrayze(p)
    N_stim = tm.arrayze(N_stim)
    if N_impl is not None:
        N_impl = tm.arrayze(N_impl)
    if debug:
        print('_data_check() returning \n{} \n{} \n{}'.format(p, N_stim, N_impl))
    return p, N_stim, N_impl

def get_reliable_prob_mask(p, N_stim, N_impl=None, alpha=0.05, min_n_suc=10, min_n_fail=10, min_n_impl=3, max_ci=0.2, debug=False):
    """
    It applies our heuristics and returns a mask saying which cases should be kept
    Returns:
        the mask for proportions to select
    """
    p, N_stim, N_impl = _data_check(p, N_stim, N_impl, debug)
    ci = ts.get_binom_ci(p, N_stim, alpha, ac_corr_flag=True)
    if debug:
        print('ci: \n {}'.format(ci))
    mask = ((N_stim * (1. - p) >= min_n_fail) | (N_stim * p >= min_n_suc))  & (ci <= max_ci)
    if N_impl is not None:
        mask = mask & (N_impl >= min_n_impl)
    return mask, ci

def get_reliable_feat_mask(N_above, N_impl=None, min_n_feat=20, min_n_impl=3, debug=False):
    mask = (N_above >= min_n_feat)
    if N_impl is not None:
        mask = mask & (N_impl >= min_n_impl)
    return mask 

