import numpy as np
import os.path

def arrayze(v):
    if isinstance(v, np.ndarray):
        return v
    else:
        return np.array([v])

def isin(a, v):
    """
    Works like numpy.isin() but a has to be sorted.
    Just that here on 1-dim
    Returns a mask of the same len as v
    CAREFUL numpy.isin() would return the mask of len a
    """
    assert a.ndim == 1
    assert v.ndim == 1
    indeces = np.searchsorted(a, v)
    return a[indeces % len(a)] == v

def _remove_nans(a1, a2):
    """
    If a1 has a NaN on position i, we remove this element in a1 and a2
    and update element a2[i+1] so that a2[i+1] := a2[i] + a2[i+1],
    in other words, the element to the right increases by the removed one
    Returns updated a1, a2
    """
    assert a1.ndim == 1
    assert a1.shape == a2.shape
    # now we need to loop from the end and remove nans
    while True:
        if np.isnan(a1[-1]):
            # The below should create new a1, a2, not manipulate the originals
            a1 = a1[:-1]
            a2 = a2[:-1]
        else:
            break
    idx_nan = np.isnan(a1)
    idx_notnan = ~idx_nan
    idx_nan_w = np.where(idx_nan)[0]
    r2 = a2.copy()
    r2[idx_nan_w + 1] += r2[idx_nan_w]    
    r1 = a1[idx_notnan]
    r2 = r2[idx_notnan]
    return r1, r2
    
def remove_nans(a1, a2):
    assert a1.ndim == 1
    assert a1.shape == a2.shape
    # Below it is nancumsum not to return any nans
    a2_cumsum = np.nancumsum(a2)
    idx = ~np.isnan(a1)
    return a1[idx], a2_cumsum[idx]

def get_extreme_indeces(a, k):
    """
    a can contain nans
    k is the number of entries to return, e.g. 2 is top 2 (largest), -2 is bottom 2 (smallest)
    Will returs indeces to be used on the original array a[*argwhere.T]
    Always the first entry will be the extreme
    """
    # It will first work on flattened array without nans, to find top/bottom values
    # and them find entries larger / smaller then kth in the original array
    af = a.flatten()
    afn = af[~np.isnan(af)]
    argsort = afn.argsort()
    if k > 0:
        idx = argsort[-k:][::-1]
        argwhere = np.argwhere(a >= afn[idx[k-1]])
        # by now argwhere would be in order in which elements ocurr in the original matrix
        # but we want them sorted
        idx2 = a[*argwhere.T].argsort()[::-1]
    elif k < 0:
        idx = argsort[:-k]
        argwhere = np.argwhere(a <= afn[idx[-k-1]])
        # see above
        # idx2 = a[argwhere[:,0], argwhere[:,1]].argsort()
        idx2 = a[*argwhere.T].argsort()
    argwhere_sorted = argwhere[idx2]
    return argwhere_sorted
    
    

def group(ax, ay, l, over, f, debug=False):
    """
    Arguments:
        ax, ay are a 1-dim arrays of same length
        l is the length of the bin
        over is the length of the lapping between bins
        f grouping function
    Returns:
        array of len n of aggregated results of f
        array of len n of number of values aggregated in each bin
        array of centres of bins
    """
    assert ax.ndim == 1
    assert ax.shape == ay.shape
    ax_max = np.nanmax(ax)
    if debug:
        print('group() debug: ax_max', ax_max)
    ax_min = 0#np.min(ax)
    start = ax_min
    end = ax_min + l
    shift = l - over
    r1, r2, r3 = [], [], []
    while start < ax_max:
        idx = (ax >= start) & (ax < end)
        vs = ay[idx]
        r1.append(f(vs))
        r2.append(len(vs))
        r3.append((start + end) / 2.)
        if debug:
            print('group() debug:', r1, r2, r3)
        start += shift
        end += shift
    return r1, np.array(r2), np.array(r3)


def merge_in_array(a, N, rows_merge_list, cols_merge_list, N_imp=None, propagate_nan_flag=False, debug=False):
    """
    Merges rows / columns in an array. Merged cells are scaled by their respective N and added. 
    This method will only return merged rows / columns, not the whole array.
    It is to avoid potential confusion with later indexing.
    Arguments:
        a is an array which can be rectangular
        N is the merge scaling factor, should have the same shape as a, does not have to be int
        rows_merge_list should be of the form [(a_idx_i, a_idx_j), (a_idx_k, a_idx_l),]
            where a_idx_i means parcel i index. The above intruction means that stim parcels
            i and j should be merged, as well as k and l. 
            If None, no merging is performed.
        cols_merge_list - as above, but for columns
        N_imp is the number of implanted patiens, should be merged as N. 
        propagate_nan_flag = False is needed to merge cells where one is NaN. 
            If it is True, NaN will propagate in merging
    Returns:
        pm - 'probability merged' of shape (len(rows_merge_list), len(cols_merge_list)),
            unless any of the lists were None, then the original shape of p is used.
        Nm - 'N merged' (number of stimulations per merged parcels) of the same shape as pm
        Nim - N_imp merged 
    """
    def flatten(l):
        return [e2 for e in l for e2 in e]
    def sum_rows(a, rows_merge_list):
        if rows_merge_list is None:
            return a
        # rs is rows summed
        ars = np.full((len(rows_merge_list), a.shape[1]), np.nan)
        last_row_idx = 0
        for i, rows_to_merge in enumerate(rows_merge_list):
            row = np.sum(a[last_row_idx:last_row_idx + len(rows_to_merge)], axis=0)
            if debug:
                print('merge_in_array(): rows_to_merge:\n', rows_to_merge)
                print('row', row)
            ars[i] = row
            last_row_idx += len(rows_to_merge)
        if debug:
            print('merge_in_array().sum_rows(): ars:\n', ars)
        return ars
    def sum_cols(a, cols_merge_list):
        if cols_merge_list is None:
            return a
        if rows_merge_list is None:
            x_len = len(a)
        else:
            x_len = len(rows_merge_list)
        ar = np.full((x_len, len(cols_merge_list)), np.nan)
        last_col_idx = 0
        for i, cols_to_merge in enumerate(cols_merge_list):
            col = np.sum(a[:, last_col_idx:last_col_idx + len(cols_to_merge)], axis=1)
            if debug:
                print('merge_in_array(): cols_to_merge:\n', cols_to_merge)
                print('col', col)
            ar[:, i] = col
            last_col_idx += len(cols_to_merge)
        if debug:
            print('merge_in_array().sum_cols(): ar:\n', ar)
        return ar
    def validate_merge_list(l):
        """validates the passed list"""
        if l is None:
            return
        for e in l:
            assert isinstance(e, tuple) or isinstance(e, list), type(e)
            for i in range(len(e)):
                isinstance(e[i], int)
    # entry assertions
    assert a.shape == N.shape
    if N_imp is not None : assert N.shape == N_imp.shape
    assert np.sum(np.isnan(N)) == 0
    if N_imp is not None : assert np.sum(np.isnan(N_imp)) == 0 
    assert np.all(N >= 0), 'some N smaller than zero or NaN found'
    # If the below fails, maybe you should have passed p*N instead of N
    assert np.all(np.isnan(a) == (N == 0))
    if N_imp is not None:
        assert np.all(N_imp >= 0)
    validate_merge_list(rows_merge_list)
    validate_merge_list(cols_merge_list)
    # computation
    a_upscaled = a * N
    # this is the only place when we condition on this flag
    # before returning the result it is not necessary to check which zeros should be 
    # changed to NaNs, because we require here that N is 0, whenever there is a NaN in data,
    # thus, the result will be 0/0 = NaN
    if not propagate_nan_flag:
        a_upscaled[np.isnan(a)] = 0
    if debug:
        print('merge_in_array(): a:\n', a)
        print('merge_in_array(): N:\n', N)
        if N_imp is not None : print('merge_in_array(): N_imp:\n', N_imp)
        print('merge_in_array(): a_upscaled:\n', a_upscaled)
    # take slices of the a and N matrices 
    if rows_merge_list is not None:
        rows_idxs = flatten(rows_merge_list)
    else:
        rows_idxs = np.arange(a.shape[0])
    if cols_merge_list is not None:
        cols_idxs = flatten(cols_merge_list)
    else:
        cols_idxs = np.arange(a.shape[1])
    idxs = np.ix_(rows_idxs, cols_idxs)
    if debug:
        print('merge_in_array(): idxs:\n', idxs)
    # below _s is 'sliced' (a subarray)
    a_upscaled_s = a_upscaled[idxs]
    N_s = N[idxs]
    if N_imp is not None : N_imp_s = N_imp[idxs]
    if debug:
        print('merge_in_array(): a_upscaled_s:\n', a_upscaled_s)
        print('merge_in_array(): N_s:\n', N_s)
        if N_imp is not None : print('merge_in_array(): N_imp_s:\n', N_imp_s)
    # submatrices are taken, now we need to merge rows / cols
    # sum rows
    a_upscaled_s_rm = sum_rows(a_upscaled_s, rows_merge_list)
    N_s_rm = sum_rows(N_s, rows_merge_list)
    if N_imp is not None : N_imp_s_rm = sum_rows(N_imp_s, rows_merge_list)
    # sum columns
    a_upscaled_s_rm_cs = sum_cols(a_upscaled_s_rm, cols_merge_list)
    N_s_rm_cs = sum_cols(N_s_rm, cols_merge_list)
    if N_imp is not None : N_imp_s_rm_cs = sum_cols(N_imp_s_rm, cols_merge_list)
    # downscale
    a_merged_dowscaled = a_upscaled_s_rm_cs / N_s_rm_cs
    if debug:
        print('merge_in_array(): a_merged_dowscaled:\n', a_merged_dowscaled)
    if rows_merge_list:
        assert a_merged_dowscaled.shape[0] == len(rows_merge_list)
        assert N_s_rm_cs.shape[0] == len(rows_merge_list), str(N.shape)
        if N_imp is not None : assert N_imp_s_rm_cs.shape[0] == len(rows_merge_list), str(N_imp.shape)
    if cols_merge_list:
        assert a_merged_dowscaled.shape[1] == len(cols_merge_list), str(a_merged_dowscaled.shape[1]) + " " + str(len(cols_merge_list))
        assert N_s_rm_cs.shape[1] == len(cols_merge_list)
        if N_imp is not None : assert N_imp_s_rm_cs.shape[1] == len(cols_merge_list)
    if N_imp is None : return a_merged_dowscaled, N_s_rm_cs
    else : return a_merged_dowscaled, N_s_rm_cs, N_imp_s_rm_cs

def copy_q2q(a, indeces=None, q1q4_flag=True, q2q3_flag=True):
    """
    Works on a square numpy array
    Might copy q1 to q4 or q2 to q3 
    """
    assert a.shape[0] == a.shape[1]
    if indeces is None:
        start_index = 0
        end_index = a.shape[0]
    else:
        start_index, end_index = indeces[0], indeces[1]
    ei = end_index - start_index
    assert ei % 2 == 0
    ei2 = int(ei / 2)
    if q1q4_flag:
        a[start_index+ei2:end_index, start_index+ei2:end_index] = a[start_index:start_index+ei2, start_index:start_index+ei2]
    if q2q3_flag:
        a[start_index+ei2:end_index, start_index:start_index+ei2] = a[start_index:start_index+ei2, start_index+ei2:end_index]
    
def overwrite_upper_with_lower(a, indeces=None, N=None, multiply=1., q1q4_flag=True, q2q3_flag=True, propagate_nan_flag=False):
    """
    Works on a square numpy array
    If q1q4_flag it writes into q1 the weighted mean of q1 and q4 
    If q2q3_flag it writes into q2 the weighted mean of q2 and q3 
    N and multiply are for weighting
    If there is a NaN in a, it will be propagated, unless propagate_nan_flag is set to False
    In such case, only two NaNs should give a NaN. NaN + X => X
    """
    assert a.shape[0] == a.shape[1]
    if N is None:
        N = np.ones(a.shape)
    else:
        assert a.shape == N.shape
    if indeces is None:
        start_index, end_index = 0, a.shape[0]
    else:
        start_index, end_index = indeces[0], indeces[1]
    ei = end_index - start_index
    assert ei % 2 == 0
    ei2 = int(ei / 2)
    q1_idx1 = start_index
    q1_idx2 = start_index + ei2
    q4_idx1 = start_index + ei2
    q4_idx2 = end_index
    def subst(a, s1, s2):
        if not propagate_nan_flag:
            nan_mask_1 = np.isnan(a[s1])
            nan_mask_2 = np.isnan(a[s2])
            a[s1][nan_mask_1] = 0
            a[s2][nan_mask_2] = 0
            nan_mask = nan_mask_1 & nan_mask_2
        a[*s1] = multiply * (  N[*s1] * a[*s1] + N[*s2] * a[*s2]) / (N[*s1] + N[*s2])
        if not propagate_nan_flag and (np.any(nan_mask) or np.any(nan_mask_2)):
            a[*s1][nan_mask] = np.nan
            a[*s2][nan_mask_2] = np.nan
    if q1q4_flag:
        s1 = (slice(q1_idx1, q1_idx2), slice(q1_idx1, q1_idx2))
        s2 = (slice(q4_idx1, q4_idx2), slice(q4_idx1, q4_idx2))
        subst(a, s1, s2)
    if q2q3_flag:
        s1 = (slice(q1_idx1, q1_idx2), slice(q4_idx1, q4_idx2))
        s2 = (slice(q4_idx1, q4_idx2), slice(q1_idx1, q1_idx2))
        subst(a, s1, s2)

