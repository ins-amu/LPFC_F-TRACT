import numpy as np
import scipy.special
import scipy.stats as ss


def _corrcoef_partial(rxy, rxz, ryz):
    return (rxy - rxz * ryz) / (np.sqrt(1 - rxz * rxz) * np.sqrt(1 - ryz * ryz)) 

def corrcoef_partial(x, y, z):
    """
    z gets regressed out
    """
    assert x.shape == y.shape == z.shape
    mask_notnan_x = ~np.isnan(x)
    mask_notnan_y = ~np.isnan(y)
    mask_notnan_z = ~np.isnan(z)
    mask_notnan = mask_notnan_x & mask_notnan_y & mask_notnan_z
    rxy = np.corrcoef(x[mask_notnan], y[mask_notnan])[0, 1]
    rxz = np.corrcoef(x[mask_notnan], z[mask_notnan])[0, 1]
    ryz = np.corrcoef(y[mask_notnan], z[mask_notnan])[0, 1]
    return _corrcoef_partial(rxy, rxz, ryz), rxy, rxz, ryz

def corrcoef_regout(x, y, z):
    # I confirmed it gives the same as partial
    rx = ss.linregress(z, x)
    slope_x, intercept_x = rx[0], rx[1]
    xc = x - z * slope_x - intercept_x
    ry = ss.linregress(z, y)
    slope_y, intercept_y = ry[0], ry[1]
    yc = y - z * slope_y - intercept_y
    xt = slope_x * z + intercept_x
    return np.corrcoef(xc, yc)

def z_binom_test(p1, p2, n1, n2, abs_flag=False):
    """
    Returns z value for two proportions, works for vectors
    If abs_flag is False, fFirst (p1) goes the one we suppose
    to be greater, otherwise it will return a negative value
    see https://www.itl.nist.gov/div898/handbook/prc/section3/prc33.htm
    """
    # below commented out because they might be scalars too
    # assert p1.shape == p2.shape == n1.shape == n2.shape, str(p1.shape) + " " + str(p2.shape) + " " + str(n1.shape) + " " + str(n2.shape)
    p = ((n1 * p1) + (n2 * p2)) / (n1 + n2)
    p1_minus_p2 = p1 - p2
    if abs_flag:
        p1_minus_p2 = np.abs(p1_minus_p2)
    z = p1_minus_p2 / (np.sqrt(p * (1 - p) * (1/n1 + 1/n2)))
    return z

def pv_binom_test(p1, p2, n1, n2, two_sided_flag=True):
    """
    For two-sided it returns p-value that p1 and p2 are different
    For one-sided p1 should be greater than p2
    """
    z = z_binom_test(p1, p2, n1, n2, abs_flag=two_sided_flag)
    # -1 below, because we want the small surface to the left, not the large one to the right
    z = -1 * np.abs(z)
    if two_sided_flag:
        sided_factor = 2
    else:
        sided_factor = 1    
    r = sided_factor * ss.norm.cdf(z)
    isnan_r = np.isnan(r)
    # NOTE here I allow for pv = 1
    assert np.all((r[~isnan_r] >= 0) & (r[~isnan_r] <= 1)), "All pv :\n" + str(r)
    return r

def nan_false_discovery_control(ps):
    r = np.full(ps.shape, np.nan)
    isnan_mask = np.isnan(ps)
    p_corr = ss.false_discovery_control(ps[~isnan_mask])
    r[~isnan_mask] = p_corr
    return r

def get_binom_ci(p, N, alpha, ac_corr_flag=False, debug=False):
    """
    Returns the Wald or Agresti–Coull confidence interval
    Arguments: 
        p - prob for success
        N - total number of tries
        alpha - confidence level (alpha/2 will be used for z lookup)
        ac_corr_flag - should the Agresti–Coull correction be applied to the Wald interval
    Returns:
        CI, so that it will be p +/- CI
    """
    z_alpha = ss.norm.ppf(1 - alpha/2)
    if ac_corr_flag:
        ns = p * N
        N = N + z_alpha * z_alpha
        p = (1. / N) * (ns + z_alpha * z_alpha / 2.)       
    return z_alpha * np.sqrt(p * (1. - p) / N)

def test_get_binom_ci(corr):
    p = 0.1
    N = 11
    res = get_binom_ci(p, N, 0.05, corr)
    print(res)
    # return
    M = 31
    p = np.linspace(0., 1., M)
    N = np.linspace(1, 100, 2 * M).astype(int)
    pm, Nm = np.meshgrid(p, N)
    print(p)
    print(N)
    np.set_printoptions(linewidth=np.nan)
    np.set_printoptions(threshold=np.inf)
    res = get_binom_ci(pm, Nm, 0.05, corr)
    # print(res)
    # print(res/pm)
    print(((Nm * (1. - pm) > 10) | (Nm * pm > 10)) & (res < 0.15))
    # print(res/pm < 0.5)
    # print(res < 0.1)
    # print((Nm * (1. - pm) > 10))


def _check_np(p, N):
    if isinstance(p, np.ndarray):
        pass
    else:
        p = np.array(p)
    if isinstance(N, np.ndarray):
        pass
    else:
        N = np.array(N)
    nanmask = np.isnan(p)
    assert np.all(nanmask ^ ((p >= 0) & (p <= 1.)))
    assert p.shape == N.shape

def get_mask_ci_max(p, N, alpha, CI_max=None, ratio_max=None, debug=False):
    """
    Computes CI and then mask of entries to be kept
    ratio_max is CI/p max
    """
    _check_np(p, N)
    CI = get_binom_ci(p, N, alpha, ac_corr_flag=True)
    mask = np.full(p.shape, False)
    if CI_max is not None:
        mask[CI <= CI_max] = True
    if ratio_max is not None:
        mask[CI <= p * ratio_max] = True
    return mask, CI
    
def get_mask_binom_norm_approx(p, N, limit=5, suc=True, debug=False):
    """
    If suc is False it will be for failures
    Check if binomial approximation can be used (for limit = 5)
    """
    mask = np.full(p.shape, False)
    if suc:
        multiplier = p
    else:
        multiplier = (1. - p)
    mask[(multiplier * N) > limit] = True
    return mask

def get_mask_mult_comp_ci_zero(p, N, alpha, debug=False):
    """
    Arguments:
        p - proportion (probability)
        N - the number of trials, can be an array, has to be the same shape as p
        alpha - significance level
    Returns:
        A boolean array of the same shape as k and N having False for cases that should be rejected
        The idea is here to check for CIs reaching to zero
    """
    _check_np(p, N)
    z = np.sqrt((N * p) / (1. - p)) 
    if debug:
        print('z = ', z)
    p_value = 1. - scipy.special.erf(z * np.sqrt(0.5))
    if debug:
        print('p_value = ', p_value)
    p_value_corrected = nan_false_discovery_control(p_value)
    if debug:
        print('p_value corrected = ', p_value_corrected)
    mask = np.full(p.shape, False)
    mask[p_value_corrected < alpha] = True
    return mask, p_value_corrected

def cohen_d(x, y):
    """
    x and y are two distributions
    """
    n1, n2 = np.sum(~np.isnan(x)), np.sum(~np.isnan(y))
    s1, s2 = np.nanvar(x, ddof=1), np.nanvar(y, ddof=1)
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    u1, u2 = np.nanmean(x), np.nanmean(y)
    return (u1 - u2) / s

def cohen_h(p1, p2):
    """
    It will return a signed (directed) value
    """
    def phi(p):
        return 2 * np.arcsin(np.sqrt(p)) 
    return phi(p2) - phi(p1)   

if __name__ == "__main__":
    # test_get_binom_ci(False)
    # test_get_binom_ci(True)
    print(pv_binom_test(0.3, 0.0, 196, 134, two_sided_flag=False))
    
