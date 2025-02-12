import os
import numpy as np
# import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
from matplotlib_curly_brace import curlyBrace

from collections import OrderedDict
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize

our_roi_od = OrderedDict([
                        ('lh-dlPFC',  (0,     7)),
                        ('lh-IFG',    (7,     10)),
                        ('rh-dlPFC', (10,    18)),
                        ('rh-IFG',   (18,    22)),
                        ])

lobes_od_1h = OrderedDict([
                        ('Occipital',   (0,     4)),
                        ('Temporal',    (4,     13)),
                        ('Parietal',    (13,    18)),
                        ('Cingulate',   (18,    22)),
                        ('Frontal',     (22,    33)),
                        ('Other',       (33,    41)),
                        ])

lobes_od_sofi_manual = OrderedDict([
                        ('lh_Occipital',   (0,     4)), 
                        ('lh_Temporal',    (4,     13)), 
                        ('lh_Parietal',    (13,    18)), 
                        ('lh_Cingulate',   (18,    22)), 
                        ('lh_Frontal',     (22,    33)), 
                        ('lh_Other',       (33,    36)),
                        ('rh_Occipital',   (36,    40)), 
                        ('rh_Temporal',    (40,    49)), 
                        ('rh_Parietal',    (49,    54)), 
                        ('rh_Cingulate',   (54,    59)), 
                        ('rh_Frontal',     (59,    69)), 
                        ('rh_Other',       (69,    72)),
                        ])

def get_lobes_od_2h(od):
    lobes_od_2h = OrderedDict()
    for h, i in zip(('lh', 'rh'), (0, 41)):
        for k, v in od.items():
            lobes_od_2h[h + '-' + k] = (v[0] + i, v[1] + i)
    return lobes_od_2h

def idx_reader(filename):
    with open(filename, 'r') as f:
        r = [line.strip().replace('(', '').replace(')', '').replace(',', '') for line in f]
    return np.array(r).astype(int)

def translate_od(od, ref):
    """
    ref is reference, so indeces indicating which ones (from all) we use in ours
    """
    def find_idx(v, ref_l, first_or_last):
        """
        v is a tuple
        ref_l is the reference list 
        first_or_last is 0 or -1
        """
        values_in_range = list(range(*v))
        while True:
            value = values_in_range[first_or_last]
            try:
                idx = ref_l.index(value)
            except ValueError:
                values_in_range.remove(value)
            else:
                return idx
    od_r = OrderedDict()
    ref_l = ref.tolist()
    for k, v in od.items():
        first = find_idx(v, ref_l, 0)
        second = find_idx(v, ref_l, -1) + 1
        od_r[k] = (first, second)
    return od_r


def imshow_with_braces(a, od, figsize=(16, 6), suptitle="", filename=None, **kwargs):
    """
    a is the matrix to be plotted
    od is an ordered dict od lobes indeces
    if filename is not None, the figure will be saved
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    # img = ax.imshow(a, interpolation='None', cmap='coolwarm', aspect=0.75, norm=LogNorm(), **kwargs)
    
    img = ax.imshow(a, interpolation='None', aspect=0.75, **kwargs)
    fig.colorbar(img, ax=ax, location='right', shrink=0.5)
    index_half = int(list(od.values())[-1][1] / 2)
    brace_offset = 3
    half_row = 0.5
    text_offset_var_c = -1.5
    text_offset = 4.5
    ax.set_xlim(-half_row, a.shape[1] - half_row)
    ax.set_ylim(a.shape[0] - half_row, -half_row)
    yticks_labels = list(range(0, a.shape[0], 5))
    ax.set_yticks(yticks_labels) 
    ax.set_yticklabels(yticks_labels) 
    xticks_labels = list(range(0, a.shape[1], 5))
    ax.set_xticks(xticks_labels) 
    ax.set_xticklabels(xticks_labels) 
    cb_kwargs =  {'str_text':'', 'bool_auto':False, 'color':'k', 'lw':2, 'clip_on':False}
    text_kwargs = {'verticalalignment':'center', 'horizontalalignment':'center'}
    assert list(our_roi_od.values())[-1][1] == a.shape[0]
    for (lobe, text_offset_var) in zip(od.keys(), [text_offset_var_c, 0, text_offset_var_c, 0] + [text_offset_var_c, 0, text_offset_var_c, 0, text_offset_var_c, 0] * 2):
        # below, the first dim is horizontal
        p1 = np.array([od[lobe][0] - half_row, a.shape[0] + brace_offset])
        p2 = np.array([od[lobe][1] - half_row, a.shape[0] + brace_offset]) 
        mid = (p1 + p2) / 2
        curlyBrace.curlyBrace(fig, ax, p1, p2, **cb_kwargs)
        plt.text(mid[0], mid[1] + text_offset + text_offset_var, s=lobe, **text_kwargs)
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0.15, hspace = 0, wspace = 0)
    if suptitle:
        fig.suptitle(suptitle, fontsize='xx-large')
    if filename:
        fig.savefig(filename)
    return fig, ax

def plot_aff_eff_with_braces(rootpath,output_path, var_name, cmap, norm, cb_max = None):
    """
    Hardcoded... 
    """
    ######################
    # lobes OD translation
    ######################
    ref = idx_reader(os.path.join(rootpath, "Lausanne2008-33__Lausanne2008-125", "index_y_33to125_0_100ms.txt"))
    print("ref", ref)
    lobes_od = get_lobes_od_2h(lobes_od_1h)
    print("lobes_od", lobes_od)
    lobes_od_r = translate_od(lobes_od, ref)
    increase = list(our_roi_od.values())[-1][1]
    for k in lobes_od_r.keys():
        lobes_od_r[k] = (lobes_od_r[k][0] + increase, lobes_od_r[k][1] + increase)
    print("lobes_od_r", lobes_od_r)
    print("lobes_od_sofi_manual", lobes_od_sofi_manual)
    od = OrderedDict()
    od.update(our_roi_od)
    od.update(lobes_od_r)
    print("od", od)
    
    ###################
    # data reading
    ###################
   
    a_square = np.loadtxt(os.path.join(rootpath, "Lausanne2008-125__Lausanne2008-33", var_name + "_125to125_0_100ms.txt"))
    a_eff = np.loadtxt(os.path.join(rootpath, "Lausanne2008-125__Lausanne2008-33", var_name + "_125to33_0_100ms.txt"))
    a_aff = np.loadtxt(os.path.join(rootpath, "Lausanne2008-33__Lausanne2008-125", var_name + "_33to125_0_100ms.txt"))
    tot_eff = np.hstack((a_square, a_eff))
    print('tot_eff.shape', tot_eff.shape)
    tot_aff = np.vstack((a_square, a_aff)).T
    print('tot_aff.shape', tot_aff.shape)
    #
    if norm is None and cb_max is not None:
        norm = Normalize(vmin=0, vmax=cb_max) 
    elif norm is None and cb_max is None: 
        norm = Normalize(vmin=0, vmax=1) 
    #
    os.chdir(output_path)
    imshow_with_braces(tot_eff, od, figsize=(16, 8), suptitle="Efferent connectivity", filename="eff_{}.svg".format(var_name), cmap=cmap, norm=norm)
    imshow_with_braces(tot_aff, od, figsize=(16, 8), suptitle="Afferent connectivity", filename="aff_{}.svg".format(var_name), cmap=cmap, norm=norm)
    plt.show()
    
def main():
    rootpath = r"F:\Test_Results"
    output_path = r'F:'
    plot_aff_eff_with_braces(rootpath,output_path, 'p', 'plasma', None, cb_max = 0.5)
    plot_aff_eff_with_braces(rootpath,output_path, 'N', 'coolwarm', LogNorm())
   

if __name__ == "__main__":
    main()
