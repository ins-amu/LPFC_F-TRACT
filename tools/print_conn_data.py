#    Author: Maciej Jedynak, <maciej.jedynak@protonmail.com>
#
#    Developed for the F-TRACT project
#    Copyright: 2024-, Maciej Jedynak, Aix-Marseille Univeriste
#    Adapted dir ind to new config SAA

import numpy as np
import os
import matplotlib.pyplot as plt
import common_lpfc as cl
from plotting import plot_bars as pb


def get_connection(parc_name_lpfc, parc_name_33, res_read_data, direction):
    """direction can be either 'eff' or 'aff'"""
    names_125, names_33, p_125_33, N_125_33, CI_125_33, p_33_125, N_33_125, CI_33_125 = res_read_data
    if parc_name_lpfc is None:
        idx_lpfc = slice(None)
    else:
        idx_lpfc = names_125.tolist().index(parc_name_lpfc)
    idx_33 =    names_33.tolist().index(parc_name_33)
    if direction == "eff":
        idxs = (idx_lpfc, idx_33)
        return p_125_33[idxs], CI_125_33[idxs]
    elif direction == "aff":
        idxs = (idx_33, idx_lpfc)
        return p_33_125[idxs], CI_33_125[idxs]
    else:
        raise Exception(f"Direction can be eff or aff, not {direction}")

def print_for_large():
    # read names and data
    res_read_data = cl.read_data()
    names_125, names_33, p_125_33, N_125_33, CI_125_33, p_33_125, N_33_125, CI_33_125 = res_read_data
    # If you want to know parcel index of something:
    # print(names_125.shape)
    # print(names_125.tolist().index('lh.caudalmiddlefrontal_1'))
    # print(names_33.shape)
    # print(names_33.tolist().index('rh.parstriangularis')) #changed ctx_ preffix
    # print(p_125_33.shape)
    # print(p_33_125.shape)
    # print('names_125')
    # print("---")
    # now query data
    # for detailed connections
    # example: print(get_connection('lh.caudalmiddlefrontal_1', 'ctx-rh-parstriangularis', res_read_data, 'eff'))
    # Below print efferent for a number of recording parcels (in 33)
    # parcel_names = ('ctx-lh-superiorfrontal', 'ctx-lh-precuneus', 'ctx-lh-caudalanteriorcingulate', 'ctx-lh-insula', 'ctx-lh-inferiortemporal', 'ctx-lh-middletemporal', 'ctx-lh-superiortemporal', 'ctx-lh-bankssts')
    #parcel_names = ('Left-Amygdala', 'Right-Amygdala')
    parcel_names = ( 'lh.parstriangularis','rh.parstriangularis' ) #'lh.caudalanteriorcingulate','rh.caudalanteriorcingulate',)#('lh.insula', 'rh.insula', )#('lh.precentral','lh.caudalmiddlefrontal',)#('lh.bankssts', )#'lh.caudalanteriorcingulate','lh.rostralanteriorcingulate',
    for parcel_name in parcel_names:
        print(f'33 parcel:   {parcel_name}')
        for mode in ('eff', 'aff'):
            print(f'Mode:{mode}')
            values = get_connection(None, parcel_name, res_read_data, mode)
            # 0 below is probs
            argsort = np.argsort(values[0])
            idxs = slice(None, None, -1)
            for name, prob, CI in zip(names_125[argsort][idxs], values[0][argsort][idxs], values[1][argsort][idxs]):
                print(name.ljust(30), round(prob, 3), round(CI, 3))
        
def print_dir_indir():
    names_125, names_33 = cl.read_data()[0:2]
    suffix = "to_33.txt"
    hemi = "l"
    print(f"Hemisphere: {hemi}")

    row_labels = ['lh.dlpfc', 'rh.dlpfc', 'lh.ifg', 'rh.ifg']

    for dirindir_dirname in ("Direct_connectivity_0_50_ms", "Indirect_connectivity_100_400_ms"):
        print(f"Directedness: {dirindir_dirname}")
        dirname = os.path.join(cl.root_path, "Mean_ROI_Eff", dirindir_dirname)

        probability_arr = np.loadtxt(os.path.join(dirname, "p_roi_" + suffix))[:, :]
        CI_arr = np.loadtxt(os.path.join(dirname, "CI_roi_" + suffix))[:, :]
        delay_arr = np.loadtxt(os.path.join(dirname, "peak_delay_med_roi_" + suffix))[:, :]

        for row_idx, label in enumerate(row_labels):
            print(f"\nSource ROI: {label}")

            prob_row = probability_arr[row_idx]
            CI_row = CI_arr[row_idx]
            delay_row = delay_arr[row_idx]

            # sorted p descendant
            argsort_idx = np.argsort(prob_row)[::-1]

            for name, prob, CI, delay in zip(names_33[:][argsort_idx],
                                             prob_row[argsort_idx],
                                             CI_row[argsort_idx],
                                             delay_row[argsort_idx]):
                print(name.ljust(30), prob, CI, delay)
        print()


def print_networks_coarse():
    dv = pb.get_coarse_data(cl.root_path)
    print("probability of lpfc segments (ant, post, sup, inf, ifg) to networks (VN, SMN, DAN, VAN, LN, FPN, DN)")
    print('Probabilities')
    print('ant         post       sup        inf         ifg')
    print(dv['p']['l'].reshape(-1, 5))
    #plot heatmap for visual
    arr = dv['p']['l'].reshape(7, 5)
    plt.imshow(arr, cmap='viridis', aspect='auto')
    plt.colorbar(label='probability')
    plt.title('Heatmap p segments - networks')
    plt.yticks(np.arange(7), ['VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DN'])
    plt.xticks(np.arange(5), ['ant', 'post', 'sup', 'inf', 'ifg'])
    plt.ylabel('functional network')
    plt.xlabel('segment')
    for i in range(arr.shape[0]):  # filas
        for j in range(arr.shape[1]):  # columnas
            plt.text(j, i, f'{round(arr[i, j], 2):.2f}',
                     ha='center', va='center',
                     color='white' if arr[i, j] < arr.max() / 2 else 'black')
    plt.show()

    print('CIs')
    arr = dv['CI']['l'].reshape(-1, 5)
    plt.imshow(arr, cmap='plasma', aspect='auto')
    plt.colorbar(label='CI')
    plt.title('Heatmap CI segments - networks')
    plt.yticks(np.arange(7), ['VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DN'])
    plt.xticks(np.arange(5), ['ant', 'post', 'sup', 'inf', 'ifg'])
    plt.ylabel('functional network')
    plt.xlabel('segment')
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            plt.text(j, i, f'{round(arr[i, j], 2):.2f}',
                     ha='center', va='center',
                     color='white' if arr[i, j] < arr.max() / 2 else 'black')
    plt.show()
    print('ant         post       sup        inf         ifg')
    print(dv['CI']['l'].reshape(-1, 5))
    # e.g. dv['p']['l'] will have values in the same order as bars in the plot
    

def print_networks_fine():
    names_125 = cl.read_data()[0]
    p = np.loadtxt(os.path.join(cl.root_path, 'Functional_comb', 'p_Eff_roi_L_nets_0_100_ms.txt'))
    # it has shape (10, 7)

def main():
    print_for_large()
    # print_dir_indir()
    #print_networks_coarse()
    # print_networks_fine()

if __name__ == '__main__':
    main()

