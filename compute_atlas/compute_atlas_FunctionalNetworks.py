"""
LPFC Effective connectivity to Functional Networks - Computed for LPFC roi in Lausanne2008-125 and Functional networks recording parcels as defined in Definitions (Lausanne2008-250)
Merging base aggregation file created with 'compute_atlas_LPFC_Eff_Aff'
"""

import os.path
import numpy as np
import Definitions
from stimulation_essentials.stimulation_essentials import StimulationEssentials
import compute_atlas_functions as af

def main():
    se = StimulationEssentials(Definitions.SE_path, use_indeces=False)
    se.zero_negative_times()
    output_directory_name_base= r'F:\FTRACT\Data_LPFC_FTRACT'
    #########################
    # Efferent connectivity - For functional network analysis 125-250 (merged)
    #########################
    stim_parcellation_name = 'Lausanne2008-125'
    rec_parcellation_name = 'Lausanne2008-250'
    agg_name = "agg125_250"
    output_directory_name = os.path.join(output_directory_name_base, stim_parcellation_name + '_' + rec_parcellation_name)
    os.chdir(output_directory_name)
    agg = np.load(agg_name + '.pkl', allow_pickle=True)
    zth = 5
    #########################
    # Merge functional networks
    #########################
    output_directory_name_fn = os.path.join(output_directory_name_base, 'Functional_Networks');os.makedirs(output_directory_name_fn, exist_ok=True)
    agg_fn = agg.get_merged([(i,) for i in Definitions.index_dlpfc_ifg], Definitions.idx_fn_L + Definitions.idx_fn_R)
    agg_fn.save(os.path.join(output_directory_name_fn, agg_name + '_fn'))
    np.savetxt(os.path.join(output_directory_name_fn, 'stim_parcels.txt'), Definitions.roi_dlpfc_ifg_125, fmt='%s')
    all_labels = Definitions.fn_L_labels + Definitions.fn_R_labels
    all_labels_array = np.array(all_labels).reshape(-1, 1)  #column
    np.savetxt(os.path.join(output_directory_name_fn, 'rec_parcels_merged.txt'), all_labels_array, fmt='%s')
    coll, prob= af.collapse(se, agg_fn, 0, 100, z_th=zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_fn)
    #mat_folder_directory_name = os.path.join(output_directory_name_fn, 'atlas_mats'); os.makedirs(mat_folder_directory_name, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_collapser, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")

    #################
    # Merged segments
    #################
    output_directory_name_fn_segments = os.path.join(output_directory_name_base, 'Functional_Networks', 'Segmented_roi'); os.makedirs(output_directory_name_fn_segments, exist_ok=True)
    # segments_to_merge :
    rows_merge = ['lh.dlpfc_ant', 'lh.dlpfc_post', 'lh.dlpfc_sup', 'lh.dlpfc_inf', 'lh.ifg', 'rh.dlpfc_ant', 'rh.dlpfc_post', 'rh.dlpfc_sup', 'rh.dlpfc_inf',  'rh.ifg']
    # in merge : rows left segments, right segments. in columns: left fn, right fn
    agg_fn_seg = agg.get_merged((tuple(Definitions.index_ant_dlpfc_L), tuple(Definitions.index_post_dlpfc_L),
                                 tuple(Definitions.index_sup_dlpfc_L), tuple(Definitions.index_inf_dlpfc_L),
                                 tuple(Definitions.index_ifg_L),

                                 tuple(Definitions.index_ant_dlpfc_R), tuple(Definitions.index_post_dlpfc_R),
                                 tuple(Definitions.index_sup_dlpfc_R), tuple(Definitions.index_inf_dlpfc_R),
                                 tuple(Definitions.index_ifg_R)
                                 ),
                                 Definitions.idx_fn_L + Definitions.idx_fn_R)
    agg_fn_seg.save(os.path.join(output_directory_name_fn_segments, agg_name + '_fn_segments'))


    np.savetxt(os.path.join(output_directory_name_fn_segments, 'stim_parcels_merged.txt'), rows_merge, fmt='%s')
    all_labels = Definitions.fn_L_labels + Definitions.fn_R_labels
    all_labels_array = np.array(all_labels).reshape(-1, 1)  # column
    np.savetxt(os.path.join(output_directory_name_fn_segments, 'rec_parcels_merged.txt'), all_labels_array, fmt='%s')
    coll, prob = af.collapse(se, agg_fn_seg, 0, 100, z_th= zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_fn_segments)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_collapser, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")


if __name__ == "__main__":
    main()
