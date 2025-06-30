"""
Atlas computation of Efferent and Afferent LPFC normal connectivity.
Lausanne parcellations : 125 for stimulated and 33 for recording parcels.
                         125 for stimulated and 125 for recording parcels.
                         125 for stimulated and 250 for recording parcels.
                         33  for stimulated and 33 for recording parcels.
0-100 ms
"""

import os.path
import numpy as np
import json
import Definitions
import compute_atlas_functions as af
from stimulation_essentials.stimulation_essentials import StimulationEssentials
from stimulation_essentials import filter_se

def main():
    se = StimulationEssentials(Definitions.SE_path, use_indeces=False)

    se.zero_negative_times()
    ############
    # Filter
    ############
    filter_kwargs = {
                'age_range':(0, np.inf),
                'spike_count_max':np.inf,
                'spike_count_min':-1,
                'spike_count_frac':0,
                'include_with_no_spike_count':True,
                'segmentations_stim': [0, 1, 2], #1 grey 2 white
                'segmentations_rec': [0, 1, 2],
                'dist_min':0,
                'only_symmetric_flag':False,
                'crf_exclusive_list':None,
                'allowed_centres_list':None,
                'fibre_condition':None,
                'intensity_range':(0, np.inf),
                'frequency_range':(0, np.inf),
                'pulse_width_range':(0, np.inf),
                'only_different_electrodes_flag':False,
                'min_x_coordinate':2,
                'allowed_contact_ids':None,

    }
    filter_idxs = filter_se.filter_se(se, **filter_kwargs)
    output_directory_name_base = r'F:\FTRACT\Data_LPFC_FTRACT'
    stim_parcellation_name = 'Lausanne2008-125'
    rec_parcellation_name = 'Lausanne2008-33'
    agg_name = "agg125_33"
    zth = 5
    output_directory_name = os.path.join(output_directory_name_base, stim_parcellation_name + '_' + rec_parcellation_name);os.makedirs(output_directory_name, exist_ok=True)
    with open(os.path.join(output_directory_name, 'filter_conf.json'), 'w', encoding="utf-8", newline='\r\n') as fp: json.dump(filter_kwargs, fp, indent=4,)
    agg = af.aggregate(se,
                    stim_parcellation_name,
                    rec_parcellation_name,
                    output_directory_name,
                    agg_name,
                    filter_idxs,
                    True)
    print('The aggregate:')
    print(agg)
    coll, prob= af.collapse(se, agg, 0, 100)
    output_directory_name_collapser= coll.save_collapse_result(output_directory_name)
    #mat_folder_directory_name = os.path.join(output_directory_name, 'atlas_mats');os.makedirs(mat_folder_directory_name, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    #Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_collapser, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")

    ####################
    #Local connectivity#
    ####################
    agg_name = "agg125_125"
    output_directory_name = os.path.join(output_directory_name_base, stim_parcellation_name + '_' + stim_parcellation_name); os.makedirs(output_directory_name, exist_ok=True)
    with open(os.path.join(output_directory_name, 'filter_conf.json'), 'w', encoding="utf-8", newline='\r\n') as fp:json.dump(filter_kwargs, fp, indent=4, )
    agg = af.aggregate(se,
                    stim_parcellation_name,
                    stim_parcellation_name,
                    output_directory_name,
                    agg_name,
                    filter_idxs,
                    True)
    print('The aggregate:')
    print(agg)
    coll, prob= af.collapse(se, agg, 0, 100, z_th=zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name)
    #mat_folder_directory_name = os.path.join(output_directory_name, 'atlas_mats'); os.makedirs(mat_folder_directory_name, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_collapser, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")

    #########################
    # Afferent connectivity
    #########################
    stim_parcellation_name = 'Lausanne2008-33'
    rec_parcellation_name = 'Lausanne2008-125'
    agg_name = "agg33_125"
    output_directory_name = os.path.join(output_directory_name_base, stim_parcellation_name + '_' + rec_parcellation_name);os.makedirs(output_directory_name, exist_ok=True)
    with open(os.path.join(output_directory_name, 'filter_conf.json'), 'w', encoding="utf-8", newline='\r\n') as fp: json.dump(filter_kwargs, fp, indent=4, )
    agg = af.aggregate(se,
                    stim_parcellation_name,
                    rec_parcellation_name,
                    output_directory_name,
                    agg_name,
                    filter_idxs,
                    True)
    print('The aggregate:')
    print(agg)
    coll, prob= af.collapse(se, agg, 0, 100, z_th = zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name)
    #mat_folder_directory_name = os.path.join(output_directory_name, 'atlas_mats'); os.makedirs(mat_folder_directory_name, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_collapser, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")
    #########################
    # Efferent connectivity - For functional network analysis 125-250 (merged)
    #########################
    stim_parcellation_name = 'Lausanne2008-125'
    rec_parcellation_name = 'Lausanne2008-250'
    agg_name = "agg125_250"
    output_directory_name = os.path.join(output_directory_name_base, stim_parcellation_name + '_' + rec_parcellation_name);os.makedirs(output_directory_name, exist_ok=True)
    with open(os.path.join(output_directory_name, 'filter_conf.json'), 'w', encoding="utf-8", newline='\r\n') as fp:json.dump(filter_kwargs, fp, indent=4, )
    agg = af.aggregate(se,
                    stim_parcellation_name,
                    rec_parcellation_name,
                    output_directory_name,
                    agg_name,
                    filter_idxs,
                    True)
    print('The aggregate:')
    print(agg)
    coll, prob= af.collapse( se, agg, 0, 100, z_th=zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name)
    #mat_folder_directory_name = os.path.join(output_directory_name, 'atlas_mats');os.makedirs(mat_folder_directory_name, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_collapser, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")

    # ####################
    # # Resolution analysis#
    # ####################
    agg_name = "agg33_33"
    parcellation_name = 'Lausanne2008-33'
    output_directory_name = os.path.join(output_directory_name_base, parcellation_name + '_' + parcellation_name);os.makedirs(output_directory_name, exist_ok=True)
    with open(os.path.join(output_directory_name, 'filter_conf.json'), 'w', encoding="utf-8", newline='\r\n') as fp: json.dump(filter_kwargs, fp, indent=4, )
    agg = af.aggregate(se,
                       parcellation_name,
                       parcellation_name,
                       output_directory_name,
                       agg_name,
                       filter_idxs,
                       True)
    print('The aggregate:')
    print(agg)
    coll, prob = af.collapse(se, agg, 0, 100, z_th=zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name)
    # mat_folder_directory_name = os.path.join(output_directory_name, 'atlas_mats'); os.makedirs(mat_folder_directory_name, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_collapser, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")

    #500
    agg_name = "agg500_500"
    parcellation_name = 'Lausanne2008-500'
    output_directory_name = os.path.join(output_directory_name_base, parcellation_name + '_' + parcellation_name);os.makedirs(output_directory_name, exist_ok=True)
    with open(os.path.join(output_directory_name, 'filter_conf.json'), 'w', encoding="utf-8", newline='\r\n') as fp: json.dump(filter_kwargs, fp, indent=4, )
    agg = af.aggregate(se,
                       parcellation_name,
                       parcellation_name,
                       output_directory_name,
                       agg_name,
                       filter_idxs,
                       True)
    print('The aggregate:')
    print(agg)
    coll, prob = af.collapse(se, agg, 0, 100, z_th=zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name)
    # mat_folder_directory_name = os.path.join(output_directory_name, 'atlas_mats'); os.makedirs(mat_folder_directory_name, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_collapser, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")


if __name__ == "__main__":
    main()
