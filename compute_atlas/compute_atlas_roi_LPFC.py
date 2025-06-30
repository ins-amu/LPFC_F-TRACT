"""
Atlas computation of Efferent connectivity of DLPFC and IFG merged parcels.
Lausanne 33 recording parcels
direct : 0-50ms and indirect : 100-400ms
"""

import os.path
import numpy as np
import Definitions
import compute_atlas_functions as af
from stimulation_essentials.stimulation_essentials import StimulationEssentials
from parcellation2.parcellation2 import Parcellation2


def main():
    se = StimulationEssentials(Definitions.SE_path, use_indeces=False)
    se.zero_negative_times()
    output_directory_name_base= r'F:\FTRACT\Data_LPFC_FTRACT'
    stim_parcellation_name = 'Lausanne2008-125'
    rec_parcellation_name = 'Lausanne2008-33'
    #Read aggregation file of base atlas
    agg_name = "agg125_33"
    output_directory_name = os.path.join(output_directory_name_base, stim_parcellation_name + '_' + rec_parcellation_name)
    os.chdir(output_directory_name)
    agg = np.load(agg_name + '.pkl', allow_pickle=True)
    zth = 5
    #############
    # Merging
    #############
    #LPFC Definition in Definitions.py
    #get indexes of interest
    os.chdir(Definitions.parcellation_path)
    # labels of rec parcellation (33)
    labels_all_rec = [line.strip() for line in open(rec_parcellation_name + '.txt')]
    p = Parcellation2(Definitions.parcellation_path, rec_parcellation_name)
    index_rec_parcellation = [p.get_index_by_name(i) for i in labels_all_rec]
    #Save labels (in order) in use (in this case it would be the same as normal Lausanne2008-33.txt)
    output_directory_name_merged = os.path.join(output_directory_name_base, 'DLPFC_IFG_to_' + rec_parcellation_name); os.makedirs(output_directory_name_merged, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_merged, 'rec_parcels_' + rec_parcellation_name + '.txt'),[p.get_name_by_index(i) for i in index_rec_parcellation],  fmt='%s')
    #Merging all dlpfc in one, and all ifg in one; all valid recording parcels
    agg_merged = agg.get_merged((tuple(Definitions.index_dlpfc_L), tuple(Definitions.index_dlpfc_R),
                                 tuple(Definitions.index_ifg_L), tuple(Definitions.index_ifg_R)), [(i,) for i in index_rec_parcellation])
    agg_merged.save(os.path.join(output_directory_name_merged, agg_name + '_merged'))
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_merged, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")

    #Save labels (in order) in use -  merged parcels, new parcels, give names to them.
    np.savetxt(os.path.join(output_directory_name_merged, 'stim_parcels_merged.txt'), ['lh.dlpfc','rh.dlpfc', 'lh.ifg', 'rh.ifg'], fmt='%s')
    #Direct Atlas
    coll, prob = af.collapse(se, agg_merged, 0,50, z_th= zth)
    output_directory_name_merged_direct = os.path.join(output_directory_name_merged, 'Direct'); os.makedirs(output_directory_name_merged_direct, exist_ok=True)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged_direct)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)

    #Indirect Atlas
    coll, prob = af.collapse(se, agg_merged, 100,400, z_th= zth)
    output_directory_name_merged_indirect = os.path.join(output_directory_name_merged, 'Indirect'); os.makedirs(output_directory_name_merged_indirect, exist_ok=True)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged_indirect)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)

    ############################
    #Merged local connectivity
    ############################
    agg_name = "agg125_125"
    output_directory_name = os.path.join(output_directory_name_base, stim_parcellation_name + '_' + stim_parcellation_name)
    os.chdir(output_directory_name)
    agg = np.load(agg_name + '.pkl', allow_pickle=True)
    #DLPFC-DLPFC ; IFG-IFG
    output_directory_name_merged_local = os.path.join(output_directory_name_base, 'DLPFC_IFG_to_DLPFC_IFG');os.makedirs(output_directory_name_merged_local, exist_ok=True)
    agg_merged_local = agg.get_merged((tuple(Definitions.index_dlpfc_L), tuple(Definitions.index_dlpfc_R),
                                       tuple(Definitions.index_ifg_L), tuple(Definitions.index_ifg_R)),
                                      (tuple(Definitions.index_dlpfc_L), tuple(Definitions.index_dlpfc_R),
                                       tuple(Definitions.index_ifg_L), tuple(Definitions.index_ifg_R)))
    agg_merged_local.save(os.path.join(output_directory_name_merged_local, agg_name + '_merged'))
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_merged_local, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")
    np.savetxt(os.path.join(output_directory_name_merged_local, 'stim_parcels_merged.txt'),['lh.dlpfc', 'rh.dlpfc', 'lh.ifg', 'rh.ifg'], fmt='%s')
    np.savetxt(os.path.join(output_directory_name_merged_local, 'rec_parcels_merged.txt'), ['lh.dlpfc', 'rh.dlpfc', 'lh.ifg', 'rh.ifg'], fmt='%s')

    # Direct Atlas
    coll, prob = af.collapse(se, agg_merged_local, 0, 50, z_th= zth)
    output_directory_name_merged_direct = os.path.join(output_directory_name_merged_local, 'Direct');os.makedirs(output_directory_name_merged_direct, exist_ok=True)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged_direct)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Indirect Atlas
    coll, prob = af.collapse(se, agg_merged_local, 100, 400, z_th= zth)
    output_directory_name_merged_indirect = os.path.join(output_directory_name_merged_local, 'Indirect');os.makedirs(output_directory_name_merged_indirect, exist_ok=True)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged_indirect)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)

if __name__ == "__main__":
    main()