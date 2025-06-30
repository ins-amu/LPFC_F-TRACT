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
    rec_parcellation_name = 'Lausanne2008-125'
    #Read aggregation file of base atlas
    agg_name = "agg125_125"
    output_directory_name = os.path.join(output_directory_name_base, stim_parcellation_name + '_' + rec_parcellation_name)
    os.chdir(output_directory_name)
    agg = np.load(agg_name + '.pkl', allow_pickle=True)
    zth = 5
    #LPFC Definition in Definitions.py
    #get indexes of interest
    p = Parcellation2(Definitions.parcellation_path, Definitions.res125)
    lh_LPFC = [l for l in Definitions.roi_dlpfc_ifg_125 if l.startswith('lh.')]
    index_lLPFC = [p.get_index_by_name(i) for i in lh_LPFC]
    rh_LPFC = [l for l in Definitions.roi_dlpfc_ifg_125 if l.startswith('rh.')]
    index_rLPFC = [p.get_index_by_name(i) for i in rh_LPFC]

    # All brain without LPFC
    labels_all_left = [l for l in Definitions.labels_125 if ((l.startswith('lh.') or l.startswith('Left')) and l not in Definitions.roi_dlpfc_ifg_125)]
    index_all_left = [p.get_index_by_name(i) for i in labels_all_left]
    labels_all_right = [l for l in Definitions.labels_125 if ((l.startswith('rh.') or l.startswith('Right')) and l not in Definitions.roi_dlpfc_ifg_125)]
    index_all_right = [p.get_index_by_name(i) for i in labels_all_right]

    # #############
    # # Merging    -- Merge ALL LPFC (Left // Right) and ALL brain (excluded LPFC)
    # #############
    #Left - eff
    output_directory_name_merged = os.path.join(output_directory_name_base, 'lh_AVG_LPFC_AVG_all_eff'); os.makedirs(output_directory_name_merged, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_merged, 'rec_parcels_merged.txt'),  index_all_left, fmt='%s')
    np.savetxt(os.path.join(output_directory_name_merged, 'stim_parcels_merged.txt'), index_lLPFC, fmt='%s')
    #Merging all brain in one and all LPFC in one; all valid recording parcels excluding LPFC
    agg_merged = agg.get_merged((index_lLPFC,),(index_all_left,))
    agg_merged.save(os.path.join(output_directory_name_merged, 'lh_avg_lpfc_all_eff'))
    # General atlas
    coll, prob = af.collapse(se, agg_merged, 0, 100, z_th= zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_merged, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")

    #Left - aff
    output_directory_name_merged = os.path.join(output_directory_name_base, 'lh_AVG_LPFC_AVG_all_aff'); os.makedirs(output_directory_name_merged, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_merged, 'rec_parcels_merged.txt'),  index_lLPFC, fmt='%s')
    np.savetxt(os.path.join(output_directory_name_merged, 'stim_parcels_merged.txt'), index_all_left, fmt='%s')
    #Merging all brain in one and all LPFC in one; all valid recording parcels excluding LPFC
    agg_merged = agg.get_merged((index_all_left,), (index_lLPFC,))
    agg_merged.save(os.path.join(output_directory_name_merged, 'lh_avg_lpfc_all_aff'))
    # General atlas
    coll, prob = af.collapse(se, agg_merged, 0, 100, z_th= zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_merged, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")

    #Right - eff
    output_directory_name_merged = os.path.join(output_directory_name_base, 'rh_AVG_LPFC_AVG_all_eff');os.makedirs(output_directory_name_merged, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_merged, 'rec_parcels_merged.txt'), index_all_right, fmt='%s')
    np.savetxt(os.path.join(output_directory_name_merged, 'stim_parcels_merged.txt'), index_rLPFC, fmt='%s')
    agg_merged = agg.get_merged((index_rLPFC,), (index_all_right,))
    agg_merged.save(os.path.join(output_directory_name_merged, 'rh_avg_lpfc_all_eff'))
    # General atlas
    coll, prob = af.collapse(se, agg_merged, 0, 100, z_th= zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_merged, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")

    # Right - aff
    output_directory_name_merged = os.path.join(output_directory_name_base, 'rh_AVG_LPFC_AVG_all_aff'); os.makedirs(output_directory_name_merged, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_merged, 'rec_parcels_merged.txt'),index_rLPFC , fmt='%s')
    np.savetxt(os.path.join(output_directory_name_merged, 'stim_parcels_merged.txt'), index_all_right, fmt='%s')
    agg_merged = agg.get_merged((index_all_right,),(index_rLPFC,))
    agg_merged.save(os.path.join(output_directory_name_merged, 'rh_avg_lpfc_all_aff'))
    # General atlas
    coll, prob = af.collapse(se, agg_merged, 0, 100, z_th= zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_merged, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")

    #############
    # Merging    -- Individual 125 LPFC parcels (Left // Right) and ALL brain (excluded LPFC)
    #############
    #%LEFT - eff
    output_directory_name_merged = os.path.join(output_directory_name_base, 'lh_AVG_eff');os.makedirs(output_directory_name_merged, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_merged, 'rec_parcels_merged.txt'),index_all_left, fmt='%s')
    np.savetxt(os.path.join(output_directory_name_merged, 'stim_parcels.txt'), index_lLPFC , fmt='%s')
    #Merging all brain in one and all LPFC in one; all valid recording parcels excluding LPFC
    agg_merged = agg.get_merged( [(i,) for i in index_lLPFC], (index_all_left,))
    agg_merged.save(os.path.join(output_directory_name_merged, 'lh_lpfc_all_eff'))
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_merged, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")
    # General atlas
    coll, prob = af.collapse(se, agg_merged, 0, 100, z_th= zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)

    # %LEFT -aff
    output_directory_name_merged = os.path.join(output_directory_name_base, 'lh_AVG_aff');os.makedirs(output_directory_name_merged, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_merged, 'rec_parcels.txt'), index_lLPFC, fmt='%s')
    np.savetxt(os.path.join(output_directory_name_merged, 'stim_parcels_merged.txt'), index_all_left, fmt='%s')
    # Merging all brain in one and all LPFC in one; all valid recording parcels excluding LPFC
    agg_merged = agg.get_merged( (index_all_left,), [(i,) for i in index_lLPFC])
    agg_merged.save(os.path.join(output_directory_name_merged, 'lh_lpfc_all_aff'))
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_merged, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")
    # General atlas
    coll, prob = af.collapse(se, agg_merged, 0, 100, z_th= zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)

    #%RIGHT - eff
    output_directory_name_merged = os.path.join(output_directory_name_base, 'rh_AVG_eff');os.makedirs(output_directory_name_merged, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_merged, 'rec_parcels_merged.txt'),index_all_right, fmt='%s')
    np.savetxt(os.path.join(output_directory_name_merged, 'stim_parcels.txt'), index_rLPFC , fmt='%s')
    #Merging all brain in one and all LPFC in one; all valid recording parcels excluding LPFC
    agg_merged = agg.get_merged( [(i,) for i in index_rLPFC], (index_all_right,))
    agg_merged.save(os.path.join(output_directory_name_merged,'rh_lpfc_all_eff'))
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_merged, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")
    # General atlas
    coll, prob = af.collapse(se, agg_merged, 0, 100, z_th= zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)

    #%RIGHT - aff
    output_directory_name_merged = os.path.join(output_directory_name_base, 'rh_AVG_aff');os.makedirs(output_directory_name_merged, exist_ok=True)
    np.savetxt(os.path.join(output_directory_name_merged, 'rec_parcels.txt'),index_rLPFC, fmt='%s')
    np.savetxt(os.path.join(output_directory_name_merged, 'stim_parcels_merged.txt'), index_all_right  , fmt='%s')
    #Merging all brain in one and all LPFC in one; all valid recording parcels excluding LPFC
    agg_merged = agg.get_merged( (index_all_right,), [(i,) for i in index_rLPFC])
    agg_merged.save(os.path.join(output_directory_name_merged,'rh_lpfc_all_aff'))
    # Write txt for info on SE used
    np.savetxt(os.path.join(output_directory_name_merged, 'SE_used.txt'), [Definitions.SE_path], fmt="%s")
    # General atlas
    coll, prob = af.collapse(se, agg_merged, 0, 100, z_th= zth)
    output_directory_name_collapser = coll.save_collapse_result(output_directory_name_merged)
    np.savetxt(os.path.join(output_directory_name_collapser, 'probability.txt'), prob)

if __name__ == "__main__":
    main()