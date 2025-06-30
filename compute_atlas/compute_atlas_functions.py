"""
Functions for atlas computation with ft_atlas project.
Specificities for FTRACT_LPFC project
"""

import os.path
import numpy as np
import time
from stimulation_essentials.stimulation_essentials import StimulationEssentials
import stimulation_essentials.constants
from ft_atlas.aggregate.aggregate import Aggregate
from ft_atlas.collapser import Collapser
from ft_atlas import extract_functions as ef
from ft_atlas import collapse_functions as colf
from ft_atlas.condition_functions import collapse as condcf
from ft_atlas.probability import compute_probability


def aggregate(  se:StimulationEssentials,
                stim_parcellation_name:str,
                rec_parcellation_name:str,
                output_directory_name:str,
                agg_name:str,
                filter_idxs:np.ndarray,
                really_aggregate_flag:bool)->Aggregate:

    if really_aggregate_flag:
        print("Generating aggregate")
        t = time.time()
        agg = Aggregate(
                            se=se,
                            stim_parcellation_name=stim_parcellation_name,
                            rec_parcellation_name=rec_parcellation_name,
                            filter_idxs=filter_idxs,
                        )
        agg.save(os.path.join(output_directory_name, agg_name))
        print('Aggregate generation took', time.time() - t, 'seconds')
    else:
        agg = Aggregate.load(os.path.join(output_directory_name, agg_name + '.pkl'))
    return agg

def collapse(se:StimulationEssentials, agg:Aggregate, min_pd, max_pd, z_th = 5)->dict:

    coll = Collapser()
    ##############
    # conditions
    #############
    coll.add_collapse_condition(condcf.max_peak_delay_condition_factory(max_peak_delay=max_pd, z_th=z_th))
    coll.add_collapse_condition(condcf.min_peak_delay_condition_factory(min_peak_delay=min_pd, z_th=z_th))
    ##############
    # conf: extract and collapse functions
    #############
    # To have implantation stats
    coll.add_collapse_conf_entry(ef.extract_name_factory('implantation'), [colf.count_unique_str])
    # To have patient stats
    coll.add_collapse_conf_entry(ef.extract_name_factory('patient'), [colf.count_unique_str])
    # stimulation parameters
    for param in ('intensity', 'pulse_width', 'frequency', 'mode'):
        coll.add_collapse_conf_entry(ef.extract_stim_param_factory(param), colf.standard_functions)
    # From this below, we will get probability
    coll.add_collapse_conf_entry(ef.extract_peak_feature_factory('ampl', z_th=z_th), [colf.N_total])
    # peak charac
    for peak_charac_prefix in stimulation_essentials.constants.PEAK_CHARAC_PREFIXES_LOW:
        coll.add_collapse_conf_entry(ef.extract_peak_feature_factory(peak_charac_prefix, z_th=z_th), colf.standard_functions) #here just amp, peak delays
    # component charac
    for comp_charac_prefix in stimulation_essentials.constants.COMP_CHARAC_PREFIXES_LOW:
        coll.add_collapse_conf_entry(ef.extract_component_feature_factory(comp_charac_prefix), colf.standard_functions) #keep patients data necesary only

    # euclidean distance
    euclidean_distance_extract_function = ef.extract_euclidean_distance_factory('pat')
    coll.add_collapse_conf_entry(euclidean_distance_extract_function, colf.standard_functions)
    # for probability
    coll.add_collapse_conf_entry(ef.extract_peak_feature_factory('ampl', z_th=z_th), [colf.N_total])
    #Peak delay - Already saved due with standard functions for peak characteristics
                # #Mean
                # coll.add_collapse_conf_entry(ef.extract_peak_feature_factory('peak_delay', z_th=z_th), [colf.nanmean])
                # #Median
                # coll.add_collapse_conf_entry(ef.extract_peak_feature_factory('peak_delay', z_th=z_th), [colf.nanquantile_factory(0.5)])


    print('Collapsing...')
    t = time.time()
    collapse_result = coll.collapse(agg, se)
    # saved with collapser method
                # N_recordings = collapse_result['feature_ampl_zth' + str(z_th)]['N_total']
                # N_implantations_unique = collapse_result['implantation_name']['count_unique_str']
                # N_patients_unique = collapse_result['patient_name']['count_unique_str']
                # mean_peak_delay = collapse_result['feature_peak_delay_zth' + str(z_th)]['nanmean']
                # median_peak_delay = collapse_result['feature_peak_delay_zth' + str(z_th)]['nan_median_abs_deviation']

    prob, mask_prob = compute_probability(collapse_result['feature_ampl_zth' + str(z_th)]['N_with_value'], collapse_result['feature_ampl_zth' + str(z_th)]['N_total'], N_tot_min=1)

    # mask the results
    mask_feature = coll.mask_results()
    collapse_result_masked = coll.collapse_result
    print("...it took", time.time() - t)
    print("Collapse result:")
    print(collapse_result_masked)
    return (coll, prob)

