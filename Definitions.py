"""
@author: avalos-alais.s
"""
# Definitions of parcels of interest in Lausanne2008-125
# Definition of network components in Lausanne2008-250
# Paths in use for them and indexes of labels in specific parcellations
import os
from parcellation2.parcellation2 import Parcellation2
from collections import OrderedDict

SE_path = r'F:\FTRACT\SE\2024-01-10-10-52_filter_index_0_fibre_dist_mni_False_parc_reindexed'
#Not specific to stimulation essentials. Just to get parcellation definitions.
parcellation_path = os.path.join(SE_path, 'parcellation_definitions') # All labels parcellations
#Definitions in Lausanne125
#%% All labels parcellations
L_s    =  'Lausanne2008'
res33  = 'Lausanne2008-33'
res125 = 'Lausanne2008-125'
res250 = 'Lausanne2008-250' #FN analysis
res500 = 'Lausanne2008-500'

os.chdir(parcellation_path)
labels_33  = [line.strip() for line in open(res33  + '.txt')]
#Lausanne2008-33 has labels starting with 'ctx-lh-' or 'ctx-rh'. It is the only one and it causes problems for plotting. So correct and save it that way for further use of the parcellation.
labels_33 = [s.replace("ctx-lh-", "lh.").replace("ctx-rh-", "rh.") for s in labels_33]
with open(res33 + '.txt', 'w') as file: file.write('\n'.join(labels_33) + '\n')
labels_125 = [line.strip() for line in open(res125 + '.txt')]
labels_250 = [line.strip() for line in open(res250 + '.txt')]
labels_500 = [line.strip() for line in open(res500 + '.txt')]

#%%LPFC ROIs
#TODO: I could implement a class roi with the parcels of interest and their colors. Then methods to get name and to get color.
roi_dlpfc_125 = [
    "lh.superiorfrontal_6",
    "lh.rostralmiddlefrontal_1",
    "lh.rostralmiddlefrontal_2",
    "lh.rostralmiddlefrontal_3",
    "lh.caudalmiddlefrontal_1",
    "lh.caudalmiddlefrontal_2",
    "lh.caudalmiddlefrontal_3",

    "rh.superiorfrontal_4",
    "rh.superiorfrontal_8",
    "rh.rostralmiddlefrontal_1",
    "rh.rostralmiddlefrontal_2",
    "rh.rostralmiddlefrontal_3",
    "rh.caudalmiddlefrontal_1",
    "rh.caudalmiddlefrontal_2",
    "rh.caudalmiddlefrontal_3",
]
roi_ifg_125 = [
    "lh.parstriangularis_1",
    "lh.parsopercularis_1",
    "lh.parsopercularis_2",

    "rh.parstriangularis_1",
    "rh.parstriangularis_2",
    "rh.parsopercularis_1",
    "rh.parsopercularis_2",
    ]
roi_dlpfc_ifg_125 = [
        # lh DLPFC
        *[i for i in roi_dlpfc_125 if i.startswith("lh.")],

        # lh IFG
        *[i for i in roi_ifg_125 if i.startswith("lh.")],

        # rh DLPFC
        *[i for i in roi_dlpfc_125 if i.startswith("rh.")],

        # rh IFG
        *[i for i in roi_ifg_125 if i.startswith("rh.")],
    ]


colors_dlpfc_ifg_dict = OrderedDict([
        ('lh.superiorfrontal_6',            '#000271'),
        ('lh.rostralmiddlefrontal_1',       '#ff9233'),
        ('lh.rostralmiddlefrontal_2',       '#5233ff'),
        ('lh.rostralmiddlefrontal_3',       '#92f3ff'),
        ('lh.caudalmiddlefrontal_1',        '#ffdd1d'),
        ('lh.caudalmiddlefrontal_2',        '#9c4079'),
        ('lh.caudalmiddlefrontal_3',        '#749e33'),
        ('lh.parstriangularis_1',           '#cc3271'),
        ('lh.parsopercularis_1',            '#ec9fff'),
        ('lh.parsopercularis_2',            '#ff4cfc'),
        ('rh.superiorfrontal_4',            '#ff3731'),
        ('rh.superiorfrontal_8',            '#0e5fff'),
        ('rh.rostralmiddlefrontal_1',       '#2eaff9'),
        ('rh.rostralmiddlefrontal_2',       '#a633ff'),
        ('rh.rostralmiddlefrontal_3',       '#9aee1a'),
        ('rh.caudalmiddlefrontal_1',        '#b88dff'),
        ('rh.caudalmiddlefrontal_2',        '#0f6e0e'),
        ('rh.caudalmiddlefrontal_3',        '#ffd6f5'),
        ('rh.parstriangularis_1',           '#feff33'),
        ('rh.parstriangularis_2',           '#079599'),
        ('rh.parsopercularis_1',            '#ff780e'),
        ('rh.parsopercularis_2',            '#c2161e'),
    ])

#ROI segments for network analysis
roi_ifg_L = [i for i in roi_ifg_125 if i.startswith('lh.')]
roi_ifg_R = [i for i in roi_ifg_125 if i.startswith('rh.')]
#LEFT
roi_ant_dlpfc_L = [
    "lh.rostralmiddlefrontal_1",
    "lh.rostralmiddlefrontal_2",
    "lh.rostralmiddlefrontal_3",]
roi_post_dlpfc_L = [
    "lh.caudalmiddlefrontal_1",
    "lh.caudalmiddlefrontal_2",
    "lh.caudalmiddlefrontal_3",
    "lh.superiorfrontal_6",]
roi_sup_dlpfc_L = [
    "lh.superiorfrontal_6",
    "lh.caudalmiddlefrontal_1",
    "lh.caudalmiddlefrontal_3",
    "lh.rostralmiddlefrontal_2",
    ]
roi_inf_dlpfc_L = [
    "lh.rostralmiddlefrontal_3",
    "lh.rostralmiddlefrontal_1",
    "lh.caudalmiddlefrontal_2",
    ]
#RIGHT
roi_ant_dlpfc_R = [
    "rh.superiorfrontal_4",
    "rh.rostralmiddlefrontal_1",
    "rh.rostralmiddlefrontal_3",]
roi_post_dlpfc_R = [
    "rh.caudalmiddlefrontal_1",
    "rh.caudalmiddlefrontal_2",
    "rh.caudalmiddlefrontal_3",
    "rh.superiorfrontal_8",
    "rh.rostralmiddlefrontal_2",] #check this rostral here
roi_sup_dlpfc_R = [
    "rh.superiorfrontal_4",
    "rh.superiorfrontal_8",
    "rh.caudalmiddlefrontal_3",
    "rh.rostralmiddlefrontal_1",
    ]
roi_inf_dlpfc_R = [
    "rh.rostralmiddlefrontal_3",
    "rh.rostralmiddlefrontal_2",
    "rh.caudalmiddlefrontal_1",
    "rh.caudalmiddlefrontal_2",
    ]
#Idx roi
os.chdir(parcellation_path)
p = Parcellation2(parcellation_path, 'Lausanne2008-125')
index_dlpfc_ifg = [p.get_index_by_name(i) for i in roi_dlpfc_ifg_125]  # used in functional networks
index_dlpfc_ifg_L = [p.get_index_by_name(i) for i in roi_dlpfc_ifg_125 if i.startswith('lh.')]  # used in functional networks
index_dlpfc_ifg_R = [p.get_index_by_name(i) for i in roi_dlpfc_ifg_125 if i.startswith('rh.')]  # used in functional networks
index_dlpfc_L = [p.get_index_by_name(i) for i in roi_dlpfc_125 if i.startswith('lh.')]
index_dlpfc_R = [p.get_index_by_name(i) for i in roi_dlpfc_125 if i.startswith('rh.')]
index_ifg_L = [p.get_index_by_name(i) for i in roi_ifg_L ]
index_ifg_R = [p.get_index_by_name(i) for i in roi_ifg_R ]
#segments
index_ant_dlpfc_L = [p.get_index_by_name(i) for i in roi_ant_dlpfc_L]
index_post_dlpfc_L = [p.get_index_by_name(i) for i in roi_post_dlpfc_L]
index_sup_dlpfc_L = [p.get_index_by_name(i) for i in roi_sup_dlpfc_L]
index_inf_dlpfc_L = [p.get_index_by_name(i) for i in roi_inf_dlpfc_L]

index_ant_dlpfc_R = [p.get_index_by_name(i) for i in roi_ant_dlpfc_R]
index_post_dlpfc_R = [p.get_index_by_name(i) for i in roi_post_dlpfc_R]
index_sup_dlpfc_R = [p.get_index_by_name(i) for i in roi_sup_dlpfc_R]
index_inf_dlpfc_R = [p.get_index_by_name(i) for i in roi_inf_dlpfc_R]
#%% Definition Functional networks Yeo et al.
#Parcels from Lausanne2008-250
functional_networks_left = {
    "lh.Visual": [
        'lh.cuneus_1', 'lh.cuneus_2', 'lh.cuneus_3',
        'lh.lateraloccipital_1', 'lh.lateraloccipital_2', 'lh.lateraloccipital_3',
        'lh.lateraloccipital_4', 'lh.lateraloccipital_5', 'lh.lateraloccipital_6',
        'lh.lateraloccipital_7', 'lh.lateraloccipital_10', 'lh.lateraloccipital_11',
        'lh.lingual_1', 'lh.lingual_2', 'lh.lingual_3', 'lh.lingual_4',
        'lh.lingual_5', 'lh.lingual_6', 'lh.lingual_7', 'lh.lingual_8',
        'lh.pericalcarine_1', 'lh.pericalcarine_2', 'lh.pericalcarine_3',
        'lh.superiorparietal_14', 'lh.superiorparietal_12',
        'lh.fusiform_1', 'lh.fusiform_2', 'lh.fusiform_3', 'lh.fusiform_5',
        'lh.precuneus_11', 'lh.isthmuscingulate_3'
    ],
    "lh.Somatomotor": [
        'lh.postcentral_1','lh.postcentral_2','lh.postcentral_3',
        'lh.postcentral_4','lh.postcentral_5','lh.postcentral_6',
        'lh.postcentral_7','lh.postcentral_8','lh.postcentral_9',
        'lh.postcentral_10','lh.postcentral_11','lh.postcentral_12',
        'lh.postcentral_13','lh.postcentral_14','lh.insula_1','lh.insula_2',
        'lh.insula_3','lh.insula_4','lh.superiortemporal_1','lh.superiortemporal_2',
        'lh.superiortemporal_3','lh.superiortemporal_4','lh.superiortemporal_5',
        'lh.superiortemporal_8','lh.superiortemporal_9','lh.transversetemporal_1',
        'lh.transversetemporal_2','lh.paracentral_1','lh.paracentral_2','lh.paracentral_3',
        'lh.paracentral_4','lh.precuneus_1'
    ],
    "lh.DorsalAttention": [
        'lh.superiorparietal_1','lh.superiorparietal_2','lh.superiorparietal_3',
        'lh.superiorparietal_4','lh.superiorparietal_5','lh.superiorparietal_7',
        'lh.superiorparietal_8','lh.superiorparietal_9','lh.superiorparietal_10',
        'lh.superiorparietal_11','lh.superiorparietal_13','lh.inferiorparietal_1',
        'lh.inferiortemporal_6','lh.inferiortemporal_7','lh.inferiortemporal_8',
        'lh.fusiform_4','lh.precentral_5','lh.precuneus_4','lh.supramarginal_1',
        'lh.supramarginal_4','lh.supramarginal_9','lh.lateraloccipital_8','lh.lateraloccipital_9'
         # 'lh.precentral_4.label',#touches lpfc
         # 'lh.precentral_9.label',#touchs lpfc
         # 'lh.precentral_10.label'#touchs lpfc
    ],
    "lh.VentralAttention": [
        'lh.insula_5','lh.insula_7','lh.bankssts_3','lh.superiorfrontal_8',
        'lh.superiorfrontal_14','lh.superiorfrontal_15','lh.superiorfrontal_16',
        'lh.superiorfrontal_18','lh.supramarginal_2','lh.supramarginal_3',
        'lh.supramarginal_5','lh.supramarginal_6','lh.precuneus_2','lh.paracentral_5'

        # 'lh.insula_6.label'
        # 'lh.precentral_13.label'
        # 'lh.precentral_14.label'
        # 'lh.precentral_15.label'
        # 'lh.precentral_16.label'
        # 'lh.superiorfrontal_17.label'

    ],
    "lh.Limbic": [
        'lh.parahippocampal_1','lh.parahippocampal_2','lh.parahippocampal_3',
        'lh.inferiortemporal_1','lh.inferiortemporal_2','lh.inferiortemporal_3',
        'lh.inferiortemporal_4','lh.middletemporal_7','lh.superiortemporal_11',
        'lh.entorhinal_1','lh.fusiform_6','lh.fusiform_7','lh.fusiform_8','lh.frontalpole_1',
        'lh.lateralorbitofrontal_1','lh.lateralorbitofrontal_2','lh.lateralorbitofrontal_3',
        'lh.lateralorbitofrontal_5','lh.lateralorbitofrontal_6','lh.lateralorbitofrontal_7',
        'lh.medialorbitofrontal_1','lh.medialorbitofrontal_3','lh.medialorbitofrontal_4',
        'lh.medialorbitofrontal_5','lh.temporalpole_1',
    ],
    "lh.FrontoParietal":[
        'lh.caudalanteriorcingulate_1','lh.caudalanteriorcingulate_2','lh.inferiorparietal_8',
        'lh.inferiorparietal_9','lh.inferiorparietal_10','lh.isthmuscingulate_1',
        'lh.posteriorcingulate_1','lh.posteriorcingulate_2','lh.precuneus_5','lh.precuneus_8',
        'lh.superiorfrontal_7','lh.superiorparietal_6','lh.supramarginal_8','lh.supramarginal_10',
        'lh.middletemporal_1','lh.middletemporal_2',

        # 'lh.lateralorbitofrontal_4.label'
        # 'lh.superiorfrontal_9.label'
        # 'lh.superiorfrontal_13.label'
        # 'lh.rostralmiddlefrontal_6.label'
        # 'lh.rostralmiddlefrontal_7.label' #Touchin LPFC

    ],
    "lh.Default": [
        'lh.inferiorparietal_5','lh.inferiorparietal_6','lh.inferiorparietal_7',
        'lh.inferiortemporal_5','lh.middletemporal_3','lh.middletemporal_4',
        'lh.middletemporal_5','lh.middletemporal_6','lh.bankssts_1','lh.bankssts_2',
        'lh.isthmuscingulate_2','lh.medialorbitofrontal_2','lh.posteriorcingulate_3',
        'lh.posteriorcingulate_4','lh.precuneus_3','lh.precuneus_6','lh.precuneus_7',
        'lh.precuneus_9','lh.precuneus_10','lh.rostralanteriorcingulate_1',
        'lh.rostralanteriorcingulate_2','lh.rostralmiddlefrontal_9','lh.rostralmiddlefrontal_11',
        'lh.rostralmiddlefrontal_12','lh.superiorfrontal_1','lh.superiorfrontal_2',
        'lh.superiorfrontal_3','lh.superiorfrontal_4','lh.superiorfrontal_5',
        'lh.superiortemporal_6','lh.superiortemporal_7','lh.superiortemporal_10',
        'lh.supramarginal_7','lh.inferiorparietal_2','lh.inferiorparietal_3','lh.inferiorparietal_4'

        # 'lh.rostralmiddlefrontal_10.label'
        # 'lh.superiorfrontal_6.label'
        # 'lh.superiorfrontal_10.label'
        # 'lh.superiorfrontal_11.label'
        # 'lh.superiorfrontal_12.label'
        # 'lh.parsorbitalis_1.label'
        # 'lh.parsorbitalis_2.label'
        # 'lh.parstriangularis_1.label'
        # 'lh.parstriangularis_2.label'
        # 'lh.rostralmiddlefrontal_8.label' #Touchin LPFC

    ]
}

functional_networks_right = {
    "rh.Visual": [
        'rh.cuneus_1','rh.cuneus_2','rh.cuneus_3','rh.cuneus_4',
        'rh.lateraloccipital_1','rh.lateraloccipital_2','rh.lateraloccipital_3',
        'rh.lateraloccipital_4','rh.lateraloccipital_5','rh.lateraloccipital_6',
        'rh.lateraloccipital_7','rh.lateraloccipital_8','rh.lateraloccipital_9',
        'rh.lateraloccipital_10','rh.lingual_1','rh.lingual_2','rh.lingual_3',
        'rh.lingual_4','rh.lingual_5','rh.lingual_6','rh.lingual_7','rh.pericalcarine_1',
        'rh.pericalcarine_2','rh.pericalcarine_3','rh.pericalcarine_4','rh.superiorparietal_12',
        'rh.superiorparietal_13','rh.fusiform_1','rh.fusiform_2','rh.fusiform_3',
        'rh.fusiform_5', 'rh.precuneus_1','rh.precuneus_2'
    ],
    "rh.Somatomotor": [
        'rh.precentral_3','rh.precentral_6','rh.precentral_7','rh.precentral_8',
        'rh.precentral_11','rh.precentral_15','rh.precentral_16','rh.postcentral_1',
        'rh.postcentral_2','rh.postcentral_3','rh.postcentral_4','rh.postcentral_5',
        'rh.postcentral_6','rh.postcentral_7','rh.postcentral_8','rh.postcentral_9',
        'rh.postcentral_10','rh.postcentral_11','rh.postcentral_12','rh.insula_1','rh.insula_2',
        'rh.insula_4','rh.superiortemporal_1','rh.superiortemporal_2','rh.superiortemporal_4',
        'rh.superiortemporal_5','rh.superiortemporal_8','rh.superiortemporal_9',
        'rh.transversetemporal_1','rh.paracentral_1','rh.paracentral_2','rh.paracentral_3',
        'rh.paracentral_4','rh.paracentral_5','rh.supramarginal_7','rh.supramarginal_9'

    # 'rh.precentral_1.label'
    # 'rh.precentral_2.label'
    # 'rh.insula_3.label' #limit LPFC
    ],
    "rh.DorsalAttention": [
        'rh.superiorparietal_1','rh.superiorparietal_2','rh.superiorparietal_3',
        'rh.superiorparietal_4','rh.superiorparietal_5','rh.superiorparietal_6',
        'rh.superiorparietal_7','rh.superiorparietal_8','rh.superiorparietal_9',
        'rh.superiorparietal_10','rh.superiorparietal_11','rh.inferiortemporal_6',
        'rh.inferiortemporal_7','rh.fusiform_4','rh.precentral_12','rh.precentral_14',
        'rh.precuneus_8','rh.precuneus_9','rh.inferiorparietal_11','rh.inferiorparietal_12',
        'rh.supramarginal_1','rh.supramarginal_2'

        # 'rh.precentral_4.label'
        # 'rh.precentral_5.label'
        # 'rh.precentral_9.label'
        # 'rh.precentral_10.label'
        # 'rh.precentral_13.label'
    ],
    "rh.VentralAttention": [
        'rh.insula_5','rh.insula_6','rh.superiorfrontal_11','rh.superiorfrontal_16',
        'rh.supramarginal_3','rh.supramarginal_4','rh.supramarginal_5','rh.supramarginal_6',
        'rh.supramarginal_8','rh.paracentral_6','rh.middletemporal_1','rh.bankssts_1',
        'rh.bankssts_2','rh.bankssts_3','rh.caudalanteriorcingulate_2'

        # 'rh.insula_7.label'
        # 'rh.superiorfrontal_8.label'
    ],
    "rh.Limbic":[
        'rh.parahippocampal_1','rh.parahippocampal_2', 'rh.parahippocampal_3',
        'rh.inferiortemporal_1','rh.inferiortemporal_2','rh.inferiortemporal_3',
        'rh.inferiortemporal_4','rh.superiortemporal_11','rh.entorhinal_1',
        'rh.fusiform_6','rh.fusiform_7','rh.fusiform_8', 'rh.frontalpole_1',
        'rh.lateralorbitofrontal_1','rh.lateralorbitofrontal_3','rh.lateralorbitofrontal_5',
        'rh.lateralorbitofrontal_6','rh.lateralorbitofrontal_7', 'rh.medialorbitofrontal_3',
        'rh.medialorbitofrontal_4','rh.medialorbitofrontal_5','rh.temporalpole_1',
        'rh.middletemporal_8', 'rh.middletemporal_9','rh.superiortemporal_10',
        # 'rh.lateralorbitofrontal_2', limit lpfc
    ],
    "rh.FrontoParietal" : [
        'rh.caudalanteriorcingulate_3','rh.lateralorbitofrontal_4', 'rh.posteriorcingulate_1',
        'rh.posteriorcingulate_2','rh.posteriorcingulate_3','rh.precuneus_5',
        'rh.middletemporal_2','rh.rostralmiddlefrontal_12','rh.rostralmiddlefrontal_13',
        'rh.inferiortemporal_5','rh.inferiorparietal_1','rh.inferiorparietal_2',
        'rh.inferiorparietal_3','rh.inferiorparietal_5','rh.inferiorparietal_6',
        # 'rh.rostralmiddlefrontal_6', dlpfc
        # 'rh.rostralmiddlefrontal_7',limite dlpfc
        # 'rh.rostralmiddlefrontal_8', limite dlpfc
        # 'rh.rostralmiddlefrontal_9', limite dlpfc
        # 'rh.superiorfrontal_7', limite dlpfc
        # 'rh.superiorfrontal_13', limite dlpfc
        # 'rh.superiorfrontal_14', limite dlpfc
        # 'rh.superiorfrontal_15', limite dlpfc
        # 'rh.superiorfrontal_17', limite dlpfc
    ],

    "rh.Default": [
        'rh.inferiorparietal_4', 'rh.inferiorparietal_7','rh.middletemporal_3',
        'rh.middletemporal_4','rh.middletemporal_5','rh.middletemporal_6',
        'rh.middletemporal_7','rh.isthmuscingulate_1', 'rh.isthmuscingulate_2',
        'rh.medialorbitofrontal_1','rh.medialorbitofrontal_2','rh.posteriorcingulate_4',
        'rh.precuneus_3','rh.precuneus_4','rh.precuneus_6','rh.precuneus_7',
        'rh.precuneus_10','rh.caudalanteriorcingulate_1','rh.rostralanteriorcingulate_1',
        'rh.rostralanteriorcingulate_2','rh.rostralmiddlefrontal_11','rh.superiorfrontal_1',
        'rh.superiorfrontal_2','rh.superiorfrontal_3','rh.superiorfrontal_5',
        'rh.superiortemporal_3', 'rh.superiortemporal_6','rh.superiortemporal_7',
        'rh.inferiorparietal_8','rh.inferiorparietal_9', 'rh.inferiorparietal_10',

        # 'rh.rostralmiddlefrontal_10', limite dlpfc
        # 'rh.superiorfrontal_4.label',
        # 'rh.superiorfrontal_6.label', limite dlpfc
        # 'rh.superiorfrontal_9.label', limite dlpfc
        # 'rh.superiorfrontal_10.label', limite dlpfc
        # 'rh.superiorfrontal_12.label', limite dlpfc
        # 'rh.parsorbitalis_1.label', limite dlpfc
        # 'rh.parsorbitalis_2.label', limite dlpfc
        # 'rh.parstriangularis_1.label', limite dlpfc
        # 'rh.parstriangularis_2.label',


    ]
}
p = Parcellation2(parcellation_path, 'Lausanne2008-250')
fn_L_labels = list(functional_networks_left.keys())
idx_fn_L = [tuple(p.get_index_by_name(parcel) for parcel in parcel_list)
            for parcel_list in functional_networks_left.values()]
fn_R_labels = list(functional_networks_right.keys())
idx_fn_R = [tuple(p.get_index_by_name(parcel) for parcel in parcel_list)
            for parcel_list in functional_networks_right.values()]
colors_fn = ['#89259a', '#64aff3', '#2ecd00', '#f566d6', '#ffe496', '#ff972b','#ff1010']  # purple,blue, green, pink, yellow, orange, red


#%% For Plotting, Labels paths
labels_path = r'C:\Users\avalos-alais.s\ft_LPFC\FTRACT_LPFC\MNE-data\subjects\cvs_avg35_inMNI152'
labels_path_125 = os.path.join(labels_path + '-Lausanne125', 'label')  # for plots
labels_path_33 = os.path.join(labels_path + '-Lausanne33', 'label')  # for plots
labels_path_500 = os.path.join(labels_path + '-Lausanne500', 'label')  # for resolutions
labels_path_250 = os.path.join(labels_path + '-Lausanne250', 'label')  # for networks
meshdirname = labels_path #Default in plotting functions, change parameter direction for implementation