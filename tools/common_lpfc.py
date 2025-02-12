# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 17:12:15 2024

@author: avalos-alais.s
"""

from collections import OrderedDict
import numpy as np
import os


colors_lpfc_dict = OrderedDict([
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

colors_lpfc = list(colors_lpfc_dict.values())
names_lpfc  = list(colors_lpfc_dict.keys())

def idx_reader(filename):
    with open(filename, 'r') as f:
        r = [line.strip().replace('(', '').replace(')', '').replace(',', '') for line in f]
    return np.array(r).astype(int)

def read_data():

    names_path = r'C:\Users\avalos-alais.s\ft_dti_integration/parcellation/parcellation_definitions'
    all_names_33 = np.loadtxt(os.path.join(names_path, 'Lausanne2008-33.txt'), dtype=str)
    all_names_125 = np.loadtxt(os.path.join(names_path, 'Lausanne2008-125.txt'), dtype=str)
    root_path = r'F:\Results\Res_30oct'
    s = 'Lausanne2008-125__Lausanne2008-33/{}_0_100ms.txt'
    N_lpfc =  np.loadtxt(os.path.join(root_path, s.format('N_125to125')))
    p_dlpfc =  np.loadtxt(os.path.join(root_path, s.format('p_125to125')))
    N_125_33 = np.loadtxt(os.path.join(root_path, s.format('N_125to33')))
    p_125_33 = np.loadtxt(os.path.join(root_path, s.format('p_125to33')))
    CI_125_33 = np.loadtxt(os.path.join(root_path, s.format('CI_125to33')))
    index_x_125_33 = idx_reader(os.path.join(root_path, s.format('index_x_125to33')))
    names_125 = all_names_125[index_x_125_33]
    index_y_125_33 = idx_reader(os.path.join(root_path, s.format('index_y_125to33')))
    names_33 = all_names_33[index_y_125_33]
    s = 'Lausanne2008-33__Lausanne2008-125/{}_0_100ms.txt'
    N_33_125 = np.loadtxt(os.path.join(root_path, s.format('N_33to125')))
    p_33_125 = np.loadtxt(os.path.join(root_path, s.format('p_33to125')))
    CI_33_125 = np.loadtxt(os.path.join(root_path, s.format('CI_33to125')))
    return names_125, names_33, p_125_33, N_125_33, CI_125_33, p_33_125, N_33_125, CI_33_125
    
