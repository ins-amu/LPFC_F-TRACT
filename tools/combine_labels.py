# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 17:00:26 2024

@author: avalos-alais.s
"""

#Combine label objects, save the combined label in indicated path. Return new label object
# labels2combine : list of label names ex : [['lh.caudalmiddlefrontal_1.label', 'lh.frontalpole_1.label',] 
# labels_path : path to folder containing label objects to combine 

import mne
import os
def combine_labels(labels2combine, labels_path, new_label_name = None,  labels_output_path = None ):
    
    labels2combine = [label if label.endswith('.label') else label + '.label' for label in labels2combine]
 
    os.chdir(labels_path)
    if len(labels2combine) >= 2 :
        label_combined = mne.read_label(labels2combine[0]).__add__(mne.read_label(labels2combine[1]))
        for l in range(2,len(labels2combine)) :
           label_combined =  label_combined.__add__(mne.read_label(labels2combine[l]))
           
        if new_label_name is not None and labels_output_path is not None:
             os.makedirs(labels_output_path, exist_ok=True)  
             output_file = os.path.join(labels_output_path, f"{new_label_name}.label")
             label_combined.save(output_file)  # Save the combined label
             print(f"Combined label saved as {output_file}")
    
        return label_combined
    elif len(labels2combine) == 1 :
        label = mne.read_label(labels2combine[0])
        output_file = os.path.join(labels_output_path, f"{new_label_name}.label")
        label.save(output_file)
        print('No combination was done, just one label in the list to combine')
        return label
    else: 
        print('Empty list, no labels to combine')
  