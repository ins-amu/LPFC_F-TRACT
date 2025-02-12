# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:36:00 2025

@author: avalos-alais.s
"""

#Index management
import os 
import numpy as np
def get_idx (path_list, text_labels, labels) :
    
    os.chdir(path_list)
    with open(text_labels, 'r') as file:  
        labels_list= [line.strip() for line in file]
        
    index = [labels_list.index(x) for x in labels]
    return index 

def save_idx_txt(data, filename):
    if isinstance(data, np.ndarray):  # Case where data is a NumPy array
        np.savetxt(filename, data)
    elif isinstance(data, (list, tuple)):  # Case where data is a list or tuple
        with open(filename, 'w') as file:
            if isinstance(data[0], (list, tuple)):  # List of tuples/lists
                file.writelines(f"{tup}\n" for tup in data)
            else:  # Just a list or a single tuple
                file.write(f"{data}\n")
                