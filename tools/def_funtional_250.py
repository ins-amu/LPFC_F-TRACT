# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 17:00:58 2024

@author: avalos-alais.s
"""

#DEF FUNCTIONAL NETWORKS (merge 250) to DLPFC 125

import os
def functional_networks_250( parcel_defs, plot_flag = False) : 
    os.chdir( parcel_defs)
    labels = open('Lausanne2008-250.txt','r')
    labels =labels.read()
    labels = labels.split()
    labels_list = [s + '.label' for s in labels] 
    
    labels_dlpfc = open('Lausanne2008-125.txt','r')
    labels_dlpfc = labels_dlpfc.read()
    labels_dlpfc = labels_dlpfc.split()
 
    #% Networ definitions : Visual; Somatomotor; Dorsal Attention ; Ventral Attention ; Limbic ; Frontoparietal ; Default
    
    Visual_L = [
                'lh.cuneus_1.label',
                'lh.cuneus_2.label',
                'lh.cuneus_3.label',
                'lh.lateraloccipital_1.label',
                'lh.lateraloccipital_2.label',
                'lh.lateraloccipital_3.label',
                'lh.lateraloccipital_4.label',
                'lh.lateraloccipital_5.label',
                'lh.lateraloccipital_6.label',
                'lh.lateraloccipital_7.label',
                'lh.lateraloccipital_10.label',
                'lh.lateraloccipital_11.label',
                
     
                'lh.lingual_1.label',
                'lh.lingual_2.label',
                'lh.lingual_3.label',
                'lh.lingual_4.label',
                'lh.lingual_5.label',
                'lh.lingual_6.label',
                'lh.lingual_7.label',
                'lh.lingual_8.label',
                'lh.pericalcarine_1.label',
                'lh.pericalcarine_2.label',
                'lh.pericalcarine_3.label',
                'lh.superiorparietal_14.label',
                'lh.superiorparietal_12.label',
                
                'lh.fusiform_1.label',
                'lh.fusiform_2.label',
                'lh.fusiform_3.label',
                
                'lh.fusiform_5.label',
            
                'lh.precuneus_11.label',
                'lh.isthmuscingulate_3.label', #or visual
    ]
    Somatomotor_L = [
                    'lh.precentral_1.label',
                    'lh.precentral_2.label',
                    'lh.precentral_3.label',
                    'lh.precentral_6.label',
                    'lh.precentral_7.label',
                    'lh.precentral_8.label',
                    'lh.precentral_11.label',
                    'lh.precentral_12.label',
                    
                  
                    
                    'lh.postcentral_1.label',
                    'lh.postcentral_2.label',
                    'lh.postcentral_3.label',
                    'lh.postcentral_4.label',
                    'lh.postcentral_5.label',
                    'lh.postcentral_6.label',
                    'lh.postcentral_7.label',
                    'lh.postcentral_8.label',
                    'lh.postcentral_9.label',
                    'lh.postcentral_10.label',
                    'lh.postcentral_11.label',
                    'lh.postcentral_12.label',
                    'lh.postcentral_13.label',
                    'lh.postcentral_14.label',
                    
                    'lh.insula_1.label',
                    'lh.insula_2.label',
                    'lh.insula_3.label',
                    'lh.insula_4.label',
                  
                  
                    'lh.superiortemporal_1.label',
                    'lh.superiortemporal_2.label',
                    'lh.superiortemporal_3.label',
                    'lh.superiortemporal_4.label',
                    'lh.superiortemporal_5.label',
                    'lh.superiortemporal_8.label',
                    'lh.superiortemporal_9.label',
                    
                    'lh.transversetemporal_1.label',
                    'lh.transversetemporal_2.label',
                  
                    'lh.paracentral_1.label',
                    'lh.paracentral_2.label',
                    'lh.paracentral_3.label',
                    'lh.paracentral_4.label',
                    
                 
                    'lh.precuneus_1.label',
              

        ]
    
    Dorsal_Attention_L = [
        'lh.superiorparietal_1.label',
        'lh.superiorparietal_2.label',
        'lh.superiorparietal_3.label',
        'lh.superiorparietal_4.label',
        'lh.superiorparietal_5.label',
        'lh.superiorparietal_7.label',
        'lh.superiorparietal_8.label',
        'lh.superiorparietal_9.label',
        'lh.superiorparietal_10.label',
        'lh.superiorparietal_11.label',
        'lh.superiorparietal_13.label',
    
        'lh.inferiorparietal_1.label',
       
        
        'lh.inferiortemporal_6.label',
        'lh.inferiortemporal_7.label',
        'lh.inferiortemporal_8.label',
        
        
        'lh.fusiform_4.label',
        
        
        # 'lh.precentral_4.label',#touches dlpfc
        'lh.precentral_5.label',
        # 'lh.precentral_9.label',#touchs dlpfc
        # 'lh.precentral_10.label'#touchs dlpfc,
        
        'lh.precuneus_4.label',
        
        'lh.supramarginal_1.label',
        'lh.supramarginal_4.label',
        'lh.supramarginal_9.label',
        'lh.lateraloccipital_8.label',
        'lh.lateraloccipital_9.label',
    
        
        ]
    Ventral_Attention_L = [
        'lh.insula_5.label',
        # 'lh.insula_6.label', #touchs dlpfc
        'lh.insula_7.label',
        
        
        'lh.bankssts_3.label',
        
        
        # 'lh.precentral_13.label',#touchs dlpfc
        # 'lh.precentral_14.label',#touchs dlpfc
        # 'lh.precentral_15.label',#touchs dlpfc
        # 'lh.precentral_16.label',#touchs dlpfc
        
         
        'lh.superiorfrontal_8.label',
        'lh.superiorfrontal_14.label',
        'lh.superiorfrontal_15.label',
        'lh.superiorfrontal_16.label',
        # 'lh.superiorfrontal_17.label',#touches dlpfc
        'lh.superiorfrontal_18.label',
        
        'lh.supramarginal_2.label',
        'lh.supramarginal_3.label',
        'lh.supramarginal_5.label',
        'lh.supramarginal_6.label',
         
        'lh.precuneus_2.label',
        'lh.paracentral_5.label',
        ]
    Limbic_L = [
       
        'lh.parahippocampal_1.label',
        'lh.parahippocampal_2.label',
        'lh.parahippocampal_3.label',
        
        'lh.inferiortemporal_1.label',
        'lh.inferiortemporal_2.label',
        'lh.inferiortemporal_3.label',
        'lh.inferiortemporal_4.label',
        
        
        'lh.middletemporal_7.label',
       
        'lh.superiortemporal_11.label',
            
        'lh.entorhinal_1.label',
        
        'lh.fusiform_6.label',
        'lh.fusiform_7.label',
        'lh.fusiform_8.label',
        
        'lh.frontalpole_1.label',
        
        'lh.lateralorbitofrontal_1.label',
        'lh.lateralorbitofrontal_2.label',
        'lh.lateralorbitofrontal_3.label',
        'lh.lateralorbitofrontal_5.label',
        'lh.lateralorbitofrontal_6.label',
        'lh.lateralorbitofrontal_7.label', #no veo
        
        'lh.medialorbitofrontal_1.label',
        'lh.medialorbitofrontal_3.label',
        'lh.medialorbitofrontal_4.label',
        'lh.medialorbitofrontal_5.label',
        'lh.temporalpole_1.label',
       
        
        ]
    Frontoparietal_L = [
           
     
           
            
            
            'lh.caudalanteriorcingulate_1.label',
            'lh.caudalanteriorcingulate_2.label',
            
            'lh.inferiorparietal_8.label',
            'lh.inferiorparietal_9.label',
            'lh.inferiorparietal_10.label',
            
            'lh.isthmuscingulate_1.label',
            
            # 'lh.lateralorbitofrontal_4.label', #touches dlpfc
            
            'lh.posteriorcingulate_1.label',
            'lh.posteriorcingulate_2.label',
            
            'lh.precuneus_5.label',
            'lh.precuneus_8.label',
            
            'lh.superiorfrontal_7.label',
            # 'lh.superiorfrontal_9.label',#touches dlpfc
            # 'lh.superiorfrontal_13.label'#touches dlpfc,
            
            'lh.superiorparietal_6.label',
            
            'lh.supramarginal_8.label',
            'lh.supramarginal_10.label',
            
            'lh.middletemporal_1.label', 
            'lh.middletemporal_2.label',
            
            # 'lh.rostralmiddlefrontal_6.label',#touchs dlpfc
            # 'lh.rostralmiddlefrontal_7.label',#touchs dlpfc
          #estas estan pegadas a la dlpfc podria ser default
            
            
        ]
    Default_L = [
        
        'lh.inferiorparietal_5.label',
        'lh.inferiorparietal_6.label',
        'lh.inferiorparietal_7.label',
        
        
        'lh.inferiortemporal_5.label',
        'lh.middletemporal_3.label',
        'lh.middletemporal_4.label',
        'lh.middletemporal_5.label',
        'lh.middletemporal_6.label',
        
        'lh.bankssts_1.label',
        'lh.bankssts_2.label',
        
        'lh.isthmuscingulate_2.label',
    
        
        'lh.medialorbitofrontal_2.label',
        
        'lh.posteriorcingulate_3.label',
        'lh.posteriorcingulate_4.label',
        
        'lh.precuneus_3.label',
        'lh.precuneus_6.label',
        'lh.precuneus_7.label',
        'lh.precuneus_9.label',
        'lh.precuneus_10.label',
        
        'lh.rostralanteriorcingulate_1.label',
        'lh.rostralanteriorcingulate_2.label',
        
        'lh.rostralmiddlefrontal_9.label',
        # 'lh.rostralmiddlefrontal_10.label',#touches dlpfc
        'lh.rostralmiddlefrontal_11.label',
        'lh.rostralmiddlefrontal_12.label',
        
        'lh.superiorfrontal_1.label',
        'lh.superiorfrontal_2.label',
        'lh.superiorfrontal_3.label',
        'lh.superiorfrontal_4.label',
        'lh.superiorfrontal_5.label',
        # 'lh.superiorfrontal_6.label',#touches dlpfc
        
        # 'lh.superiorfrontal_10.label', #touches dlpfc
        # 'lh.superiorfrontal_11.label', #touches dlpfc
        # 'lh.superiorfrontal_12.label', #touches dlpfc
     
        
        'lh.superiortemporal_6.label',
        'lh.superiortemporal_7.label',
        'lh.superiortemporal_10.label',
        
        'lh.supramarginal_7.label',
        'lh.inferiorparietal_2.label',
        'lh.inferiorparietal_3.label',
        'lh.inferiorparietal_4.label',
        
        # 'lh.parsorbitalis_1.label',#touches dlpfc
        # 'lh.parsorbitalis_2.label',#touches dlpfc
        # 'lh.parstriangularis_1.label',  #touchs dlpfc
        # 'lh.parstriangularis_2.label',
        # 'lh.rostralmiddlefrontal_8.label', 
        
        ]
    
    
    #RH
    Visual_R = [
    'rh.cuneus_1.label',
    'rh.cuneus_2.label',
    'rh.cuneus_3.label',
    'rh.cuneus_4.label',
    'rh.lateraloccipital_1.label',
    'rh.lateraloccipital_2.label',
    'rh.lateraloccipital_3.label',
    'rh.lateraloccipital_4.label',
    'rh.lateraloccipital_5.label',
    'rh.lateraloccipital_6.label',
    'rh.lateraloccipital_7.label',
    'rh.lateraloccipital_8.label',
    'rh.lateraloccipital_9.label',
    'rh.lateraloccipital_10.label',
   
    'rh.lingual_1.label',
    'rh.lingual_2.label',
    'rh.lingual_3.label',
    'rh.lingual_4.label',
    'rh.lingual_5.label',
    'rh.lingual_6.label',
    'rh.lingual_7.label',

    'rh.pericalcarine_1.label',
    'rh.pericalcarine_2.label',
    'rh.pericalcarine_3.label',
    'rh.pericalcarine_4.label',
    'rh.superiorparietal_12.label',
    'rh.superiorparietal_13.label',
    'rh.fusiform_1.label',
    'rh.fusiform_2.label',
    'rh.fusiform_3.label',
    'rh.fusiform_5.label',
    'rh.precuneus_1.label',
    'rh.precuneus_2.label',
]

    Somatomotor_R = [
        # 'rh.precentral_1.label', limite dlpfc
        # 'rh.precentral_2.label', limite dlpfc
        'rh.precentral_3.label',
        'rh.precentral_6.label',
        'rh.precentral_7.label',
        'rh.precentral_8.label',
        'rh.precentral_11.label',
        'rh.precentral_15.label',
        'rh.precentral_16.label',
        
        'rh.postcentral_1.label',
        'rh.postcentral_2.label',
        'rh.postcentral_3.label',
        'rh.postcentral_4.label',
        'rh.postcentral_5.label',
        'rh.postcentral_6.label',
        'rh.postcentral_7.label',
        'rh.postcentral_8.label',
        'rh.postcentral_9.label',
        'rh.postcentral_10.label',
        'rh.postcentral_11.label',
        'rh.postcentral_12.label',
        'rh.insula_1.label',
        'rh.insula_2.label',
        # 'rh.insula_3.label', limite dlpfc
        'rh.insula_4.label',
        'rh.superiortemporal_1.label',
        'rh.superiortemporal_2.label',
        'rh.superiortemporal_4.label',
        'rh.superiortemporal_5.label',
        'rh.superiortemporal_8.label',
        'rh.superiortemporal_9.label',
        'rh.transversetemporal_1.label',
        'rh.paracentral_1.label',
        'rh.paracentral_2.label',
        'rh.paracentral_3.label',
        'rh.paracentral_4.label',
        'rh.paracentral_5.label',
        
        
        'rh.supramarginal_7.label',
        'rh.supramarginal_9.label',
    ]
    
    Dorsal_Attention_R = [
        'rh.superiorparietal_1.label',
        'rh.superiorparietal_2.label',
        'rh.superiorparietal_3.label',
        'rh.superiorparietal_4.label',
        'rh.superiorparietal_5.label',
        'rh.superiorparietal_6.label',
        'rh.superiorparietal_7.label',
        'rh.superiorparietal_8.label',
        'rh.superiorparietal_9.label',
        'rh.superiorparietal_10.label',
        'rh.superiorparietal_11.label',
        
       
        'rh.inferiortemporal_6.label',
        'rh.inferiortemporal_7.label',
        'rh.fusiform_4.label',
        # 'rh.precentral_4.label', limite dlpfc
        # 'rh.precentral_5.label', limite dlpfc
        # 'rh.precentral_9.label', limite dlpfc 
        # 'rh.precentral_10.label',limite dlpfc
        'rh.precentral_12.label',
        # 'rh.precentral_13.label', limite dlpfc 
        'rh.precentral_14.label',

        
        'rh.precuneus_8.label',
        'rh.precuneus_9.label',
        
        
        'rh.inferiorparietal_11.label',
        'rh.inferiorparietal_12.label',
        'rh.supramarginal_1.label',
        'rh.supramarginal_2.label',
        
    ]
    
    Ventral_Attention_R = [
        'rh.insula_5.label',
        'rh.insula_6.label',
        # 'rh.insula_7.label', limite dlpfc

        # 'rh.superiorfrontal_8.label',  limite dlpfc
        'rh.superiorfrontal_11.label',
        'rh.superiorfrontal_16.label',
        
        
        'rh.supramarginal_3.label',
        'rh.supramarginal_4.label',
        'rh.supramarginal_5.label',
        'rh.supramarginal_6.label',
        'rh.supramarginal_8.label',
        
        
        'rh.paracentral_6.label',
    
        'rh.middletemporal_1.label',
        
        'rh.bankssts_1.label',
        'rh.bankssts_2.label',        
        'rh.bankssts_3.label',
        
        'rh.caudalanteriorcingulate_2.label',
    ]
    
    Limbic_R = [
        'rh.parahippocampal_1.label',
        'rh.parahippocampal_2.label',
        'rh.parahippocampal_3.label',
        'rh.inferiortemporal_1.label',
        'rh.inferiortemporal_2.label',
        'rh.inferiortemporal_3.label',
        'rh.inferiortemporal_4.label',
        
        'rh.superiortemporal_11.label',
        'rh.entorhinal_1.label',
        'rh.fusiform_6.label',
        'rh.fusiform_7.label',
        'rh.fusiform_8.label',
        'rh.frontalpole_1.label',
        'rh.lateralorbitofrontal_1.label',
        # 'rh.lateralorbitofrontal_2.label', limite dlpfc
        'rh.lateralorbitofrontal_3.label',
        'rh.lateralorbitofrontal_5.label',
        'rh.lateralorbitofrontal_6.label',
        'rh.lateralorbitofrontal_7.label', #no veo
        
        'rh.medialorbitofrontal_3.label',
        'rh.medialorbitofrontal_4.label',
        'rh.medialorbitofrontal_5.label',
        'rh.temporalpole_1.label',
        
        'rh.middletemporal_8.label',
        'rh.middletemporal_9.label',
        
        'rh.superiortemporal_10.label',
    ]
    
    Frontoparietal_R = [
        'rh.caudalanteriorcingulate_3.label',
      
       
        'rh.lateralorbitofrontal_4.label',
        
        'rh.posteriorcingulate_1.label',
        'rh.posteriorcingulate_2.label',
        'rh.posteriorcingulate_3.label',
        
        'rh.precuneus_5.label',

        
        'rh.middletemporal_2.label',
        # 'rh.rostralmiddlefrontal_6.label', dlpfc
         # 'rh.rostralmiddlefrontal_7.label',limite dlpfc
        # 'rh.rostralmiddlefrontal_8.label', limite dlpfc
        # 'rh.rostralmiddlefrontal_9.label', limite dlpfc
        'rh.rostralmiddlefrontal_12.label',
        'rh.rostralmiddlefrontal_13.label',
        
        'rh.inferiortemporal_5.label',
        'rh.inferiorparietal_1.label',
        'rh.inferiorparietal_2.label',
        'rh.inferiorparietal_3.label',
        'rh.inferiorparietal_5.label',
        'rh.inferiorparietal_6.label',
        
        # 'rh.superiorfrontal_7.label', limite dlpfc
        # 'rh.superiorfrontal_13.label', limite dlpfc
        # 'rh.superiorfrontal_14.label', limite dlpfc
        # 'rh.superiorfrontal_15.label', limite dlpfc
        # 'rh.superiorfrontal_17.label', limite dlpfc 
    ]
    
    Default_R = [
        
        'rh.inferiorparietal_4.label',
        'rh.inferiorparietal_7.label',
        
       
        'rh.middletemporal_3.label',
        'rh.middletemporal_4.label',
        'rh.middletemporal_5.label',
        'rh.middletemporal_6.label',
        'rh.middletemporal_7.label',
        
        'rh.isthmuscingulate_1.label',
        'rh.isthmuscingulate_2.label',
        
        'rh.medialorbitofrontal_1.label',
        'rh.medialorbitofrontal_2.label',
        
        
        'rh.posteriorcingulate_4.label',
        
        'rh.precuneus_3.label',
        'rh.precuneus_4.label',
        'rh.precuneus_6.label',
        'rh.precuneus_7.label',
        
        'rh.precuneus_10.label',
        
        'rh.caudalanteriorcingulate_1.label',
        
        'rh.rostralanteriorcingulate_1.label',
        'rh.rostralanteriorcingulate_2.label',
        
        # 'rh.rostralmiddlefrontal_10.label', limite dlpfc
        'rh.rostralmiddlefrontal_11.label',
       
        
        'rh.superiorfrontal_1.label',
        'rh.superiorfrontal_2.label',
        'rh.superiorfrontal_3.label',
        # 'rh.superiorfrontal_4.label',
        'rh.superiorfrontal_5.label',
        # 'rh.superiorfrontal_6.label', limite dlpfc
        # 'rh.superiorfrontal_9.label', limite dlpfc
        
        # 'rh.superiorfrontal_10.label', limite dlpfc
        # 'rh.superiorfrontal_12.label', limite dlpfc
        
        
        'rh.superiortemporal_3.label',
        'rh.superiortemporal_6.label',
        'rh.superiortemporal_7.label',
        
        # 

        
        'rh.inferiorparietal_8.label',
        'rh.inferiorparietal_9.label', 
        'rh.inferiorparietal_10.label',
        
        # 'rh.parsorbitalis_1.label', limite dlpfc
        # 'rh.parsorbitalis_2.label', limite dlpfc 
        # 'rh.parstriangularis_1.label', limite dlpfc 
        # 'rh.parstriangularis_2.label',
        
        
        ]
    

    #%% DLPFC
    
    dlpfc_roi = [
        "lh.caudalmiddlefrontal_1",
        "lh.caudalmiddlefrontal_2",
        "lh.caudalmiddlefrontal_3",
        "lh.parsopercularis_1",
        "lh.parsopercularis_2",
        "lh.parstriangularis_1",
        "lh.rostralmiddlefrontal_1",
        "lh.rostralmiddlefrontal_2",
        "lh.rostralmiddlefrontal_3",
        "lh.superiorfrontal_6",
        
        "rh.superiorfrontal_4",
        "rh.superiorfrontal_8",
        "rh.rostralmiddlefrontal_1",
        "rh.rostralmiddlefrontal_2",
        "rh.rostralmiddlefrontal_3",
        "rh.caudalmiddlefrontal_1",
        "rh.caudalmiddlefrontal_2",
        "rh.caudalmiddlefrontal_3",
        "rh.parstriangularis_1",
        "rh.parstriangularis_2",
        "rh.parsopercularis_1",
        "rh.parsopercularis_2"
        ]
    index_dlpfc = [labels_dlpfc.index(l) for l in dlpfc_roi]

    
    names_labels_lobes = [ 'lh.Visual', 'lh.Somatomotor', 'lh.DorsalAttention', 'lh.VentralAttention' , 'lh.Limbic', 'lh.Frontoparietal', 'lh.Default' ,
                            'rh.Visual', 'rh.Somatomotor', 'rh.DorsalAttention', 'rh.VentralAttention' , 'rh.Limbic', 'rh.Frontoparietal', 'rh.Default' 
                            ]

    #%%
    index_Visual_L = []
    for l in Visual_L:
        index_Visual_L.append(labels_list.index(l))
    index_Visual_L=tuple(index_Visual_L)
    index_Somatomotor_L = []
    for l in Somatomotor_L:
        index_Somatomotor_L.append(labels_list.index(l))
    index_Somatomotor_L=tuple(index_Somatomotor_L)
    index_DorsalAttention_L = []
    for l in Dorsal_Attention_L:
        index_DorsalAttention_L.append(labels_list.index(l))
    index_DorsalAttention_L=tuple(index_DorsalAttention_L)
    index_VentralAttention_L = []
    for l in Ventral_Attention_L:
        index_VentralAttention_L.append(labels_list.index(l))
    index_VentralAttention_L=tuple(index_VentralAttention_L)
    index_Limbic_L = []
    for l in Limbic_L:
        index_Limbic_L.append(labels_list.index(l))
    index_Limbic_L=tuple(index_Limbic_L)
    index_Frontoparietal_L = []
    for l in Frontoparietal_L:
        index_Frontoparietal_L.append(labels_list.index(l))
    index_Frontoparietal_L=tuple(index_Frontoparietal_L)
    index_Default_L = []
    for l in Default_L:
        index_Default_L.append(labels_list.index(l))
    index_Default_L=tuple(index_Default_L)
    
    
    
    index_Visual_R = []
    for l in Visual_R:
        index_Visual_R.append(labels_list.index(l))
    index_Visual_R = tuple(index_Visual_R)

    index_Somatomotor_R = []
    for l in Somatomotor_R:
        index_Somatomotor_R.append(labels_list.index(l))
    index_Somatomotor_R = tuple(index_Somatomotor_R)
    
    index_DorsalAttention_R = []
    for l in Dorsal_Attention_R:
        index_DorsalAttention_R.append(labels_list.index(l))
    index_DorsalAttention_R = tuple(index_DorsalAttention_R)
    
    index_VentralAttention_R = []
    for l in Ventral_Attention_R:
        index_VentralAttention_R.append(labels_list.index(l))
    index_VentralAttention_R = tuple(index_VentralAttention_R)
    
    index_Limbic_R = []
    for l in Limbic_R:
        index_Limbic_R.append(labels_list.index(l))
    index_Limbic_R = tuple(index_Limbic_R)
    
    index_Frontoparietal_R = []
    for l in Frontoparietal_R:
        index_Frontoparietal_R.append(labels_list.index(l))
    index_Frontoparietal_R = tuple(index_Frontoparietal_R)
    
    index_Default_R = []
    for l in Default_R:
        index_Default_R.append(labels_list.index(l))
    index_Default_R = tuple(index_Default_R)


    #%% Index
    index_n = [index_Visual_L, index_Somatomotor_L, index_DorsalAttention_L,  index_VentralAttention_L, index_Limbic_L, index_Frontoparietal_L, index_Default_L,  
               index_Visual_R, index_Somatomotor_R, index_DorsalAttention_R,  index_VentralAttention_R, index_Limbic_R, index_Frontoparietal_R, index_Default_R
             ]

    return  names_labels_lobes, index_n , index_dlpfc 

