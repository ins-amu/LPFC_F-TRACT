# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 16:57:51 2024

@author: avalos-alais.s
"""

#FTRACT-DLPFC&IFG --- PLOT
import mne
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from PIL import Image
import pandas as pd
from tools.def_funtional_250 import functional_networks_250
from collections import OrderedDict
import sys
sys.path.append('ENIGMA')
from enigmatoolbox.plotting import plot_subcortical
#TODO : automat figure formating

#Path to MNE data mesh 
meshdirname = r'C:\Users\avalos-alais.s\FTRACT_LPFC\MNE-data\subjects\cvs_avg35_inMNI152'
#% Definition of LPFC 
def plot_lpfc_definition(output_path, labels_path, labels_all_33,  meshdirname = meshdirname ) : 
    
    
    figR = mne.viz.create_3d_figure(size=(400, 400), bgcolor= None, smooth_shading=None, handle=None, scene=True, show=False )
    figL = mne.viz.create_3d_figure(size=(400, 400), bgcolor= None, smooth_shading=None, handle=None, scene=True, show=False )
    Brain = mne.viz.get_brain_class()
    brainR = Brain("", hemi='rh', surf='inflated', subjects_dir=meshdirname, figure=figR,cortex = "low_contrast", title = 'LPFC Definition')
    brainL = Brain("", hemi='lh', surf='inflated', subjects_dir=meshdirname, figure=figL,cortex = "low_contrast", title = 'LPFC Definition')
    
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
    
    labels_path125 =  os.path.join(labels_path + '-Lausanne125' , 'label')  #for plots
    labels_path33 =  os.path.join(labels_path + '-Lausanne33', 'label'   )  #for plots
    # os.chdir(output_path)
    # labels_all_33 = [line.strip() for line in open('-Lausanne36.txt')]
    # labels_all_33 = [s.replace("ctx-lh-", "lh.").replace("ctx-rh-", "rh.") for s in labels_all_33] #special for Lau33
    os.chdir(labels_path33)
    for l in labels_all_33 : 
        if 'lh.' in l : 
            label_temp = mne.read_label(l +'.label')
            brainL.add_label(label_temp, color = 'k', borders = 0.00001) 
        if 'rh.' in l : 
            label_temp = mne.read_label(l +'.label')
            brainR.add_label(label_temp, color = 'k', borders = 0.00001) 
            
    os.chdir(labels_path125)
    for i in range(0, len(colors_lpfc_dict)) :  
        if 'lh.' in names_lpfc[i] : 
            label_temp = mne.read_label(names_lpfc[i] +'.label')
            brainL.add_label(label_temp, color = colors_lpfc[i]) 
        if 'rh.' in names_lpfc[i] : 
            label_temp = mne.read_label(names_lpfc[i] +'.label')
            brainR.add_label(label_temp, color = colors_lpfc[i]) 
    os.chdir(output_path)
    brainL.show_view(azimuth = -200, elevation = 90,  roll =90)  
    brainL.save_image(filename='lh_LPFC_def.png', mode='rgba')  
    brainR.show_view(azimuth = 200, elevation =-90,  roll = -90)   
    brainR.save_image(filename='rh_LPFC_def.png', mode='rgba')   
    
    brainR.close()
    brainL.close()
#%%Average connectivity plot
# fig: lateral part of both hemispheres ipsilateral avg connectivity of roi parcels. 
def plot_mean_p(mat_p, labels, labels_path, output_fig_paths, output_fig_name, title,  vmax = 0.2, vmin = 0 , cbar = 'plasma', nan_color = '#64646400', meshdirname =  meshdirname ): 
    #MAT P SHOULD BE [LEFT, RIGHT] /  vector
    figR = mne.viz.create_3d_figure(size=(400, 400), bgcolor=(255,255,255), smooth_shading=None, handle=None, scene=True, show=False)
    figL = mne.viz.create_3d_figure(size=(400, 400), bgcolor=(255,255,255), smooth_shading=None, handle=None, scene=True, show=False)
    Brain = mne.viz.get_brain_class()
    cmap = plt.get_cmap(cbar)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    values= np.linspace(vmin, vmax, 255)
    
    #LEFT
    brainL = Brain("", hemi='lh', surf='inflated', subjects_dir=meshdirname, figure=figL,cortex = nan_color)
    brainR = Brain("", hemi='rh', surf='inflated', subjects_dir=meshdirname, figure=figR,cortex = nan_color)  
    
    j_colors = np.zeros(len(mat_p), dtype=int)
    for i, x in enumerate(mat_p):
        j = np.argmin(np.abs(values - x))
        
        j_colors[i] = j
    colors =  [cmaplist[j] for j in j_colors]
    
    os.chdir(labels_path)
    c = 0 
    for r in labels : #roi labels: aff/eff we plot over the roi 
        
        if 'lh.' in r : 
             label_r = mne.read_label(r  +  '.label')
             if not np.isnan(mat_p[c]): 
                brainL.add_label(label_r, color = colors[c] ) 
        if 'rh.' in r :  
            label_r = mne.read_label(r  +  '.label')
            if  not np.isnan(mat_p[c]): 
                brainR.add_label(label_r, color = colors[c] ) 
        c = c+1;
   
    fig, axs = plt.subplots(1, 3, figsize=(10,7), gridspec_kw={'width_ratios': [1, 1, 0.2]})
   
    brainR.show_view(view="lat")
    img_R_lat = brainR.screenshot(time_viewer=True)
    img_R_lat = img_R_lat[50:-50, :, :]
    brainL.show_view(view="lat")
    img_L_lat = brainL.screenshot(time_viewer=True)
    img_L_lat = img_L_lat[50:-50, :, :]
    
    axs[0].imshow(img_L_lat)
    axs[0].axis('off')  
    axs[1].imshow(img_R_lat)
    axs[1].axis('off')
    axs[2].axis('off')
    axs[2].axis('off')
    

    fig.suptitle(title , fontsize=30)
    cm = mpl.colormaps[cbar]
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
   
    cb_ax = fig.add_axes([0.88, 0.25, 0.025, 0.4])
    cbar1 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cm), cax=cb_ax, orientation='vertical')
    cbar1.ax.tick_params(labelsize=25)
    ticks = np.linspace(vmin, vmax, 5)  
    cbar1.set_ticks(ticks)
    cbar1.ax.set_yticklabels([f"{t:.2f}" for t in ticks])

    plt.tight_layout()
    os.chdir(output_fig_paths)
    plt.savefig(output_fig_name)

    plt.show()
    brainL.close()
    brainR.close()
#%%Directed connectivity of combined resolutions (+subcortical) 
# original mat or merge in roi 
def plot_efferent (mat_all, mat_roi, labels_roi, labels_all , labels_path_roi, labels_path_all,  output_fig_paths, vmax = 0.5, vmin = 0 , cbar = 'plasma', border_print = False, flag_subcort = 'True', labels_fig_names = None, nan_color = '#64646400', meshdirname = meshdirname ):
    
    #AFF just get transpose mat 
    #fig names is title, labels roi is fig name
    
    labels_roi = [label if label.endswith('.label') else label + '.label' for label in labels_roi]
    labels_all = [label if label.endswith('.label') else label + '.label' for label in labels_all]
    
    if labels_fig_names is None : 
        labels_fig_names = [parcels.replace('rh.',  'Right ').replace('lh.',  'Left ').replace('_',' ').replace('.label', '') for parcels in labels_roi]
    
    cmap = plt.get_cmap(cbar)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    values= np.linspace(vmin, vmax, 255)
    
    #Subcortical
    if flag_subcort :  #Assumes ROI is not subcort
        subcort_path = os.path.join(output_fig_paths, 'subcort')
        os.makedirs(subcort_path, exist_ok=True)
        
        idx_l_amyg  = labels_all.index('Left-Amygdala.label')
        idx_r_amyg  = labels_all.index('Right-Amygdala.label')
        idx_l_hipp  = labels_all.index('Left-Hippocampus.label')
        idx_r_hipp  = labels_all.index('Right-Hippocampus.label')
        
        val_l_amyg = mat_all[:, idx_l_amyg]
        val_r_amyg = mat_all[:, idx_r_amyg]
        val_l_hipp = mat_all[:, idx_l_hipp]
        val_r_hipp = mat_all[:, idx_r_hipp]
        
    for s  in range(0,len(mat_all)):
        figR = mne.viz.create_3d_figure(size=(400, 400), bgcolor=(255,255,255), smooth_shading=None, handle=None, scene=True, show=False)
        figL = mne.viz.create_3d_figure(size=(400, 400), bgcolor=(255,255,255), smooth_shading=None, handle=None, scene=True, show=False)
        Brain = mne.viz.get_brain_class()
        # brain window size in pixels
        # empty string below not to go into any subdir when looking for meshes
        brainR = Brain("", hemi='rh', surf='inflated', subjects_dir=meshdirname, figure=figR,cortex = nan_color)
        brainL = Brain("", hemi='lh', surf='inflated', subjects_dir=meshdirname, figure=figL,cortex = nan_color)
        
        x_data = mat_all[s,:]
        j_colors = np.zeros_like(x_data, dtype=int)
        for i, x in enumerate(x_data):# Encontrar el índice más cercano para cada valor de x_data en values
            j = np.argmin(np.abs(values - x))
            j_colors[i] = j
        colors =  [cmaplist[j] for j in j_colors]
        
        os.chdir(labels_path_all)
        c = 0
        for r in range(0, mat_all.shape[1]): # Plot Cortical parcels Recording (starting with lh. or rh)
            if 'lh.' in labels_all[r] : 
                 label_temp = mne.read_label(labels_all[r])
                 if not np.isnan(x_data[r]): 
                    brainL.add_label(label_temp, color = colors[c] ) 
                 if border_print : brainL.add_label(label_temp, color = 'k', borders = 0.00001) 
            if 'rh.' in labels_all[r] :  
                label_temp = mne.read_label(labels_all[r])
                if not np.isnan(x_data[r]):
                    brainR.add_label(label_temp, color = colors[c] ) 
                if border_print : brainR.add_label(label_temp, color = 'k', borders = 0.00001) 

            c = c+1;
        #ROI : Overlap finer parcellation only for the ipsilateral hemisphere of the stimulation. 
        # transparent for NaN (if not enough stats for finer parcellation let see if bigger parcellation has the stats)
        x_data = mat_roi[s,:]
        j_colors = np.zeros_like(x_data, dtype=int)
        for i, x in enumerate(x_data):
            j = np.argmin(np.abs(values - x))
            j_colors[i] = j
        colors =  [cmaplist[j] for j in j_colors]
        
        os.chdir(labels_path_roi)
        c = 0
        for r in labels_roi: 
            label_temp = mne.read_label(r)
            if 'lh.' in r and 'lh.' in labels_roi[s]: 
                if np.isnan(x_data[c]): 
                    brainL.add_label(label_temp, color = nan_color, borders= False )
                else : 
                   brainL.add_label(label_temp, color = colors[c] ) 
                
                if border_print : brainL.add_label(label_temp, color = 'k', borders = 0.00001) 
          
            if 'rh.' in r and 'rh.' in labels_roi[s] : 
                if  np.isnan(x_data[c]): 
                    brainR.add_label(label_temp, color = nan_color, borders= False )
                else : 
                    brainR.add_label(label_temp, color = colors[c] ) 
                    
                if border_print : brainR.add_label(label_temp, color = 'k', borders = 0.00001) 
               
            c = c+1;
        # Subcortical
        if flag_subcort :
             sub_cort =  np.empty(14) ; sub_cort[:] =  np.nan
             sub_cort[1] = val_l_amyg[s] # s of the stimulated parcel 
             sub_cort[3] = val_l_hipp[s]
             sub_cort[8] = val_r_amyg[s]
             sub_cort[10] = val_r_hipp[s] 
             pd.DataFrame(sub_cort, columns=['Values'])
             
             os.chdir(subcort_path)
             
             fig_sc = plot_subcortical(array_name=sub_cort, size=(1000, 500),
             
              cmap=cbar, color_bar=False, ventricles=False, color_range=(vmin, vmax) , screenshot = True ,  filename= labels_roi[s] + '_subcort.png' )

        #FIGURES 
        
        fig_name = labels_fig_names[s]
        
        brainR.show_view(view="lat")
        img_R_lat = brainR.screenshot(time_viewer=True)
        img_R_lat = img_R_lat[40:-40, 3:-3, :]
        brainR.show_view(view="med")
        img_R_med = brainR.screenshot(time_viewer=True)
        img_R_med = img_R_med[40:-40, 3:-3, :]
        
        brainL.show_view(view="lat")
        img_L_lat = brainL.screenshot(time_viewer=True)
        img_L_lat = img_L_lat[40:-40, 3:-3, :]
        brainL.show_view(view="med")
        img_L_med = brainL.screenshot(time_viewer=True)
        img_L_med = img_L_med[40:-40, 3:-3, :]
        
        fig, axs = plt.subplots(3, 3, figsize=(20, 22), gridspec_kw={'height_ratios': [1, 0.001, 1], 'width_ratios': [1, 1, 0.1]})

        axs[0, 0].imshow(img_L_lat)
        axs[0, 0].axis('off')  
        axs[0, 1].imshow(img_R_lat)
        axs[0, 1].axis('off')
        axs[2, 0].imshow(img_L_med)
        axs[2, 0].axis('off')
        axs[2, 1].imshow(img_R_med)
        axs[2, 1].axis('off')
        
        for ax in axs[:, 2]:
            ax.axis('off')
        
        axs[1, 0].remove()  
        axs[1, 1].remove()  
        axs[1, 2].remove()  
        
        if flag_subcort:
            subcort_ax = fig.add_axes([0.08, 0.45, 0.8, 0.4]) 
        
            img_subcortical = Image.open(fig_sc)        
            subcort_ax.imshow(np.array(img_subcortical))
            subcort_ax.axis('off')
        
        fig.suptitle(fig_name, fontsize=60)
        
        cm = mpl.colormaps[cbar]
        norm = plt.Normalize(vmin=vmin, vmax=vmax)

        cb_ax = fig.add_axes([0.93, 0.4, 0.02, 0.4]) 
        cbar1 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cm), cax=cb_ax, orientation='vertical')
        cbar1.ax.tick_params(labelsize=30)
        
        plt.tight_layout(rect=[0, 0.20, 1, 1])
        os.chdir(output_fig_paths)
        fig_name = fig_name.replace(' ','_')
        plt.savefig(fig_name + '.svg')
         
        plt.show()
        brainL.close()
        brainR.close()
    

#%%Functional Networks 
#TODO: Generalize to plot any def : take out function of def, names and index as params. 
def plot_functional_net (path_defs, list_labels_roi, labels_path125, labels_nets, colors, functional_fig, meshdirname = meshdirname  ) :   
    #plot the definition of the special segmentation 
    [names_labels_regions, index_net, index_roi] = functional_networks_250(path_defs)
    labels_roi = [list_labels_roi[idx] for idx in index_roi]
    
    Brain = mne.viz.get_brain_class()
    subjects_dir = mne.datasets.sample.data_path() / "subjects"
   
    mne.datasets.fetch_hcp_mmp_parcellation(subjects_dir=subjects_dir, verbose=True)
   
    mne.datasets.fetch_aparc_sub_parcellation(subjects_dir=subjects_dir, verbose=True)
   
    mne.utils.set_config("SUBJECTS_DIR", subjects_dir, set_env=True)
   
    brainL = Brain("cvs_avg35_inMNI152-Lausanne125","lh","inflated",subjects_dir=subjects_dir,cortex='#64646400',background="white",size=(800, 600),title = 'Definition' )
    brainR = Brain("cvs_avg35_inMNI152-Lausanne125","rh", "inflated",subjects_dir=subjects_dir,cortex='#64646400',background="white",size=(800, 600),title = 'Definition' )
   
    #roi :
    os.chdir(labels_path125)
    for p in labels_roi : 
        label_p =  mne.read_label(p+'.label')
        if 'lh' in p : 
            brainL.add_label(label_p, borders = 0.00001, color = 'k' ) 
        if 'rh' in p : 
            brainR.add_label(label_p, borders = 0.00001, color = 'k' ) 
           
    #networks
    os.chdir(labels_nets)
    c = 0
    for l in names_labels_regions : 
        label_temp = mne.read_label(l+ '.label')
        if l.startswith('lh.') : 
            brainL.add_label(label_temp, color = colors[c]) 
            c=c+1
        if l.startswith('rh.') : 
            brainR.add_label(label_temp, color = colors[c]) 
            c=c+1
   
    patches = []
    nets = names_labels_regions[:len(names_labels_regions)//2 ]
    nets = [n.replace('lh.', '') for n in nets]
    for color, text in zip(colors, nets):
        patch = mpatches.Patch(color=color, label=text)
        patches.append(patch)
   
    
    brainL.show_view(view="lat")
    img_L_lat = brainL.screenshot(time_viewer=True)
    img_L_lat = img_L_lat[30:-30, 30:-30, :]
    
    brainL.show_view(view="med")
    img_L_med = brainL.screenshot(time_viewer=True)
    img_L_med = img_L_med[30:-30, 30:-30, :]
    
    brainR.show_view(view="lat")
    img_R_lat = brainR.screenshot(time_viewer=True)
    img_R_lat = img_R_lat[30:-30, 30:-30, :]
    
    brainR.show_view(view="med")
    img_R_med = brainR.screenshot(time_viewer=True)
    img_R_med = img_R_med[30:-30, 30:-30, :]

  
    fig, axs = plt.subplots(2, 2, figsize=(40, 20), gridspec_kw={'width_ratios': [1, 1] , 'height_ratios': [1,1]})

    axs[0,0].imshow(img_L_lat)
    axs[0,0].axis('off')  
    axs[1,0].imshow(img_L_med)
    axs[1,0].axis('off')
    axs[0,1].imshow(img_R_lat)
    axs[0,1].axis('off')
    axs[1,1].imshow(img_R_med)
    axs[1,1].axis('off')
 
    fig.suptitle('Functional Networks', fontsize=80)
    fig.legend(handles=patches, loc='upper right', fontsize=50) 
    plt.tight_layout()
    plt.subplots_adjust(wspace=-0.55, hspace=-0.15)  

    os.chdir(functional_fig)
    plt.savefig('Functional_Networks.svg')

    plt.show() 
    
def plot_efferent_fn(path_table, table, names_labels_regions,  roi , path_labels_regions, path_labels_roi , parcellation_file , path_output,  cbar_max ,roi_color = '#72B2F4' ,  colormap  = 'plasma', subject = 'cvs_avg35_inMNI152 -Lausanne250', meshdirname = meshdirname ):
    #TODO: generalize to all specific segmentaions
    #TODO: genrealize for afferent
    #plot directed efferent connectivity from roi to specific segmetation 
    roi = [label if label.endswith('.label') else label + '.label' for label in roi]
    names_labels_regions = [label if label.endswith('.label') else label + '.label' for label in names_labels_regions]
    
    Brain = mne.viz.get_brain_class()
    subjects_dir = mne.datasets.sample.data_path() / "subjects"
    
    mne.datasets.fetch_hcp_mmp_parcellation(subjects_dir=subjects_dir, verbose=True)
    
    mne.datasets.fetch_aparc_sub_parcellation(subjects_dir=subjects_dir, verbose=True)
    
    mne.utils.set_config("SUBJECTS_DIR", subjects_dir, set_env=True)
    X = table  
    #Normalize to the max in all data
    min_value = np.nanmin(X)
    max_value = np.nanmax(X)
    #%% Colors
    cmap = plt.get_cmap(colormap)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    if cbar_max is not None:
        values= np.linspace(0,cbar_max, 255)
    else : values= np.linspace(min_value,max_value, 255) 
    for l in range(0,len(names_labels_regions)) :
        x_data = X[:,l]  # all dlpfc - one region recording
     
        j_colors = np.zeros_like(x_data, dtype=int)
    
        for i, x in enumerate(x_data):
            j = np.argmin(np.abs(values - x))
            j_colors[i] = j

        if 'lh' in names_labels_regions[l]:
            brain = Brain(
                subject,
                "lh",
                "inflated",
                subjects_dir=subjects_dir,
                cortex="low_contrast",
                background="white",
                size=(800, 600),
                title = 'Efferent connectivity')
        else :
            brain = Brain(
                subject,
                "rh",
                "inflated",
                subjects_dir=subjects_dir,
                cortex="low_contrast",
                background="white",
                size=(800, 600),
                title = 'Efferent connectivity' )
                           
       #plot result of roi stims in region over roi
        os.chdir(path_labels_roi) #labels125 - dlpfc

        for p in range(0, len(roi)) :
            
            if 'lh.' in  names_labels_regions[l] and 'lh.' in roi[p]:
                label_x = mne.read_label(roi[p] )
                if not np.isnan(x_data[p]):
                    brain.add_label(label_x, color =cmaplist[j_colors[p]],borders= False )
                    
            elif 'rh.' in names_labels_regions[l] and 'rh.' in roi[p] :
                label_x = mne.read_label(roi[p])
                if not  np.isnan(x_data[p]):
                    brain.add_label(label_x, color = cmaplist[j_colors[p]], borders=False )

        #plot region of interest:     
        os.chdir(path_labels_regions)
        rec_region = mne.read_label(names_labels_regions[l] )
        brain.add_label(rec_region, color = roi_color[l] , borders = False) #dlpfc stimulated
        fig_name =  names_labels_regions[l].replace('.label','')
        fig_title = fig_name.replace('lh.', 'Left ').replace('rh.', 'Right ').replace('_Lau125_dlpfc', '')
        fig, axs = plt.subplots(1, 2, figsize=(5, 5), gridspec_kw={'width_ratios': [1,  0.2]}) # saco el lado medial
       
        brain.show_view(view="lat")
        img_lat = brain.screenshot(time_viewer=True)
        img_lat = img_lat[80:-80, 80:-80, :]
    
        axs[0].imshow(img_lat)
        axs[0].axis('off') 
        axs[1].set_position([0.1, 0.2, 1, 0.5]) 
        axs[1].axis('off')
        
        fig.suptitle(fig_title , fontsize=30)
        
        norm = plt.Normalize(vmin=0, vmax=cbar_max)
        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs[1], orientation='vertical', fraction=0.5, pad=0.5)
        cbar.ax.tick_params(labelsize=10)

        plt.tight_layout()
            
        os.chdir(path_output)
        plt.savefig(fig_name + '.svg',format='svg', dpi=1000)
    
        plt.show()
        brain.close()
#%%Matrix plot
def matrix_plot(mat, parcel_stim_names, parcel_rec_names, title, fig_name, output_path, colorm = 'plasma', vmax = 0.5) : 
    #Not last version we use in the paper ( keys missing) 
    #plot horizontal  labels_list,
    plt.figure(figsize=(40, 10))  

    im = plt.imshow(mat, aspect=1, cmap=colorm) 
    
    cbar = plt.colorbar(im, fraction=0.03)
    cbar.ax.tick_params(labelsize=30)
    im.set_clim(vmin=0, vmax=vmax)

    plt.yticks(ticks=np.arange(0, len(parcel_stim_names), 5), labels=np.arange(0, len(parcel_stim_names), 5),size=30)
    plt.xticks( ticks=np.arange(0, len(parcel_rec_names), 5), labels=np.arange(0, len(parcel_rec_names), 5), rotation=90, size=30, ha='right')
    
    plt.title(title, size=40)
    
    os.chdir(output_path)
    plt.savefig(fig_name, bbox_inches='tight', dpi=1000)
    plt.show()
#%%tool: convert matrix to vector
def convert_to_vector(matrix):
    rows, columns = matrix.shape
    # Case: If the matrix has 2 rows and many columns
    if rows == 2:
        # Concatenate the two rows
        vector = np.concatenate((matrix[0, :], matrix[1, :]))
    # Case: If the matrix has 2 columns and many rows
    elif columns == 2:
        # Concatenate the two columns
        vector = np.concatenate((matrix[:, 0], matrix[:, 1]))
    else:
        raise ValueError("The matrix must have exactly 2 rows or 2 columns.")

    return vector
   


