# LPFC Connectivity Mapping Project

## Overview
This repository contains scripts and data for generating matrices and brain maps that illustrate "High Resolution and Electrophysiological Mapping of Effective Connectivity of Lateral Prefrontal Cortex" by Avalos-Alais & Jedynak et al. The analysis is based on sEEG stimulation data from the FTRACT database.

## Data Description
While raw data is not shared, we provide processed data matrices containing intracerebral evoked potentials (iEPs) parameters:
- Probability with corresponding confidence intervals
- Peak delays
The number of recordings for each connection is also provided noted by 'N'.

The matrices are organized with stimulated parcels in rows and recorded parcels in columns, using Lausanne-2008 parcellations at different resolutions.

#### Time Windows of analysis of connectivity
- General data: 0-100ms
- Direct connectivity: 0-50ms
- Indirect connectivity: 100-400ms

## Repository Structure

### Main scripts
- `compute_atlas_mats`: Computation of connectivity matriceswith the particular considerations of this project. Taking FTRACT raw connectivity matrices in Lausanne parcellations it coordinates de selection of pertinent rows and columns and the masking to avoid outliers. 
- `matrices`: Handles data processing and matrix generation - called by compute_atlas_mats
- `plotting`: Folder of visualization functions (matrices, bar_plots and brain maps)
- `figures`: Generates paper figures reading processed data (data shared in `Results`)
- 'Definitions' : Contains main definitions of ROIs, eg. which Lausanne2008-125 parcels compose the 'LPFC', or functional networks, etc. Basic functionalities of indexing.

### Modules and Data needed to generate results
- ENIGMA: module used for subcortical plotting
- MNE_data: contains mesh and '.label' object to plot brains
- matplotlib_curly_brace: module used to plot curly braces with region data in matrices

The complete environnement used is shared in a .yml and main requirements are also specified in a text file, automatically included with setup_env.sh  

### Results Directory
All results showed in the paper are provided for connectivity computed for Zth5. We also share the main matrix of efferent and afferent connectivity (Lausanne2008125 to 33 parcellation) computed for general time window 0-100ms for thresholds 3, 4, 6, and 7. 

Contains processed matrices with the following naming convention:
`parameter_stimulatedParcelResolution_recordedParcelResolution_timeWindow`

Example: `p_125to33_0_100ms` represents:
- Efferent probability of connectivity 
    - From LPFC parcels (Lausanne2008-125)
    - To brain regions (Lausanne2008-33)
- Time window: 0-100ms

For analysis including custom merging of regions other terms appear (eg. 'all' to refer to all ipsilateral hemisphere merged as one parcel; and 'nets' for brain parcellated into functional networks).
Some cases also include 'L' and 'R' referring to hemisphere sides. 

#### Results Organization

1. **AVG**
   - Average efferent and afferent connectivity between Lausanne2008-125 parcels and the ipsilateral brain averaged
   - Contains: CI, Index, N, probability matrices
   - [Figure 2B]

2. **AVG_LPFC_AVG_all**
   - Average connectivity between LPFC merged as one parcel and the ipsilateral brain merged as one parcel
   - Contains: CI, Index, N, probability matrices

3. **Functional_Networks**
   - LPFC Lausanne2008-125 parcels to ipsilateral functional networks
   - Contains: CI, Index, N, probability matrices, Labels_nets
   - Segmented ROI analysis : anterior/posterior/inferior/superior DLPFC and IFG 
	- LPFC segments to ipsilateral functional networks
 	- Contains: CI, Index, N, probability matrices
   - [Figure 4]

4. **Lausanne2008-33_125 : Afferent LPFC connectivity**
   - Connectivity from the whole brain in Lausanne2008-33 parcellation to LPFC parcels in Lausanne2008-125 parcellation.
   - Contains: CI, Index, N, probability matrices
   - [Figure 2C & 2E second row]

5. **Lausanne2008-125_33 : Efferent LPFC connectivity**
   - Connectivity from LPFC parcels in Lausanne2008-125 parcellation to the whole brain in Lausanne2008-33 parcellation.
   - Contains: CI, Index, N, probability matrices
   - [Figure 2C & 2E second row]

6. **Mean_ROI_Eff**
   - Directedness of connectivity from DLPFC and IFG as individual merged parcels to the rest of the brain in Lausanne2008-33 parcellation. 
   - Contains : DLPFC and IFG label objets for plotting
		Direct connectivity (0-50ms)
		- Contains: CI, Index, N, probability matrices
   		Indirect connectivity (100-400ms)
   		- Contains: CI, Index, N, probability matrices
   - [Figure 3]

7. **N_implanted_contacts**
   - Number of implanted contacts on the LPFC or recording LPFC stimulations. 
   - [Figure 2A]
   
8. **Resolutions**
   - Resolution comparison data, three LPFC stimulations for brain parcelled in Lausanne2008-33/125/500
   - [Figure 1DEF]

## Technical Details
- CI : binomial confidence intervals for computed probability of connectivity, alpha = 0.05
- Index : technical information for computation / plotting
- N : number of recordings used for the computation of the connectivity 
- Labels_nets : Object '.label' use to plot. 


### Matrix Generation
- ROI definitions use Lausanne2008-125 resolution (Definitions)
- LPFC consists of DLPFC + IFG for both hemispheres
- Filtering retains cortical parcels plus amygdala and hippocampus
- Combined resolution analyses have roi-roi connectivity and roi-all connectivity. Supperposition in brain plots with transparent parcel if roi-roi is nan.
- Statistical masking applied based on recording numbers and EP for delays. 

### Parcel Merge 
The merge of Lausanne2008 fine resolution parcels to create personalized regions is done previous to the computation of probability and not included in these scripts. 

### Matrix Correction Process
Using `matrices.atlas_mat`:
1. Takes original matrices (xx, xy)
2. Applies dimensional indices
4. Generates statistical masks
5. Saves corrected matrices as text files

### Figure Generation
The `figures` script generates all paper figures with corresponding supplementary materials. Outputs include:
#### Figure 1: Methodology Overview
- **D**: Lausanne2008-33 parcellation resolution
  - Demonstrates probability of connectivity for LPFC parcel 'Left_rostralmiddlefrontal'
- **E**: Lausanne2008-125 parcellation resolution
  - Demonstrates probability of connectivity for LPFC parcel 'Left_rostralmiddlefrontal_1'
- **F**: Lausanne2008-500 parcellation resolution
  - Demonstrates probability of connectivity for LPFC parcel 'Left_rostralmiddlefrontal_22'
- **G**: LPFC Definition
  - Shows left hemisphere LPFC definition
  - Includes right hemisphere LPFC definition (part of Fig4A and supplementary Fig3A)

#### Figure 2: Probabilistic Effective Connectivity from and to the LPFC
- **A**: Implanted Contacts Analysis
  - Number of contacts used for stimulation/recording of LPFC
  - Combined parcellation: LPFC in Lausanne2008-125, rest of brain in Lausanne2008-33
- **B**: Average Connectivity
  - Efferent and afferent connectivity of LPFC Lausanne2008-125 parcels
  - Towards ipsilateral brain hemisphere (merged as one)
- **C**: Probability of Effective Connectivity (Brain Plots)
  - Individual LPFC Lausanne2008-125 parcels to Lausanne2008-33 brain
  - Includes Lausanne2008-125 to Lausanne2008-125 connectivity over LPFC
  - **I**: Efferent connectivity (right hemisphere stimulation)
  - **II**: Afferent connectivity (right hemisphere stimulation)
- **E**: Probability of Effective Connectivity (Matrix Plots)
  - Similar to Figure C, but in matrix format
  - **I**: Efferent connectivity
  - **II**: Afferent connectivity
- **G**: Symmetry Analysis
  - Scatter plot for afferent and efferent connectivity symmetry

#### Figure 3: Direct and Indirect Effective Connectivity of the LPFC
Contains brain plots of right hemisphere stimulation (corresponding to supplementary Fig. 2)

##### DLPFC Section
- **A**: Number of stimulation recordings on DLPFC as one parcel
- **B**: Direct Connectivity
  - Probability of DLPFC connectivity to Lausanne2008-33 brain
  - Peak mean delays for iEPs in 0-50ms time window
- **C**: Indirect Connectivity
  - Probability of DLPFC connectivity to Lausanne2008-33 brain
  - Peak mean delays for iEPs in 100-400ms time window

##### IFG Section
- **A**: Number of stimulation recordings on IFG as one parcel
- **B**: Direct Connectivity
  - Probability of IFG connectivity to Lausanne2008-33 brain
  - Peak mean delays for iEPs in 0-50ms time window
- **C**: Indirect Connectivity
  - Probability of IFG connectivity to Lausanne2008-33 brain
  - Peak mean delays for iEPs in 100-400ms time window

#### Figure 4: Probabilistic Effective Connectivity from LPFC to Functional Networks
- **A**: Definitions
  - LPFC sub-segmentation:
    - Anterior/Posterior DLPFC
    - IFG
    - Superior/Inferior DLPFC and IFG
  - Brain segmentation by functional networks (LPFC excluded)
  - LPFC parcellation using Lausanne2008-125
  - Includes right hemisphere segments and network definitions
- **B**: LPFC Segments to Functional Networks
  - Bar plots showing:
    - Probability of connectivity
    - Number of recordings
  - Right hemisphere segments stimulation
- **C**: LPFC Lausanne2008-125 Parcels to Brain Networks
  - Right hemisphere representation
  - **I**: Brain plots
  - **II**: Bar plots

#### Supplementary: Figure S4. Probability of connectivity vs mean peak delay. 
Analysis of connectivity of Figure3.
For DLPFC and IFG direct and indirect connections, 
scatter plot of probability of effective connectivity of roi to the rest of the brain against the mean peak delay of the connections.

## Usage Notes
- Raw data processing requires the full FTRACT database (not included)
- Figure generation needs : 
	- MNE python freesurfer data (provided in folder 'MNE-data'
	- MNE python, we used version 1.8.0
	- Matrices of data provided in 'Results'
	- ENIGMA toolbox for subcortical plotting, provided
	- Matplotlib braces, toolbox provided 
