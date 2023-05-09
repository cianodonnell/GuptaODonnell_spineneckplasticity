# GuptaODonnell_spineneckplasticity
MATLAB code for computational model of dendritic spine neck and synaptic plasticity and AMPA receptor dynamics, written by Rahul Gupta.

From research paper: [Dendritic spine neck plasticity controls synaptic expression of long-term potentiation](https://www.biorxiv.org/content/10.1101/2023.01.27.525952v1)

Instruction to run the simulation:

*Step 1:* Open the script 'PLAY.m' in the MATLAB. 

*Step 2:* Set the conditions, such as  
- the location of the stimulated spine StimSpineLoc, in micron meters (um), away from the  dendritic branch origin. The branch is 100um long and the spines are located at a regular distance of 1um.
- whether the spine neck of the stimulated changes its morphology (neck morphological plasticity) upon stimulation.
  - If no neck morphological plasticity is chosen, set the variable 				NeckMorphPlast=0.
  - if the the neck undergo widening and shortening, set the variable 				NeckMorphPlast=1.
  - if the the neck undergo narrowing and elongation, set the variable 				NeckMorphPlast=2.
- If the neck morphological plasticity is considered to occur, set the time-scale tauNM of the neck morphological change.
- whether the septin7 ring membranous ring barrier undergoes a change (neck septin plasticity) upon stimulation. The change considered here only involves weakening or fading away of the barrier
  - If no neck septin plasticity is chosen, set the variable 					NeckSepPlast=0.
  - if neck septin plasticity occurs, set the variable 						NeckSepPlast=1.
- If the neck septin plasticity is considered to occur, set the time-scale tauNS of the 	neck septin change.

*Step 3:* Simply run the script. 

The outputs will appear as 4 figures:
- Figure 1: Temporal profile of the Perecent change in stimulated spine head surface area, taken from Patterson et al. 2010.
- Figure 2: Temporal profile of CaMKII-dependent AMPAR phosphorylation rate constant (k_phos) in stimulated spine head.
- Figure 3: Temporal profile of the concentration of active calcineurin [CaN], in micronMolar, in stimulated spine head and the dendritic shaft immediately underneath the stimulated spine.
- Figure 4: Temporal profile of the percent change in the total surface AMPAR fluorescence in the stimulated spine head. The experimental data from the Patterson et al. and the simulated data from the model are overlaid.

Presently, the PLAY.m is set to single spine stimulation at the spine location 50um, with no neck morphological or septin plasticity.

Note: Values of tauNM and tauNS don't intervene the simulation if the NeckMorphPlast and NeckSepPlast are set to zeros.
