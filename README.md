# varx_demo

<a href="https://doi.org/10.5281/zenodo.15127333"><img src="https://zenodo.org/badge/832732261.svg" alt="DOI"></a>

Code for analysis in: https://doi.org/10.7554/eLife.104996.1

VARX toolbox: https://github.com/lcparra/varx/tree/main

Figure data: https://doi.org/10.17605/osf.io/vc25t


## Analysis steps

+ Preprocessing
    + /data_prep/prep_rest.m preprocesses raw data.  
    + /data_prep/convert_data.m reorganized data into a simpler format. 
    + These scripts run on raw data, which is not provided, but are shared for reference.
+ Computing VARX models
    + main.m computes all versions of VARX models used in the analysis.
+ Figures
    + fig_1_filter_example.m plots examples of individual feature and data time courses as well as VARX input-response (B) filters - **Figure 1**
    + Neural mass model: /src/neurolib_brainmodel contains code for the simulation of the neural mass model - **Figures 2, S1 & S2**     
    + fig3_var_vs_varx.m compares models without external input - **Figures 3A-E & S5A-E**
        + extended_models.m evaluates the effect of individual features, as well as the effect of correlated features and intrinsic connecitons on input filters - **Figures 3F, S5F & S11** 
    + fig4_dme_vs_rest.m compares resting state and movie data - **Figures 4A-B & S7A-B**
        + movie_rest_repeated_measures.m performs corresponding statistics with a mixed-effect model - **Figures 4C-D & S7C-D**
        + connectivity_plot_movie_vs_rest.py plots connections for an example patient on a brain - **Figure 4E-G**
        + eyes_closed_rest_simple.m compares movie data with eyes-closed rest - **Figure S8**
    + fig5_input_H_vs_TRF.m compares input filters from the VARX and mTRF methods - **Figures 5, S9 & S10**
    + fig6_innovation_rest_vs_movie_hfa.m examines the effects of modeling inputs on noise in the model - **Figure 6**
    + fig7_asymetries.m compares the degreee of incoming and outgoing connections to cortical hierarchy defined by myelination - **Figures 7A&C & S14A&C**
        + plot_varx_dir_myelin_dk.py creates the corresponding brain surface plot - **Figures 7B & S14B**
    + data_summary.m summarizes demographics and the length of recordings - **Table S1**
    + figSX_compare_na6_na3.m compares examples of models with a different number of parameters - **Figure S3**
    + eye_params_movies_rest.m compares eye movements in different movies - **Figure S6**
    + fig_S6_innovation_rest_vs_movie_lfp.m examines the innovaiton process in the LFP data - **Figure S12**
    + figS8_simulate_gain_adaptation.m simulates a gain adaptation model - **Figure S13**
     

## Python Code

Brain plots are created with python code. Plots in Figure 7 use code from: https://github.com/rdgao/field-echos/tree/master. The environment.yml file provides a list of dependencies. 

hierarchy_dk_atlas.py downloads myelin maps and transforms them in the right space and needs to be run before creating the figure with plot_varx_dir_myelin_dk.py . This requires the Connectome Workbench to be installed: https://www.humanconnectome.org/software/connectome-workbench
