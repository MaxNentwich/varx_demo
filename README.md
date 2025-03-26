# varx_demo

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
    + fig_1_filter_example.m plots examples of individual feature and data time courses as well as VARX input-response (B) filters.
    + fig2_analyze_model_output.m plots results of the neural mass model simulation analysis.
    + fig3_var_vs_varx.m compares models without external input.
    + fig4_dme_vs_rest.m compares resting state and movie data.
        + movie_rest_repeated_measures.m performs corresponding statistics with a mixed-effect model.
        + connectivity_plot_movie_vs_rest.py plots connections for an example patient on a brain.
    + fig5_input_H_vs_TRF.m compares input filters from the VARX and mTRF methods.
    + fig6_innovation_rest_vs_movie_hfa.m examines the effects of modeling inputs on noise in the model.
    + fig7_asymetries.m compares the degreee of incoming and outgoing connections to cortical hierarchy defined by myelination.
        + plot_varx_dir_myelin_dk.py creates the corresponding brain surface plot     
     

## Python Code

Brain plots are created with python code. Plots in Figure 7 use code from: https://github.com/rdgao/field-echos/tree/master. The environment.yml file provides a list of dependencies. 

hierarchy_dk_atlas.py downloads myelin maps and transforms them in the right space and needs to be run before creating the figure with plot_varx_dir_myelin_dk.py . This requires the Connectome Workbench to be installed: https://www.humanconnectome.org/software/connectome-workbench
