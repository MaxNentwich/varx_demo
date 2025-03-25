# varx_demo

Code for analysis in: [https://doi.org/10.1101/2024.08.05.606665](https://doi.org/10.7554/eLife.104996.1
)

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
