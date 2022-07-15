# Sea_Clearly
Repo to run parcels simulations, and to postprocess the results, visualizing origins and destinations of marine pollution such as microplastics

Step 1):

If no Lagrangian simulation data is available yet, perform the following steps, otherwise go to step 2)

To run the Lagrangian simulations, run for example:
python sea_clearly/main_advection.py -date_start_release 2018-01-01 -date_end_release 2018-01-31 -n_days 300

Make sure that the DIR_UV, PATTERN_U and PATTERN_V variables are set correctly in the settings file, which point to the netcdf files containing the ocean current information. Also make sure DIR_INPUT and DIR_OUTPUT are set correctly.


Step 2):

To visualize the results, you can open the notebook sea_clearly/postprocess.ipynb

Set the patterns to correspond to the output files that were created in step 1. A forward and a backward analysis+animation are available
