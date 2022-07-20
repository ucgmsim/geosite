Main changes in Version 1.2 - 04-02-2022:

- The filter to handle negative values of qc and fs from the CPT data is now also valid for zero and non-existent values of qc and fs (gaps in CPTu data). 

- A tcl file containing an OpenSees site-response model is automatically generated. This tcl file is ready to be used within the workflow for conducting multiple site-response analyses implemented in Mahuika.

- An Excel file called [Site]_CPTu_trial_models.xlsx is automatically generated. This file contains the average values of the parameters estimated using CPu-based correlations, for each layer defined by the user.

- Lines to identify layers boundaries were added in the Ic plot.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Main changes in Version 1.1 - 19-11-2021:

- The only path that has to be defined now is a root directory (called rootDir). The inputPath and outputPath are defined relative to rootDir.

- The CPTu, SCPT, and SW objects are now handled by the lists CPTu_Data, SCPT_Data, and SW_Data. 
  This was done with the purpose of handling more easily multiple tests per each type (implementation still in progress).
  New class objects created (additional tests) will be appended to those lists, instead of having to define a different and arbitrary name per each object.

- A filter was added to handle a negative value of the water table, GWL. If a negative value is introduced, it is converted to a positive one. 

- A filter was added to handle negative values of qc and fs from the CPT data. 
  If there are negative values of qc or fs: 
  	(1) Those values are replaced by 0.01.
	(2) The associated gamma values are replaced by a standard value of 18 (kN/m3) to allow all the subsequent computations (vertical total stress, vertical effective stress, etc.).
	(3) The associated values of Ic are not shown (they are replaced by float('NaN')).
	(4) The associated parameters estimates of Vs, relative density, and friction angle are not shown (they are replaced by float('NaN')).

- The limit of the x-axis for the plotting of Vs is now automatically defined as: maximum value of Vs in the user-defined model + 30 (m/s).

- The limit of the y-axis for the plotting of the parameter estimates (Ic, Vs, gamma, relative density, and friction angle) is now automatically defined as: maximum value of depth in the user-defined model + 5 (m).



