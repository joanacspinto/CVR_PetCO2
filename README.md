# CVR_PetCO2

# BOLD PRE-PROCESSING (in FSL FEAT)

1st - Perform one Feat (GUI), on one illustrative subject, where we select all the options we want 
2nd - Select the design.fsf file created and alter the subject/session data using the sed function to the appropriate subject/session creating a new design.fsf. 
3rd - Call the Feat function using the new design.fsf (feat design.fsf)

In FEAT:

INPUT: BOLD.nii.gz

DATA TAB: SET input, SET output directory - e.g. BOLD.feat, SET high pass filter cutoff (e.g. 100s)

PRE-STATS TAB: SET motion correction, SET BET brain extraction, SET spatial smoothing FWHM (e.g. 5 mm), SET highpass filtering

OUTPUT OF INTEREST: filtered_func_data.nii.gz


# PETCO2 PROCESSING (in MATLAB)

INPUT: raw PETCO2 time course

OUTPUT: PETCO2_EV.txt - output time-course to be used as regressor in GLM (normalize 0 to 1, with same resolution of BOLD) and dPETCO2.txt - output value of change in PETCO2 (task - baseline, in mmHg)
