# CVR AMPLITUDE AND DELAYS USING PETCO2 MODEL FITTING
Note: CVR in %/mmHg and CVD in seconds

### If the user has the raw physio PETCO2 time course - function: CVR_compute.m (CVR_compute.m calls PETCO2_proc.m)

This is useful for time shifting, in order to perform this step without losing any BOLD MR datapoints.

*PIPELINE*

1- Pre-process BOLD data using FSL's FEAT -> to obtain filtered_func_data.nii.gz

2- Call "CVR_compute(filtered_func_data.nii.gz, raw physio dataset, shift start, shift end)"

This function call PETCO2_proc to process the PETCO2 raw data (find peaks, shift, etc.)

CVR_compute creates 3 images: 

* CVR amplitude without delay optimization (delay = 0)

* CVR amplitude with delay optimization

* CVR delay


### If the user only has the trimmed and processed PETCO2 time course - function: CVR_compute_processed.m

In this case, the time shifting will remove some BOLD MR datapoints, in order to have the same dimension as the PETCO2 regressor.

*PIPELINE*

1- Pre-process BOLD data using FSL's FEAT -> to obtain filtered_func_data.nii.gz

2- Pre-process PETCO2 data -> to obtain PETCO2_processed and delta_PETCO2_processed (the difference between PETCO2 in hypercapnia and baseline)

Call "CVR_compute_processed(filtered_func_data.nii.gz, PETCO2_processed, shift start, shift end, delta_PETCO2_processed)"

CVR_compute creates 3 images:

* CVR amplitude without delay optimization (delay = 0)

* CVR amplitude with delay optimization

* CVR delay
