# CVR AMPLITUDE AND DELAYS USING PETCO2 MODEL FITTING
Functions: CVR_compute.m and PETCO2_proc.m
(CVR_compute.m calls PETCO2_proc.m)

1- Pre-process BOLD data using FSL's FEAT -> filtered_func_data.nii.gz

2- Call "CVR_compute(filtered_func_data.nii.gz, raw physio dataset, shift start, shift end)"

CVR_compute creates 3 images: (1) CVR amplitude without delay optimization (delay = 0), (2) CVR amplitude with delay optimization, (3) CVR delay

CVR in %/mmHg and CVD in seconds

![alt tag](http://url/to/img.png)

