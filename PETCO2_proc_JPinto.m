% Joana Pinto, December 2020, Oxford; Adapted from Daniel Bulte 2018 script
% This function processes the PetCO2 timeseries to be used as regressor
% (including adding a time shift)
% INPUT: physfile_in - PETCO2 raw data; output_folder - folder where to save regressors; physio_shift - physiological time shift (time shift across brain); mechan_del - mechanical time shift (tubing, pump, etc);
% OUTPUT: out_co2 - PETCO2 processed and normalised, ready to be used as regressor; 

function [out_co2] = PETCO2_proc_v2(physfile_in, output_folder, physio_shift, mechan_shift)
addpath('/usr/local/fsl/etc/matlab');

% UNCOMMENT LINE FOR DEBUGGING
physfile_in = '/Users/joana/Documents/Work_Oxford/Data_extra/OHBA/006/2018_119_006.txt'; % Path to the physiological data file
output_folder = '/Users/joana/Documents/Work_Oxford/Data_extra/OHBA/output_test'; 
physio_shift = 0;
mechan_shift = -15;

%% SET PARAMETERS (always confirm these)
% Set here manually, if not possible to retrieve from scanner triggers 
vols = []; 
tr = []; 
% Thresholds and delays
threshold_trig = 3; % threshold to detect triggers (if necessary)
shift = mechan_shift + physio_shift; % total delay (time_shift)
% CO2 sampling 
samp_rate = 100; % Powerlab sampling rate (Hz)
air_pres = 1020; % Enter barometric pressure in mbars:

% FIND PEAKS options
MinPeakDistVal = 500;
MinPeakPromVal = 0.5;
MinThresholdVal = 2;

%% LOAD DATA

physfile = load(physfile_in);
timings = physfile(:,1);
PETCO2 = physfile(:,2); % Second column is the PETCO2 values
PETO2 = physfile(:,3); % Third column is the PETO2 values
trig = physfile(:,4); % MR triggers

% CHECK BEGINNING OF TRIGGERS (REMOVE INITIAL OUTLIERS)
begin_tri = 2; 
trim_time = timings(trig>threshold_trig);
for i = 2:size(trim_time)-1
    if trim_time(i+1)-trim_time(i)>40
        begin_tri = i;
        break
    end
end

% DETERMINE # VOLUMES AND TR FROM MR TRIGGERS (OPTIONAL)
trim_time_diff = trim_time(3:end)-trim_time(2:end-1);
tr = mean(trim_time_diff);
vols = size(trim_time,1);

figure;  plot(PETCO2);

%% PROCESSING
% TRIMMING OF END-TIDAL TIME COURSE BY INCORPORATING TIME-SHIFT (previously determined)
% Identify index of 1st volume (trigger time - shift) 
% Then removal of volumes before the 1st one
trim_time_begin_del = trim_time(begin_tri,:)-shift; % 1st trigger - time shift (all in seconds)
trim_time_end_del = trim_time(end,:)-shift; % last trigger - time shift (all in seconds)
idx_begin = find(timings >= trim_time_begin_del); % find index of (1st trigger - shift)
idx_end = find(timings >= trim_time_end_del); % find index of (1st trigger - shift)
trim_PETCO2 = PETCO2(idx_begin(1):idx_end(1)); % remove earlier volumes
trim_timings = timings(idx_begin(1):idx_end(1));

% FIND PEAKS
[PKS,LOCS]=findpeaks(trim_PETCO2, 'MinPeakDistance', MinPeakDistVal, 'MinPeakProminence',MinPeakPromVal); % values tested by trial and error
LOCS (PKS < MinThresholdVal) = [];
PKS (PKS < MinThresholdVal) = [];

figure;  plot(trim_timings,trim_PETCO2); hold on; plot(trim_timings(LOCS), PKS,'r*');

% INTERPOLATION TO TR
xi = linspace(trim_timings(LOCS(1)), trim_timings(LOCS(end)), vols);
ev_co2 = interp1(trim_timings(LOCS),PKS, xi);

% NORMALISE TIMECOURSE BETWEEN 0 and 1 TO CREATE EV 
ev_co2_norm=(ev_co2-min(ev_co2))/(max(ev_co2)-min(ev_co2));      
out_co2=round(ev_co2_norm,3,'significant');
figure; plot(out_co2); hold on

% CONVERT TO mmHg
mmHg = air_pres/1.33322387415; % pressure bar
co2_mmHg = (ev_co2.*mmHg)./100; % div by 100 as values are in percent

% SAVE text file for the co2 EV
cd (output_folder)
fid=fopen('ev_co2.txt','Wt');
fprintf(fid,'%f \r',out_co2); 
fclose(fid);

%% UNCOMMENT if you want to extract normocapnia and hypercapnia average values
% % Protocol details (set this accordingly)
% baseline = 60; % Length of initial baseline block in seconds
% blocksizeon = 120; % Length of ON block in seconds
% blocksizeoff = 120; % Length of OFF block in seconds
% nblocks = 2;
% 
% % CONVERT TO NUMBER OF VOLUMES
% baseline_tr = baseline/tr;
% blocksizeon_tr = blocksizeon/tr;
% blocksizeoff_tr = blocksizeoff/tr;        
% nblocks = 2;
% 
% % CREATE PARADIGM VECTOR
% vec_co2 = zeros(1,round(baseline_tr)+round(physio_shift/tr));
% for i = 1:nblocks
% vec_co2 = [vec_co2 ones(1,round(blocksizeon_tr)) zeros(1,round(blocksizeoff_tr))];
% end
% vec_co2 = vec_co2(1:vols);
% 
% plot(vec_co2,'r'); hold off
% 
% normocap = mean(co2_mmHg(vec_co2==0));
% hypercap = mean(co2_mmHg(vec_co2==1));
% 
% % SAVE
% cd (output_folder)
% fid=fopen('end_tidals.txt','Wt');
% fprintf(fid,'Normocapnia average in mmHg %f \r',normocap); 
% fprintf(fid,'Hypercapnia value in mmHg %f \r',hypercap);
% fclose(fid);

%
end