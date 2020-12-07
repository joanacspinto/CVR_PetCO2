% Joana Pinto, December 2020, Oxford
% Adapted from Daniel Bulte 2018 script

% This function is called in CVR_compute.m
% This function processes the PetCO2 timeseries to be used as regressor, including adding a time shift
% INPUT: physfile - physio file with PETCO2 raw data; brain_del - time shift
% OUTPUT: out_co2 - PETCO2 processed; normocap and hypercap - average PETCO2 values during those periods, respectively

function [out_co2, normocap,hypercap] = PETCO2_proc(physfile, brain_del)
addpath('/usr/local/fsl/etc/matlab');

%% SET PARAMETERS

% Protocol (this depends on the protocol)
baseline = 60; % Length of initial baseline block in seconds
blocksizeon = 120; % Length of ON block in seconds
blocksizeoff = 120; % Length of OFF block in seconds
vols = []; % possibly get this below from triggers (OPTIONAL)
tr = []; % possibly get this below from triggers (OPTIONAL)

% CO2 components 
samp_rate = 100; % Powerlab sampling rate (Hz)
air_pres = 1020; % Enter barometric pressure in mbars:

% Thresholds and delays
threshold_trig = 3; % threshold to detect triggers 
brain_del = 0;
delay = 15+brain_del; % mechanical delay (15) + time shifts in regressor

%% LOAD DATA

% UNCOMMENT LINE FOR DEBUGGING
% physfile = load('/Users/joana/Documents/Work_Oxford/Data_extra/OHBA/001/2018_119_001.txt'); % Path to the physiological data file

phys = physfile; % Load physiological data file

timings = phys(:,1);
PETCO2 = phys(:,2); % Second column is the PETCO2 values
PETO2 = phys(:,3); % Third column is the PETO2 values
trig = phys(:,4); % MR triggers

% DETERMINE # VOLUMES AND TR FROM MR TRIGGERS (OPTIONAL)
trim_time = timings(trig>threshold_trig); 
trim_time_diff = trim_time(3:end)-trim_time(2:end-1);
tr = mean(trim_time_diff);
vols = size(trim_time,1)-1;

% CONVERT TO NUMBER OF VOLUMES
baseline = baseline/tr;
blocksizeon = blocksizeon/tr;
blocksizeoff = blocksizeoff/tr;        
     
% TEMPORAL SHIFT OF END-TIDAL TIME COURSES BY 15 SECONDS (SAMPLE DELAY) +
% TIME SHIFT
trim_time_begin = trim_time(2,:); % 1st trigger
trim_time_begin_del = (trim_time_begin-delay); % 1st trigger (in sec) - delay
idx = find(timings == trim_time_begin_del); % find index of 1st trigger (in sec) - delay
trim_PETCO2 = PETCO2(idx:end);

% DETERMINE RESPIRATORY FREQUENCY DURING BASELINE AND USES THIS INFO TO 
% DETERMINE THE SIZE OF THE END-TIDAL SEARCH WINDOW

samp_period=1/samp_rate;
length=baseline*tr*samp_rate;
temper = fft(trim_PETCO2(1:length));
P2 = abs(temper/length);
P1 = P2(1:length/2+1);
P1(2:end-1)=2*P1(2:end-1);
f = samp_rate*(0:(length/2))/length;

[pk,loc] = max(P1(2:end));

pkloc = loc+1;

harm = f(pkloc);
resp_period = round(1/harm); % e.g. 8s

nsearch=(resp_period+1)*samp_rate; % search window = 1 second more than the respiratory period
windows = floor(size(trim_PETCO2,1)/nsearch); 

%% Find peaks
% CO2 trace
etc=zeros(windows);
k=1;
for i=1:windows
    localmax=trim_PETCO2((i-1)*nsearch+1);
    posmax=(i-1)*nsearch+1;
    for j=1:nsearch
        if trim_PETCO2((i-1)*nsearch+j) > localmax 
            localmax = trim_PETCO2((i-1)*nsearch+j);
            posmax=(i-1)*nsearch+j;
        end
    end  
    etc(k,1)=posmax;
    etc(k,2)=localmax;
    i=i+1;
    k=k+1;
end

% make new full sample ET time course
resamp_co2 = zeros(size(trim_PETCO2,1),1);
for x=1:size(etc,1)-1
    dist_c = etc(x+1,1) - etc(x,1);
    step_c = etc(x+1,2) - etc(x,2);
    ramp_c = step_c./dist_c;
    for g=0:dist_c
        resamp_co2(etc(x,1)+g) = etc(x,2)+(ramp_c.*g);
    end    
end   

% pad the start and end with repeats of first and last et value to maintain
% length and phase
resamp_co2(1:etc(1,1)-1)=resamp_co2(etc(1,1));
resamp_co2(etc(size(etc,1),1):size(trim_PETCO2,1))= resamp_co2(etc(size(etc,1),1));

% make new time course at the TR resolution, then normalise output
block = round(tr.*samp_rate);
ev_co2=zeros(vols,1);

for i=1:vols
    ev_co2(i)=resamp_co2(block.*i);
end

% create a timecourse of the end tidal CO2 values at the TR's for use with CVR sigmoids
co2_tc = ev_co2;

% convert to mmHg
mmHg = air_pres/1.33322387415; % pressure bar
co2_mmHg = (co2_tc.*mmHg)./100; % div by 100 as values are in percent

% average all of first baseline block
normocap = mean(co2_mmHg(1:baseline+delay));

% select 2nd half of each hypercapnic block to average
hyperblock = cat(1,co2_mmHg(baseline+delay+blocksizeon/2:baseline+delay+blocksizeon),co2_mmHg(baseline+delay+blocksizeon+blocksizeoff+blocksizeon/2:baseline+delay+blocksizeon+blocksizeoff+blocksizeon));
hypercap = mean(hyperblock);

% normalise timecourse betwwen 0 and 1 to create EV
ev_co2=(ev_co2-min(ev_co2))/(max(ev_co2)-min(ev_co2)); 

%% plot results

figure;
plot(trim_PETCO2);
hold on
plot(resamp_co2,'Color','red');

figure
plot(ev_co2);

% round to 2 sig figs
out_co2=round(ev_co2,2,'significant');

%% write text files for the o2 and co2 EV's
% fid=fopen('ev_co2.txt','Wt');
% fprintf(fid,'%f \r',out_co2); 
% fclose(fid);
% 
% fid=fopen('co2_tc.txt','Wt');
% fprintf(fid,'%f \r',co2_tc); 
% fclose(fid);
% 
% fid=fopen('end_tidals.txt','Wt');
% fprintf(fid,'Normocapnia average in mmHg %f \r',normocap); 
% fprintf(fid,'Hypercapnia value in mmHg %f \r',hypercap);
% fclose(fid);

end