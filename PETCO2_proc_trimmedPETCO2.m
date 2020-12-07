% Joana Pinto, December 2020, Oxford

% This function computes CVR amplitude (and delay) measurements
% using two different GLM strategies

% INPUT: filtered_func_data.nii.gz (BOLD pre-processed in FEAT), PETCO2 trimmed regressor (previously processed in MATLAB),
% shift_start in sec (first shift, e.g. -10), shift_end in sec (last shift, e.g. 10) (e.g. -10:10),
% dPETCO2 is the difference in PETCO2 during hypercapnia and baseline

% OUTPUT: 1. WITHOUT temporal optimization
%           1.1 CVR amplitude map
%         2. WITH temporal optimization
%           2.1 CVR amplitude map
%           2.2 CVR delay map

function [CVRamp,CVRtime,CVR_amp_nonopt] = CVR_compute(BOLD_path, PETCO2_path, shift_start, shift_end, dPETCO2)

%% Uncomment below for debugging
% BOLD_path = '/Users/joana/Documents/Work_Oxford/Data_extra/Sana_CVR/2018_119_nii/cvr.feat/filtered_func_data.nii.gz';
% PETCO2_path = load('/Users/joana/Documents/Work_Oxford/Data_extra/Sana_CVR/ev_co2.txt'); % load PETCO2 (trimmed)
% shift_start = -8;
% shift_end = 10;
% dPETCO2 = 5;

%% LOAD DATA
filtered_func_data = niftiread(BOLD_path); % load BOLD
PETCO2 = load(PETCO2_path);

%% CREATE BOLD BASELINE MAP
baseline = mean(filtered_func_data,4); % average in 4th dimension (time) - creates a baseline map
% other options can be used instead: average of a specific baseline period
% e.g. baseline = mean(filtered_func_data(:,:,:,baseline_period),4);

%% CREATE MAIN REGRESSOR (CONVOLUTION WITH HRF)
% In FEAT to add the model we select our regressor (e.g. PETCO2_EV) and
% usually we convolve it with a double gamma function (canonical
% haemodynamic response function - HRF)

t = 1:1:20; % TIME SCALE
DG_HRF = gampdf(t,6) + -.5*gampdf(t,10); % CANONIC HRF MODEL (double gamma)
DG_HRF = DG_HRF/max(DG_HRF); % SCALE HRF TO HAVE MAX AMPLITUDE OF 1

PETCO2_EV_conv = conv(PETCO2,DG_HRF,'valid'); % convolve with HRF

%% GLM - 2 approaches

% optimized shift by running several GLMs (amplitude and timing info)
shifts = shift_start:shift_end;
cope_opt = zeros(size(filtered_func_data,1), size(filtered_func_data,2), size(filtered_func_data,3),size(shifts,2));

for i = 1:size(shifts,2) % shifts
    i
    aux_pet = PETCO2_EV_conv
    
    % shift if necessary
    if shifts(i)<0
        aux_pet(1:-shifts(i)) = []; % remove initial volumes if shift is negative
    elseif shifts(i)>0
        aux_pet(end+1-shifts(i):end) = []; % remove last volumes if shift is positive
    end
    
    aux_pet = (aux_pet -(min(aux_pet)))/(max(aux_pet )-min(aux_pet )); % normalization of regressor
    
    %% GLM (*insert BASIL code here*)
    
    nTRs = size(aux_pet,1);
    X = [];
    % Generate the design matrix
    X(:,1) = PETCO2_EV_conv; X(:,2)=ones(nTRs,1); % X(:,3)=linspace(1,nTRs,nTRs)'; %3rd column models linear drift (a typical scan artifact)
    
    % plot design matrix (rescaled for plotting)
    figure(9); colormap(gray); image(X*64);
    
    for x = 1:size(filtered_func_data,1)
        x
        for y =  1:size(filtered_func_data,2)
            for z = 1:size(filtered_func_data,3)
                
                Y1 = double(squeeze(filtered_func_data(x,y,z,1:nTRs)));
                
                if sum(Y1) > 0
                    
                    beta = X\Y1;
                    cope_opt(x,y,z,i) = beta(1);
                    
                    % Other ways to perform GLM
                    %             Estimate the beta vector and error variance i.e. Solve the GLM for beta values.
                    %             beta_hat=inv(X'*X)*X'*Y1
                    %             B_pred=X*beta_hat;
                    %             ts = double(ts); % turn into normal matlab numbers
                    %           % Or
                    %             [b,dev,stats] = glmfit(X,Y1,'normal');
                    %             B_pred2=X*b([1 2 4]);
                    %             cope1=b(1);
                    
                end
            end
        end
    end
end

%% CVR AMPLITUDE (NON OPTIMIZED)
cope_nonopt = cope_opt(:,:,:,find(shifts==0)); % shift zero is the non-optimized cope
PSC = (cope_nonopt./baseline).*100; % obtain percent signal change from baseline
CVR_amp_nonopt = PSC./dPETCO2; % map of PSC normalized by value of change in PETCO2 (in %/mmHg)

%% CVR AMPLITUDE AND DELAY (OPTIMIZED)
[copeopt_amp copeopt_time] = max(cope_opt,[],4); % find maximum cope and save value and index
PSC_opt = (copeopt_amp./baseline).*100;
CVR_amp = PSC_opt./dPETCO2; % map of PSC normalized by value of change in PETCO2 (in %/mmHg)

CVR_time=copeopt_time+shift_start-1; % transform indexes to time shifts
CVR_time(isnan(PSC_opt))=NaN;

%% PLOTS
nlines = round(size(filtered_func_data,3)/8);
% CVR opt
figure (2);
for z = 1:size(filtered_func_data,3)
    subplot(nlines,8,z)
    imshow(imrotate(CVR_amp(:,:,z),90),[0 2],'Colormap', hot)
    hold on
end

% CVD
figure (3);
for z = 1:size(filtered_func_data,3)
    subplot(nlines,8,z)
    imshow(imrotate(CVR_time(:,:,z),90),[0 shift_end],'Colormap', hot)
    hold on
end

% CVR not opt
figure (4);
for z = 1:size(filtered_func_data,3)
    subplot(nlines,8,z)
    imshow(imrotate(CVR_amp_nonopt(:,:,z),90),[0 2],'Colormap', hot)
    hold on
end


end
