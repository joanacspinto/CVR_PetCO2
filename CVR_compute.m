% Joana Pinto, December 2020, Oxford

% This function computes CVR amplitude (and delay) measurements 
% using two different GLM strategies
% INPUT: filtered_func_data.nii.gz (BOLD pre processed), physio txt file with PETCO2 raw data,
% shift_start in sec (first shift, e.g. -10) and shift_end in sec 
% (last shift, e.g. 10) (e.g. -10:10)
% OUTPUT: 1. WITHOUT temporal optimization
%           1.1 CVR amplitude map 
%        2. WITH temporal optimization
%           2.1 CVR amplitude map 
%           2.2 CVR delay map       

function [CVRamp,CVRtime,CVR_amp_nonopt] = CVR_compute(BOLD_path, PETCO2_path, shift_start, shift_end)

%% Uncomment below for debugging
% BOLD_path = '/Users/joana/Documents/Work_Oxford/Data_extra/Sana_CVR/2018_119_nii/cvr.feat/filtered_func_data.nii.gz';
% PETCO2_path = '/Users/joana/Documents/Work_Oxford/Data_extra/Sana_CVR/2018_119_001.txt';
% shift_start = -10;
% shift_end = 10;
% 
%% LOAD DATA
filtered_func_data = niftiread(BOLD_path); % load BOLD
PETCO2 = load(PETCO2_path); % load PETCO2 (full, not trimmed)

%% CREATE BOLD BASELINE MAP
baseline = mean(filtered_func_data,4); % average in 4th dimension (time) - creates a baseline map
% other options can be used instead: average of a specific baseline period
% e.g. baseline = mean(filtered_func_data(:,:,:,baseline_period),4);

%% CREATE REGRESSOR
% In FEAT to add the model we select our regressor (e.g. PETCO2_EV) and
% usually we convolve it with a double gamma function (canonical
% haemodynamic response function - HRF)

t = 1:1:20; % TIME SCALE
DG_HRF = gampdf(t,6) + -.5*gampdf(t,10); % CANONIC HRF MODEL (double gamma)
DG_HRF = DG_HRF/max(DG_HRF); % SCALE HRF TO HAVE MAX AMPLITUDE OF 1

%% GLM - 2 approaches

% optimized shift by running several GLMs (amplitude and timing info)
shifts = shift_start:shift_end;
cope_opt = zeros(size(filtered_func_data,1), size(filtered_func_data,2), size(filtered_func_data,3),size(shifts,2));

for i = 1:size(shifts,2) % shifts
i  
[trim_PETCO2,normocap(i),hypercap(i)]= PETCO2_proc(PETCO2, shifts(i))
PETCO2_EV_conv = conv(trim_PETCO2,DG_HRF,'valid'); % convolve with HRF
PETCO2_EV_conv = (PETCO2_EV_conv -(min(PETCO2_EV_conv)))/(max(PETCO2_EV_conv )-min(PETCO2_EV_conv ));

nTRs = size(trim_PETCO2,1);

%% GLM (*insert BASIL code here*)
X = [];
% Generate the design matrix
X(:,1) = PETCO2_EV_conv; X(:,2)=ones(nTRs,1); % X(:,3)=linspace(1,nTRs,nTRs)'; %3rd column models linear drift (a typical scan artifact)
% plot design matrix (rescaled for plotting)
figure(4) 
% X(:,3) = X(:,3)/max(X(:,3));
colormap(gray);
image(X*64);

for x = 1:size(filtered_func_data,1)
    x
    for y =  1:size(filtered_func_data,2)
        for z = 1:size(filtered_func_data,3)
            
            Y1 = double(squeeze(filtered_func_data(x,y,z,1:nTRs)));
            
            if sum(Y1) > 0
         
            beta = X\Y1;
            cope_opt(x,y,z,i) = beta(1);
            model= X*beta;
            % residuals (x,y,z,i) = mean(Y1 - model);
                            
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

dPETCO2 = mean(hypercap-normocap);

%% CVR AMPLITUDE (NON OPTIMIZED)
cope_nonopt = cope_opt(:,:,:,find(shifts==0)); % shift zero is the non-optimized cope
PSC = (cope_nonopt./baseline).*100; % obtain percent signal change from baseline
CVR_amp_nonopt = PSC./dPETCO2; % map of PSC normalized by value of change in PETCO2 (in %/mmHg)

%% CVR AMPLITUDE AND DELAY (OPTIMIZED)
[copeopt_amp copeopt_time] = max(cope_opt,[],4); % find maximum cope and save value and index
PSC_opt = (copeopt_amp./baseline).*100;
CVR_amp = PSC_opt./dPETCO2; % map of PSC normalized by value of change in PETCO2 (in %/mmHg)

CVR_time=copeopt_time+shift_start-1; % transform indexes to shifts
CVR_time(isnan(PSC_opt))=NaN;

% end

%% PLOTS
% figure; plot(squeeze(filtered_func_data(45,45,50,:)));
% hold on; %plot(B_pred, 'k-','Linewidth', 3)
% plot(squeeze(B_pred2(i, 45,45,50,:)), 'r-','Linewidth', 3)

% CVR opt
figure (2);
for z = 1:size(filtered_func_data,3)
    subplot(8,8,z)
    imshow(imrotate(CVR_amp(:,:,z),90),[0 2],'Colormap', hot)
hold on
end

% CVD
figure (3);
for z = 1:size(filtered_func_data,3)
    subplot(8,8,z)
    imshow(imrotate(CVR_time(:,:,z),90),[shift_start shift_end],'Colormap', hot)
hold on
end

% CVR not opt
figure (4);
for z = 1:size(filtered_func_data,3)
    subplot(8,8,z)
    imshow(imrotate(CVR_amp_nonopt (:,:,z),90),[0 2],'Colormap', hot)
hold on
end

end
