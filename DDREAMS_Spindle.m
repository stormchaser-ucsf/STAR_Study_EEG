%% DDREAMS spindle analyses


clc;clear
close all

% adding paths
if ispc
    addpath(genpath('C:\Users\nikic\Documents\MATLAB\Ana EEG\hdEEG_Trauma_PTSD\SAGA_Matlab\SAGA_interface'))
    %root_path = 'F:\DATA\EEG Data Anne Richards\STAR\STAR_Pilot_data-selected\Pilot18';
    root_path='F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Participant 24028';
    cd(root_path)
    addpath('C:\Users\nikic\Documents\MATLAB\eeglab2023.1')
    %addpath(genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114'))
    addpath('C:\Users\nikic\Documents\MATLAB')


elseif isunix
    addpath(genpath('/home/user/Documents/Repositories/STAR_Study_EEG/SAGA_Matlab/SAGA_interface'))
    addpath('/home/user/Documents/MATLAB/eeglab2023.1')
    %addpath(genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114'))
    addpath('/home/user/Documents/MATLAB')


end

% load EDF
cd('/media/user/Data/Ana EEG/DDREAMS/DDREAMS_EDF_2026_v2-selected/')
eeglab

% three channels
% LEOG: reference to Fz
% REOG: reference to Fz
% E3: Fz which is difference of L and R EOG
% done

data = EEG.data;
data = data(1:3,:);


% filters
Fs=EEG.srate;


% low pass filters
lpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',30, ...
    'SampleRate',Fs);

bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',65,'HalfPowerFrequency2',100, ...
    'SampleRate',Fs);

spFilt1 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',11,'HalfPowerFrequency2',13, ...
    'SampleRate',Fs);

spFilt2 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',16, ...
    'SampleRate',Fs);

soFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.16,'HalfPowerFrequency2',1.25, ...
    'SampleRate',Fs);

deltaFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',4, ...
    'SampleRate',Fs);

% tmp stuff for spindle filtering
filt_data = filtfilt(spFilt1,data(2,:));
spn_pow = abs(hilbert(filt_data));
figure;
plot(smooth(spn_pow,100))

% so filt stuff
so_data = filtfilt(soFilt,data(2,:));
so_pow = abs(hilbert(so_data));
figure;
plot(smooth(so_pow,100))




%detect spindles
I = sleep_staging>0;
grid_sp = detect_spindles(data(1:6,:),I,soFilt,spFilt1,spFilt2,sleep_staging,Fs);


%%%


