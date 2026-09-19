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
    addpath(genpath('/home/user/Documents/Repositories/STAR_Study_EEG/'))
    addpath('/home/user/Documents/MATLAB/eeglab2023.1')
    %addpath(genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114'))
    %addpath('/home/user/Documents/MATLAB')


end

%filepath
eeg_filepath = '/media/user/Data/Ana EEG/DDREAMS/DDREAMS_EDF_2026_v2-selected';
filepath = '/media/user/Data/Ana EEG/DDREAMS/DDREAMS_EBE Sleep Stages-selected/20001_1-3_0000334000214_EpochByEpochDetail/';

% load EDF
cd('/media/user/Data/Ana EEG/DDREAMS/DDREAMS_EDF_2026_v2-selected/')
eeglab

% load
eeg_filename = '20001_N1_v2.0.edf';
eeg_filename = fullfile(eeg_filepath,eeg_filename);
EEG = pop_biosig(eeg_filename);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
eeglab redraw

% three channels
% LEOG: reference to Fz
% REOG: reference to Fz
% E3: Fz which is difference of L and R EOG
% done

data = EEG.data;
data = data(1:3,:);

% duration in seconds
dur = size(data,2)/EEG.srate;


% load CSV for sleep staging 

%filename = "/media/user/Data/Ana EEG/DDREAMS/DDREAMS_EBE Sleep Stages-selected/20001_1-3_0000334000214_EpochByEpochDetail/20001_20260124_205200____20210809_122452_0000334000214_EpochByEpochDetail.csv";
filename = "20001_N1.csv";
filename = fullfile(filepath,filename);
sleep_scores = load_csv(filename);

% for synchronization, this is the key:
% each epoch from csv is the start of a 30s window, and that 30s will have
% the associated sleep stage. Build a sleep staging score based on that.
% The EEG data will usually have some extra bit (<30s) at the end which you
% can discard. 

p = gcp('nocreate');
if isempty(p)
    parpool('threads');
end

sleep_staging=[];
parfor i = 1:length(sleep_scores.StageFinal)    
    tmp = sleep_scores.StageFinal(i);
    tmp = repmat(tmp,30*EEG.srate,1);
    sleep_staging = [sleep_staging;tmp];
end
figure;
plot(sleep_staging)

data = data(:,1:length(sleep_staging));


I = zeros(size(sleep_staging));
I(sleep_staging==2) =1;
I(sleep_staging==3) =1;

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

% notch filtering line noise
notchFilt = designfilt('bandstopiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', 59, ...
    'HalfPowerFrequency2', 61, ...
    'SampleRate', Fs);
data = filtfilt(notchFilt,data')';

% band pass filtering eeg data in 0.1-30Hz range
data = filtfilt(lpFilt,data')';

% tmp stuff for spindle filtering
% filt_data = filtfilt(spFilt1,data(2,:));
% spn_pow = abs(hilbert(filt_data));
% figure;
% plot(smooth(spn_pow,100))
% 
% % so filt stuff
% so_data = filtfilt(soFilt,data(2,:));
% so_pow = abs(hilbert(so_data));
% figure;
% plot(smooth(so_pow,100))

%detect spindles
% have to ignore the spike every 15 min so have to set those windows down
% to zero 
% get background SO epochs around spindles 

grid_sp = ...
    detect_spindles_ddreams(data(2,:),I,soFilt,spFilt1,spFilt2,sleep_staging,Fs);

% detect slow oscillations
grid_so = SO_analyses(data(2,:),I,soFilt,spFilt1,spFilt2,sleep_staging,0,Fs);


%%%


