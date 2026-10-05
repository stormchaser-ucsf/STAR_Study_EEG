% SLEEP REACTIVATION (MAIN) LOOPING OVER SUBJECTS

%% INIT
clc;clear
close all

% adding paths
if ispc
    addpath(genpath('C:\Users\nikic\Documents\MATLAB\Ana EEG\hdEEG_Trauma_PTSD\SAGA_Matlab\SAGA_interface'))
    root_path = 'F:\DATA\EEG Data Anne Richards\STAR\STAR_Pilot_data-selected\Pilot18';
    root_path='F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Snippet_Datasets_Processed';
else
    addpath(genpath('/home/user/Documents/Repositories/STAR_Study_EEG/SAGA_Matlab/SAGA_interface'))
    addpath(genpath('/home/user/Documents/Repositories/STAR_Study_EEG'))
    addpath('/home/user/Documents/MATLAB/eeglab2023.1')    
    addpath('/home/user/Documents/MATLAB')    
    root_path = '/media/user/Data/Ana EEG/STAR/Phase 2/';
end

% load dPC axes
cd('/media/user/Data/Ana EEG/STAR/Phase 2')
load dPCA_Shared_AllSubj
load dPCA_Shared_AllSubj_baseline
eeglab

%% SAVE SLEEP FILES PER SUBJECTS
Poly5toEEGlab

%% SUBJECTS

D = dir(pwd);
subj_folder={};
for j=1:length(D)
    if contains(D(j).name,'Participant')
        subj_folder = [subj_folder;fullfile(D(j).folder,D(j).name)];
    end
end
