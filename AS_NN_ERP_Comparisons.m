%% COMPARISON W/ ANTHONY AND NN PROCESSED DATASETS

%% init

clc;clear
close all


% adding paths
if ispc
    addpath(genpath('C:\Users\nikic\Documents\MATLAB\Ana EEG\hdEEG_Trauma_PTSD\SAGA_Matlab\SAGA_interface'))
    %root_path = 'F:\DATA\EEG Data Anne Richards\STAR\STAR_Pilot_data-selected\Pilot18';
    root_path='F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Participant 24011';
    cd(root_path)
    addpath('C:\Users\nikic\Documents\MATLAB\eeglab2023.1')
    %addpath(genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114'))
    addpath('C:\Users\nikic\Documents\MATLAB')
else
    root_path='/media/user/Data/Ana EEG/STAR/Phase 2/Participant 24011';
    addpath('/home/user/Documents/MATLAB/eeglab2023.1');
    addpath(genpath('/home/user/Documents/MATLAB/Ana EEG/hdEEG_Trauma_PTSD/SAGA_Matlab/SAGA_interface/'))
    cd(root_path);

    %cd('/home/user/Documents/Repositories/STAR_Study_EEG/')
end

%% % anthony ERPs

cd('/media/user/Data/Ana EEG/STAR/Phase 2/Anthony snippet ICA practice')
files= findfiles('.set',pwd);
eeglab

data_as=[];
error_files=[];
no_epochs_rej=[];

neutral_epochs= [1 3 5 7 8 11 14 17 19 21 23 24 27 29 30 33 34 37 39 40];
idx_neutral = zeros(40,1);
idx_neutral(neutral_epochs)=1;
trauma_epochs = find(idx_neutral==0);
idx_trauma = zeros(40,1);
idx_trauma(trauma_epochs)=1;
idx_trauma = logical(idx_trauma);
idx_neutral = logical(idx_neutral);

chdata_neutral_as=[];
chdata_trauma_as=[];
k=1;


for i=1:length(files)-1
    try
        EEG = pop_loadset('filename',files{i});
        eeglab redraw
        file_loaded=true;
    catch
        file_loaded=false;
        error_files=[error_files i];
    end

    if file_loaded

        % bad trial rejection
        try
            EEG = pop_autorej(EEG, 'nogui','on','threshold',200,'eegplot','on'); % automatic bad trial rejection
        catch
            EEG = pop_interp(EEG, [13], 'spherical');
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
            eeglab redraw
            EEG = pop_autorej(EEG, 'nogui','on','threshold',200,'eegplot','on'); % automatic bad trial rejection
        end


        bad_trials = find(EEG.reject.rejauto);% bad trials
        close
        data=EEG.data;
        data(:,:,bad_trials) = NaN;
        no_epochs_rej = [no_epochs_rej length(bad_trials)];

        % detrend each channel's epochs remove baseline from each channel's
        % epochs
        for j=1:size(data,1)

            % detrend
            tmp = squeeze(data(j,:,:));
            tmp = detrend(tmp);

            % remove baseline
            m = nanmean(tmp(1:1000,:),1);
            tmp = tmp - m;

            % store
            data(j,:,:) = tmp;
        end

        data_as = cat(4,data_as,data);

        % Get Trauma ERPs
        data_trauma = data(:,:,idx_trauma);
        data_trauma = nanmean(data_trauma,3);
        chdata_trauma_as = cat(3,chdata_trauma_as,data_trauma);

        % Get Neutral ERPs
        data_neutral = data(:,:,idx_neutral);
        data_neutral = nanmean(data_neutral,3);
        chdata_neutral_as = cat(3,chdata_neutral_as,data_neutral);

        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        eeglab redraw
    end
end

%% NN ERPs

data_nn=[];
error_files=[];
no_epochs_rej=[];

chdata_trauma_nn=[];
chdata_neutral_nn=[];

filepath_nn = '/media/user/Data/Ana EEG/STAR/Phase 2/Snippet_Datasets_Processed';
cd(filepath_nn);
filenames=[];
for i=1:length(files)-1
    tmp  = files{i};
    tmp = tmp(end-29:end-25);
    f = findfiles(tmp,pwd)';
    filename = f{2};
    filenames = [filenames ;filename];

    try
        EEG = pop_loadset('filename',filename);
        eeglab redraw
        file_loaded=true;
    catch
        file_loaded=false;
        error_files=[error_files i];
    end

    if file_loaded

        % bad trial rejection
        %try
        EEG = pop_autorej(EEG, 'nogui','on','threshold',200,'eegplot','on'); % automatic bad trial rejection
        % catch
        %     EEG = pop_interp(EEG, [13], 'spherical');
        %     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
        %     eeglab redraw
        %     EEG = pop_autorej(EEG, 'nogui','on','threshold',200,'eegplot','on'); % automatic bad trial rejection
        %end


        bad_trials = find(EEG.reject.rejauto);% bad trials
        close
        data=EEG.data;
        data(:,:,bad_trials) = NaN;
        no_epochs_rej = [no_epochs_rej length(bad_trials)];

        % detrend each channel's epochs remove baseline from each channel's
        % epochs
        for j=1:size(data,1)

            % detrend
            tmp = squeeze(data(j,:,:));
            tmp = detrend(tmp);

            % remove baseline
            m = nanmean(tmp(1:1000,:),1);
            tmp = tmp - m;

            % store
            data(j,:,:) = tmp;
        end

        data_nn = cat(4,data_nn,data);

         % Get Trauma ERPs
        data_trauma = data(:,:,idx_trauma);
        data_trauma = nanmean(data_trauma,3);
        chdata_trauma_nn = cat(3,chdata_trauma_nn,data_trauma);

        % Get Neutral ERPs
        data_neutral = data(:,:,idx_neutral);
        data_neutral = nanmean(data_neutral,3);
        chdata_neutral_nn = cat(3,chdata_neutral_nn,data_neutral);

        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        eeglab redraw
    end
end

%% plotting ERPs for the two

f1= squeeze(nanmean(data_as,4));
f2= squeeze(nanmean(data_nn,4));

f1 = nanmean(f1,3);
f2 = nanmean(f2,3);

figure;
plot(f1(57,:));
hold on
plot(f2(57,:));
legend({'AS','NN'})


n_nn = nanmean(chdata_neutral_nn,3);
t_nn = nanmean(chdata_trauma_nn,3);

n_as = nanmean(chdata_neutral_as,3);
t_as = nanmean(chdata_trauma_as,3);

figure;
subplot(2,1,1)
hold on
plot(n_nn(28,:))
plot(t_nn(28,:))
legend('Neutral NN','Trauma NN')
subplot(2,1,2)
hold on
plot(n_as(28,:))
plot(t_as(28,:))
legend('Neutral AS','Trauma AS')
sgtitle('Ch 28 Occipital')


figure;
ch=6;
subplot(2,1,1)
hold on
plot(n_nn(ch,:))
plot(t_nn(ch,:))
legend('Neutral NN','Trauma NN')
subplot(2,1,2)
hold on
plot(n_as(ch,:))
plot(t_as(ch,:))
legend('Neutral AS','Trauma AS')
sgtitle('Ch 6 Frontal')
