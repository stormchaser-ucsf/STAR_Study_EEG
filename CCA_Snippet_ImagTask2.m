% goal here is to look at low dimensional representations shared between
% snippet viewing and imagination task

%% load the subject matched data across both tasks


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
    %addpath(genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114'))
    addpath('/home/user/Documents/MATLAB')
    root_path='/media/user/Data/Ana EEG/STAR/Phase 2/Imagine Part 2 Preprocessed/';
end

cd(root_path)
eeglab

files=findfiles('.set',root_path)'; % all the imag2 files

%Splitting into neutral and trauma snippets
neutral_epochs= [2 4 6 7 9 11 13 14 16 18 19 22 23 25 27 30 31 33 35 38];
idx_neutral = zeros(40,1);
idx_neutral(neutral_epochs)=1;
trauma_epochs = find(idx_neutral==0);
idx_trauma = zeros(40,1);
idx_trauma(trauma_epochs)=1;
idx_trauma = logical(idx_trauma);
idx_neutral = logical(idx_neutral);

chdata_neutral=[];
chdata_trauma=[];
k=1;
error_files=[];
no_epochs_rej=[];


subj_loaded ={};
for i=1:length(files)


    try
        EEG = pop_loadset('filename',files{i});
        eeglab redraw
        file_loaded=true;
        subj_loaded = cat(1,subj_loaded,files{i}(end-32:end-28));
    catch
        file_loaded=false;
        error_files=[error_files i];
    end

    if file_loaded

        % bad trial rejection
        EEG = pop_autorej(EEG, 'nogui','on','threshold',200,'eegplot','on'); % automatic bad trial rejection
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

        % Get Trauma ERPs
        data_trauma = data(:,:,idx_trauma);
        data_trauma = nanmean(data_trauma,3);
        chdata_trauma = cat(3,chdata_trauma,data_trauma);

        % Get Neutral ERPs
        data_neutral = data(:,:,idx_neutral);
        data_neutral = nanmean(data_neutral,3);
        chdata_neutral = cat(3,chdata_neutral,data_neutral);

        %
        % chdata1=[];
        % %chdata1 = squeeze(cat(3,chdata1,data(6,:,idx_trauma),data(36,:,idx_trauma),data(37,:,idx_trauma)));
        % chdata1 = squeeze(cat(3,chdata1,data(28,:,idx_trauma),data(29,:,idx_trauma),data(30,:,idx_trauma),...
        %     data(27,:,idx_trauma),data(54,:,idx_trauma),data(55,:,idx_trauma)));
        % %chdata1 = squeeze(cat(3,chdata1,data(15,:,idx_trauma),data(43,:,idx_trauma),data(44,:,idx_trauma),...
        % %     data(47,:,idx_trauma),data(40,:,idx_trauma)));
        %
        % % ac1 = sum(abs(chdata1)>150);
        % % ac1 = find(ac1~=0);
        % % chdata1(:,ac1)=NaN;
        % chdata_trauma(k,:) = smooth(nanmean(chdata1,2),50);
        %
        % % Get neutral ERPs
        % chdata2=[];
        % %chdata2 = squeeze(cat(3,chdata2,data(6,:,idx_neutral),data(36,:,idx_neutral),data(37,:,idx_neutral)));
        % chdata2 = squeeze(cat(3,chdata2,data(28,:,idx_neutral),data(29,:,idx_neutral),data(30,:,idx_neutral),...
        %     data(27,:,idx_neutral),data(54,:,idx_neutral),data(55,:,idx_neutral)));
        % % chdata2 = squeeze(cat(3,chdata2,data(15,:,idx_neutral),data(43,:,idx_neutral),data(44,:,idx_neutral),...
        % %     data(47,:,idx_neutral),data(40,:,idx_neutral)));
        %
        % % ac2 = sum(abs(chdata2)>150);
        % % ac2 = find(ac2~=0);
        % % chdata2(:,ac2)=NaN;
        % chdata_neutral(k,:) = smooth(nanmean(chdata2,2),50);

        if i<length(files)
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
            eeglab redraw
        end
        k=k+1;

    end
end

chdata_neutral_imag = chdata_neutral;
chdata_trauma_imag = chdata_trauma;
subj_loaded_imag=subj_loaded;

delete(gcp)
do_stats_and_plot_imag(chdata_neutral_imag,chdata_trauma_imag,EEG);

% some time frequency decomposition
x = squeeze(nanmean(chdata_trauma_imag(19,:,:),3));
[pxx, f] = pwelch(x(1e3:5e3), 512, 256, 1024, EEG.srate);
[pxx1, f1] = pwelch(x(8e3:13e3), 512, 256, 1024, EEG.srate);
figure;
hold on
plot(f,log(pxx))
plot(f1,log(pxx1))
xlim([0 50])

x = squeeze(nanmean(chdata_trauma_imag(19,:,:),3)) - ...
    squeeze(nanmean(chdata_neutral_imag(19,:,:),3));
[pxx, f] = pwelch(x(1e3:5e3), 512, 256, 1024, EEG.srate);
[pxx1, f1] = pwelch(x(8e3:13e3), 512, 256, 1024, EEG.srate);
figure;
hold on
plot(f,log(pxx))
plot(f1,log(pxx1))
xlim([0 50])
legend({'open','closed'})


% clear all
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
eeglab redraw

%%%%% now getting data from snippet viewing
root_path='/media/user/Data/Ana EEG/STAR/Phase 2/Snippet_Datasets_Processed';
files=findfiles('.set',root_path)';

files1={};
for i=1:length(files)
    for j=1:length(subj_loaded_imag)
        a  = (regexp(files{i},subj_loaded_imag{j}));
        if ~isempty(a)
            files1=[files1;files(i)];
            break
        end
    end
end
files=files1;

%Splitting into neutral and trauma snippets
neutral_epochs= [1 3 5 7 8 11 14 17 19 21 23 24 27 29 30 33 34 37 39 40];
idx_neutral = zeros(40,1);
idx_neutral(neutral_epochs)=1;
trauma_epochs = find(idx_neutral==0);
idx_trauma = zeros(40,1);
idx_trauma(trauma_epochs)=1;
idx_trauma = logical(idx_trauma);
idx_neutral = logical(idx_neutral);

chdata_neutral=[];
chdata_trauma=[];
k=1;
error_files=[];
no_epochs_rej=[];
subj_loaded={};
for i=1:length(files)


    try
        EEG = pop_loadset('filename',files{i});
        eeglab redraw
        file_loaded=true;
        subj_loaded = cat(1,subj_loaded,files{i}(end-19:end-15));
    catch
        file_loaded=false;
        error_files=[error_files i];
    end

    if file_loaded

        % bad trial rejection
        EEG = pop_autorej(EEG, 'nogui','on','threshold',200,'eegplot','on'); % automatic bad trial rejection
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

        % Get Trauma ERPs
        data_trauma = data(:,:,idx_trauma);
        data_trauma = nanmean(data_trauma,3);
        chdata_trauma = cat(3,chdata_trauma,data_trauma);

        % Get Neutral ERPs
        data_neutral = data(:,:,idx_neutral);
        data_neutral = nanmean(data_neutral,3);
        chdata_neutral = cat(3,chdata_neutral,data_neutral);

        %
        % chdata1=[];
        % %chdata1 = squeeze(cat(3,chdata1,data(6,:,idx_trauma),data(36,:,idx_trauma),data(37,:,idx_trauma)));
        % chdata1 = squeeze(cat(3,chdata1,data(28,:,idx_trauma),data(29,:,idx_trauma),data(30,:,idx_trauma),...
        %     data(27,:,idx_trauma),data(54,:,idx_trauma),data(55,:,idx_trauma)));
        % %chdata1 = squeeze(cat(3,chdata1,data(15,:,idx_trauma),data(43,:,idx_trauma),data(44,:,idx_trauma),...
        % %     data(47,:,idx_trauma),data(40,:,idx_trauma)));
        %
        % % ac1 = sum(abs(chdata1)>150);
        % % ac1 = find(ac1~=0);
        % % chdata1(:,ac1)=NaN;
        % chdata_trauma(k,:) = smooth(nanmean(chdata1,2),50);
        %
        % % Get neutral ERPs
        % chdata2=[];
        % %chdata2 = squeeze(cat(3,chdata2,data(6,:,idx_neutral),data(36,:,idx_neutral),data(37,:,idx_neutral)));
        % chdata2 = squeeze(cat(3,chdata2,data(28,:,idx_neutral),data(29,:,idx_neutral),data(30,:,idx_neutral),...
        %     data(27,:,idx_neutral),data(54,:,idx_neutral),data(55,:,idx_neutral)));
        % % chdata2 = squeeze(cat(3,chdata2,data(15,:,idx_neutral),data(43,:,idx_neutral),data(44,:,idx_neutral),...
        % %     data(47,:,idx_neutral),data(40,:,idx_neutral)));
        %
        % % ac2 = sum(abs(chdata2)>150);
        % % ac2 = find(ac2~=0);
        % % chdata2(:,ac2)=NaN;
        % chdata_neutral(k,:) = smooth(nanmean(chdata2,2),50);


         if i<length(files)
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
            eeglab redraw
        end
        k=k+1;

    end
end

delete(gcp)
do_stats_and_plot(chdata_neutral,chdata_trauma,EEG);


% some time frequency decomposition
x = squeeze(nanmean(chdata_trauma(27,:,:),3));
x1 = squeeze(nanmean(chdata_neutral(27,:,:),3));
[pxx, f] = pwelch(x(2e3:5e3), 512, 256, 1024, EEG.srate);
[pxx1, f1] = pwelch(x1(2e3:5e3), 512, 256, 1024, EEG.srate);
figure;
hold on
plot(f,log(pxx))
plot(f1,log(pxx1))
legend({'Trauma','Neutral'})
xlim([0 50])

x = squeeze(nanmean(chdata_trauma(19,:,:),3)) - ...
    squeeze(nanmean(chdata_neutral(19,:,:),3));
[pxx, f] = pwelch(x(2e3:5e3), 512, 256, 1024, EEG.srate);
figure;
hold on
plot(f,log(pxx))
xlim([0 50])



%% CCA leave one subject out 

Fs=1e3;
bpFilt = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 30, 'SampleRate', Fs);

res=[];top=[];
for i=1:size(chdata_trauma,3) %leave on subject out
    test_idx=i;
    aa=ones(22,1);
    aa(test_idx)=0;
    train_idx = find(aa==1);
    %train_idx=[1:14 16:22];
    %test_idx=15;
    %ch_idx = [14 43 15 44 16 46 19 47 20 48 23 50 51 24 25 53:56 27];
    ch_idx=1:62;

    idx_snippet = 2e3:4e3;
    idx_imag = 9e3:11e3;

    snippet = chdata_trauma - chdata_neutral;
    Xa = snippet(ch_idx,idx_snippet,train_idx);
    %
    % for i=1:size(Xa,3)
    %     temp = squeeze(Xa(:,:,i));
    %     temp = temp-mean(temp);
    %     Xa(:,:,i)=temp;
    % end
    Xa= Xa(:,:);
    Xa=Xa';
    Xa= filtfilt(bpFilt,Xa);
    [c,s,l]= pca(Xa);
    Xa = s(:,1:15);



    Xa_test = snippet(ch_idx,idx_snippet,test_idx);
    Xa_test = Xa_test';
    Xa_test= filtfilt(bpFilt,Xa_test);
    %[c,s,l]= pca(Xa_test);
    %Xa_test = s(:,1:5);
    Xa_test = (Xa_test-mean(Xa_test))*c(:,1:15);

    imag_data = chdata_trauma_imag - chdata_neutral_imag;
    Xb = imag_data(ch_idx,idx_imag,train_idx);
    %
    % for i=1:size(Xb,3)
    %     temp = squeeze(Xb(:,:,i));
    %     temp = temp-mean(temp);
    %     Xb(:,:,i)=temp;
    % end

    Xb= Xb(:,:);
    Xb=Xb';
    Xb= filtfilt(bpFilt,Xb);
    [c,s,l]= pca(Xb);
    Xb = s(:,1:15);

    Xb_test = imag_data(ch_idx,idx_imag,test_idx);
    Xb_test = Xb_test';
    Xb_test= filtfilt(bpFilt,Xb_test);
    % [c,s,l]= pca(Xb_test);
    % Xb_test = s(:,1:5);
    Xb_test = (Xb_test-mean(Xb_test))*c(:,1:15);

    [Wa,Wb,S,Za,Zb] = cca(Xa,Xb);


    Za_cv = (Xa_test - mean(Xa_test)) * Wa;
    Zb_cv = (Xb_test - mean(Xb_test)) * Wb;


    figure;
    hold on
    plot(Za_cv(:,2))
    plot(Zb_cv(:,2))
    corrcoef(Za_cv(:,2),Zb_cv(:,2))

    %figure;topoplot(abs(Wb(:,1)),EEG.chanlocs)

    c=[];
    for i=1:size(Za_cv,2)
        x= corrcoef(Za_cv(:,i),Zb_cv(:,i));
        c(i) = x(1,2);
    end
    figure;
    plot(c)
    res=[res;c];

    [aa,bb]=max(c);
    top=[top bb];
    figure;
    hold on
    plot(Za_cv(:,bb))
    plot(Zb_cv(:,bb))
    corrcoef(Za_cv(:,bb),Zb_cv(:,bb))

    close all


end

figure;plot(mean(res,1))

%% CCA ON GRAND AVERAGED ERP    
% cross validate by leaving segments out 

Fs=1e3;
bpFilt = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 30, 'SampleRate', Fs);

ch_idx = [14 43 15 44 16 46 19 47 20 48 23 50 51 24 25 53:56 27];
%ch_idx=1:62;
idx_snippet = 2e3:4e3;
idx_imag = 9e3:11e3;
train_idx=1:22;

snippet = chdata_trauma - chdata_neutral;
Xa = snippet(ch_idx,idx_snippet,train_idx);
Xa = squeeze(mean(Xa,3))';
Xa = filtfilt(bpFilt,Xa);
% [c,s,l]= pca(Xa);
% Xa = s(:,1:5);

imag_data = chdata_trauma_imag - chdata_neutral_imag;
Xb = imag_data(ch_idx,idx_imag,train_idx);
Xb = squeeze(mean(Xb,3))';
Xb = filtfilt(bpFilt,Xb);
% [c,s,l]= pca(Xb);
% Xb = s(:,1:5);

train_times = [401:1200 1201:2001];
test_times = [1:400];
Xa_train = Xa(train_times,:);
Xb_train= Xb(train_times,:);
[Wa,Wb,S,Za,Zb] = cca(Xa_train,Xb_train);

Za_cv = (Xa(test_times,:)-mean(Xa(test_times,:)))*Wa;
Zb_cv = (Xb(test_times,:)-mean(Xb(test_times,:)))*Wb;

c=[];
for i=1:size(Za_cv,2)
    x= corrcoef(Za_cv(:,i),Zb_cv(:,i));
    c(i) = x(1,2);
end
figure;
plot(c)

figure;plot(Za_cv(:,1));hold on;plot(Zb_cv(:,1))

%% PCA whether the covariance structure is preserved over sig. electrodes

Fs=1e3;
bpFilt = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 30, 'SampleRate', Fs);


idx_snippet = 2e3:4e3;
idx_imag = 9e3:11e3;
train_idx=1:22;

snippet = chdata_trauma - chdata_neutral;
Xa = snippet(ch_idx,idx_snippet,train_idx);
Xa = squeeze(mean(Xa,3))';
%Xa = filtfilt(bpFilt,Xa);

imag_data = chdata_trauma_imag - chdata_neutral_imag;
Xb = imag_data(ch_idx,idx_imag,train_idx);
Xb = squeeze(mean(Xb,3))';
%Xb = filtfilt(bpFilt,Xb);

[c,s,l]= pca(Xa);
[c1,s1,l1]= pca(Xb);

[angles,Va,Vb] = principal_angles(c(:,1:5),c1(:,1:5));
angles


[Wa,Wb,S,Za,Zb] = cca(s(:,1:5),s1(:,1:5));
figure;plot(S)

% AR(p) process to estimate the stats of the cca corr?
% or using TME method to do stats end to end. 



