

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
    root_path1='/media/user/Data/Ana EEG/STAR/Phase 2/Snippet_Datasets_Processed/';
end

cd(root_path)
load('/media/user/Data/Ana EEG/STAR/Phase 2/sig_ch.mat')
eeglab

files_imag=findfiles('.set',root_path)'; % all the imag2 files
files_snippet =findfiles('.set',root_path1)'; % all the snippet files

% get the common files 
imag_idx=[];snippet_idx=[];
for i=1:length(files_imag)
    subj = files_imag{i}(end-32:end-28);
    for j=1:length(files_snippet)
        if length(regexp(files_snippet{j},subj))>0
            imag_idx=[imag_idx i];
            snippet_idx = [snippet_idx j];
            break
        end
    end
end

files_imag = files_imag(imag_idx);
files_snippet = files_snippet(snippet_idx);

% LOADING IMAGINATION FILES
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
files=files_imag;
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

% LOAD SNIPPET FILES
% clear all
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
eeglab redraw

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
files=files_snippet;
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

cd('/media/user/Data/Ana EEG/STAR/Phase 2')
save Imag_Snippet_SubjMatched chdata_neutral chdata_trauma chdata_trauma_imag ...
    chdata_neutral_imag files_snippet files_imag -v7.3


%% INTRINSIC MANIFOLD SHARED

Fs=1e3;
bpFilt = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 15, 'SampleRate', Fs);

a = chdata_trauma - chdata_neutral;
a = squeeze(mean(a,3));

b = chdata_trauma_imag - chdata_neutral_imag;
b = squeeze(mean(b,3));

idx_snippet = 2e3:5e3;
idx_imag = 8.5e3:11.5e3;

a=a(:,idx_snippet)';
b=b(:,idx_imag)';

%a=a-mean(a);
%b=b-mean(b);

[c1,s1,l1] = pca(a);
[c2,s2,l2] = pca(b);

% Q1: what is the shared dimensionality? 
figure;stem(cumsum(l1)./sum(l1))
hold on
stem(cumsum(l2)./sum(l2))
legend({'Snippet','Imag'})
hline(0.9)
axis tight
figure;
topoplot(c2(:,4),EEG.chanlocs);
axis tight


%Q2: what are the principal angles betwen the two manifolds? More sig. than
%TME stats?
dim=11;
cd('/home/user/Documents/Repositories/STAR_Study_EEG')
addpath(genpath(pwd))
angles = principal_angles(c1(:,1:dim),c2(:,1:dim));
clear dataTensor
dataTensor(:,:,1) = a;
dataTensor(:,:,2) = b;
surr_type = 'surrogate-TC';
dataTensor = double(dataTensor);
maxEntropy = run_tme(dataTensor,surr_type);
parpool('threads')
boot_angles=[];
parfor loop=1:1000
    surrTensor = simulate_time(maxEntropy);
    tmp = compute_prin_angles_PC(surrTensor,dim);
    boot_angles(:,loop) = tmp;
end

figure;
plot(angles)
hold on
plot(boot_angles,'Color',[.2 .2 .2 .02])

%Q3: cross validation ie VAF when generalizing to held out subject
vaf=[];
for i=1:size(chdata_neutral,3)
    I = ones(size(chdata_neutral,3),1);
    test_idx=i;
    I(i)=0;
    train_idx = find(I==1);

    % cross validating first on snippet viewing (how well each subject
    % matches to manifold estimated from mean of other subjects)
    a = chdata_trauma_imag(:,:,train_idx) - chdata_neutral_imag(:,:,train_idx);
    a = squeeze(mean(a,3));
    a=a(:,idx_imag)';
    [c,s,l] = pca(a);
    manifold=c(:,1:dim);
    % null
    % tmp = randn(size(c(:,1:dim)));
    % [Q,~] = qr(tmp,0);
    % manifold = Q;

    atest = chdata_trauma(:,:,test_idx) - chdata_neutral(:,:,test_idx);
    atest = atest(:,idx_snippet)';
    atest = atest-mean(atest);

    ahat = atest - atest*(manifold*manifold');
    num = norm(atest,'fro')^2 -norm(ahat,'fro')^2;
    den = norm(atest,'fro')^2;
    vaf(i) = num/den;
end
figure;
boxplot(vaf)
ylim([0 1])





