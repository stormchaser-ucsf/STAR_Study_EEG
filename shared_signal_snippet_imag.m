

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


% experimental 
idx_snippet = 2e3:5e3;
idx_imag = 8.5e3:11.5e3;

% control
% idx_snippet = [1:1e3  6e3:7e3];
% idx_imag = [1:1e3  13e3:14e3];

a = chdata_trauma - chdata_neutral;
b = chdata_trauma_imag - chdata_neutral_imag;


a = squeeze(mean(a,3));
b = squeeze(mean(b,3));


a=a(:,idx_snippet)';
b=b(:,idx_imag)';
a=double(a);
b=double(b);
%a=a-mean(a);
%b=b-mean(b);

[c1,s1,l1] = pca(a);
[c2,s2,l2] = pca(b);

% Q1: what is the shared dimensionality? 
figure;stem(cumsum(l1)./sum(l1))
hold on
stem(cumsum(l2)./sum(l2))
legend({'Snippet','Imag'})
hline(0.95)
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
%parpool('threads')
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

%Q4: comparison with baseline manifold structure?
angles_baseline=[];
% idx_snippet_baseline = [1:1e3  6e3:7e3];
% idx_imag_baseline = [1:1e3  13e3:14e3];
idx_snippet_baseline = [1:1e3 ];
idx_imag_baseline = [1:1e3 ];
pcap_baseline=[];
dim=11;
pcap_sweep_boot=[];
parfor i=1:1000
    idx = randi(size(chdata_neutral,3),size(chdata_neutral,3),1);
    a = double(chdata_trauma - chdata_neutral);
    a1 = a(:,:,idx);    
    a1 = squeeze(mean(a1,3));
    a1 = a1(:,idx_snippet_baseline)';
    %a1 = a1-mean(a1,2);

    idx = randi(size(chdata_neutral_imag,3),size(chdata_neutral_imag,3),1);
    b = double(chdata_trauma_imag - chdata_neutral_imag);
    b1 = b(:,:,idx);    
    b1 = squeeze(mean(b1,3));
    b1 = b1(:,idx_imag_baseline)';
    %b1 = b1-mean(b1,2);


    [c11,s11,l11] = pca(a1,'Centered','on');
    [c22,~,~] = pca(b1,'Centered','on');

    angles_baseline(i,:) = principal_angles(c11(:,1:dim),c22(:,1:dim));
    
    num = (c22(:,1:dim)*c22(:,1:dim)')*(c11(:,1:dim)*c11(:,1:dim)')*...
        (c22(:,1:dim)*c22(:,1:dim)');
    den = c11(:,1:dim)*c11(:,1:dim)';
    pcap_baseline(i,:) = trace(num)/trace(den);

    for j=5:20
        num = (c22(:,1:j)*c22(:,1:j)')*(c11(:,1:j)*c11(:,1:j)')*...
            (c22(:,1:j)*c22(:,1:j)');
        den = c11(:,1:j)*c11(:,1:j)';
        pcap_sweep_boot(i,j) = trace(num)/trace(den);
    end

end
pcap_sweep_boot=pcap_sweep_boot(:,5:end);

a = chdata_trauma - chdata_neutral;
a = squeeze(mean(a,3));
b = chdata_trauma_imag - chdata_neutral_imag;
b = squeeze(mean(b,3));

% experimental 
idx_snippet = 2e3:5e3;
idx_imag = 8.5e3:11.5e3;

a=a(:,idx_snippet)';
b=b(:,idx_imag)';
a=double(a);
b=double(b);

% a=a-mean(a,2);
% b=b-mean(b,2);

[c1,s1,l1] = pca(a,'Centered','on');
[c2,s2,l2] = pca(b,'Centered','on');


angles_main = principal_angles(c1(:,1:dim),c2(:,1:dim));
angles_baseline= sort(angles_baseline);
figure;
plot(angles_main)
hold on
plot(angles_baseline(25,:),'Color',[.5 .5 .5 .5])
plot(angles_baseline(975,:),'Color',[.5 .5 .5 .5])
legend({'Task specific angles','Angles in baseline period'})

num = (c2(:,1:dim)*c2(:,1:dim)')*(c1(:,1:dim)*c1(:,1:dim)')*...
    (c2(:,1:dim)*c2(:,1:dim)');
den = c1(:,1:dim)*c1(:,1:dim)';
pcap= trace(num)/trace(den);
figure;hist(pcap_baseline)
vline(pcap)
title([' p = ' num2str(1-sum(pcap>=pcap_baseline)/length(pcap_baseline))])

pcap_sweep_main=[];
for j=5:20
    num = (c2(:,1:j)*c2(:,1:j)')*(c1(:,1:j)*c1(:,1:j)')*...
        (c2(:,1:j)*c2(:,1:j)');
    den = c1(:,1:j)*c1(:,1:j)';
    pcap_sweep_main(j-4) = trace(num)/trace(den);
end
figure;plot(5:20,pcap_sweep_main)
xticks(5:20)
pcap_sweep_boot = sort(pcap_sweep_boot);
hold on
plot(5:20,pcap_sweep_boot(25,:),'Color',[.5 .5 .5 .5])
plot(5:20,pcap_sweep_boot(975,:),'Color',[.5 .5 .5 .5])
legend({'Task-specific','95% interval'})
xlabel('Dimensions')
ylabel('Percent variance captured between manifolds')
plot_beautify

% 
% % DOING ABOVE BUT STACKING ACROSS SUBEJCTS
% % SO PCA NOT ON SUBJECT AVERAGED BUT ACROSS ALL SUBJECTS STACKED
% angles_baseline=[];
% idx_snippet_baseline = [1:1e3 ];
% idx_imag_baseline = [1:1e3 ];
% pcap_baseline=[];
% dim=8;
% parfor i=1:1000
%     idx = randi(size(chdata_neutral,3),size(chdata_neutral,3),1);
%     a = double(chdata_trauma - chdata_neutral);
%     a1 = a(:,idx_snippet_baseline,idx);    
%     a1 = reshape(a1,size(a1,1),[])';
% 
%     idx = randi(size(chdata_neutral_imag,3),size(chdata_neutral_imag,3),1);
%     b = double(chdata_trauma_imag - chdata_neutral_imag);
%     b1 = b(:,idx_imag_baseline,idx);    
%     b1 = reshape(b1,size(b1,1),[])';
% 
% 
%     [c11,s11,l11] = pca(a1);
%     [c22,~,~] = pca(b1);
% 
%     angles_baseline(i,:) = principal_angles(c11(:,1:dim),c22(:,1:dim));
% 
%     num = (c22(:,1:dim)*c22(:,1:dim)')*(c11(:,1:dim)*c11(:,1:dim)')*...
%         (c22(:,1:dim)*c22(:,1:dim)');
%     den = c11(:,1:dim)*c11(:,1:dim)';
%     pcap_baseline(i,:) = trace(num)/trace(den);
% 
% end
% 
% % in actual data
% % experimental 
% idx_snippet = 2e3:5e3;
% idx_imag = 8.5e3:11.5e3;
% 
% a = chdata_trauma - chdata_neutral;
% a = a(:,idx_snippet,:);
% a= reshape(a,size(a,1),[]);
% 
% b = chdata_trauma_imag - chdata_neutral_imag;
% b = b(:,idx_imag,:);
% b= reshape(b,size(b,1),[]);
% 
% a=double(a)';
% b=double(b)';
% 
% [c1,s1,l1] = pca(a);
% [c2,s2,l2] = pca(b);
% 
% 
% angles_main = principal_angles(c1(:,1:dim),c2(:,1:dim));
% angles_baseline= sort(angles_baseline);
% figure;
% plot(angles_main)
% hold on
% plot(angles_baseline(25,:),'Color',[.5 .5 .5 .5])
% plot(angles_baseline(975,:),'Color',[.5 .5 .5 .5])
% legend({'Task specific angles','Angles in baseline period'})
% 
% num = (c2(:,1:dim)*c2(:,1:dim)')*(c1(:,1:dim)*c1(:,1:dim)')*...
%     (c2(:,1:dim)*c2(:,1:dim)');
% den = c1(:,1:dim)*c1(:,1:dim)';
% pcap= trace(num)/trace(den);
% figure;hist(pcap_baseline)
% vline(pcap)
% title([' p = ' num2str(1-sum(pcap>=pcap_baseline)/length(pcap_baseline))])
% 
