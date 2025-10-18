%% GENERATING RAW DATA TO LOOK AT ERPs B/W VIEWING NEUTRAL AND STIM IMAGE VIDEOS

%% ANALYZING THE EEG DATA FROM hDEEG PTSD TRAUMA DATA

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

% go to the specific subject's folder
subj_name = ['Participant ' '24038'];
root_path = fullfile('/media/user/Data/Ana EEG/STAR/Phase 2',subj_name);
cd(root_path)

% starting eeglab
eeglab

% load the dataset
% EEG = pop_loadset('filename','24014_snippet-20250226T125238.set','filepath',...
%     'F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Participant 24014\');
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );


% call the functions to load the SAGA datafile
Poly5toEEGlab
%eeglab redraw

%% GETTING STIM AND FILTERING, referencing THE DATA
%load into matlab workspace
data=EEG.data;

% plot the stim trigger -> the channel differs by recording?
figure;
%stim = data(83,:);
stim = data(79,:);
stem(stim)

% remove unncesary triggers at the end
% [aa bb]=ginput;
% aa=round(aa);
% vline(aa,'r')
% stim(aa:end)=0;

% extract the stim onsets from the triggers
d = [0 diff(stim)];
d1 = (find(d==1)); % these are all the transitions from 0 to 1

% identifying stim onsets
stim_onsets=[];
done=false;k=1;
while ~done

    % now go forward till you have atleast 500 consecutive zeros.
    done1=false;
    init_k = k;
    while ~done1
        if d1(k+1)-d1(k) < 500
            k=k+1;
        else
            done1=true;
            stim_onsets = [stim_onsets d1(init_k)];
            k=k+1;
        end
    end

    % chk if cycled through all indices
    if k >= length(d1)
        done = true;
    end

end
stim_onsets = stim_onsets(1:2:end); % taking into account just the two pulse flashes

% now taking into account that the stim onsets have a 40ms delay due to
% photosensor settings
stim_onsets = stim_onsets-3; % look at the screenshots from Iryne
figure;stem(stim)
vline(stim_onsets,'r')
disp(['number of detected stim is ' num2str(length(stim_onsets))]);

% rewrite the stim channel
tmp = zeros(size(stim));
tmp(stim_onsets)=1;
figure;stem(tmp)
%data(83,:) = tmp;
data(79,:) = tmp;

%% Filter the data ...ZERO PADD TO ACCOUNT FOR FILTERING ARTIFACTS
data=double(data);
d1 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',30, ...
    'SampleRate',EEG.srate);
%fvtool(d1)
data(1:68,:)=filtfilt(d1,data(1:68,:)')';

% re-reference the data to M1/M2
ref = data([13 19],:);
ref = mean(ref,1);
data(1:68,:) = data(1:68,:) - ref;

% rewrite the eeglab data structure
EEG.data=data;
eeglab redraw
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
disp('LP Filtering Done')

%% EOG bad channel correction (if required)
data=EEG.data;
figure;plot(data(65,:))
ylim([-1000 1000])
figure;plot(data(66,:))
ylim([-1000 1000])

[aa bb]=ginput;
aa=round(aa);

data(65,aa(1):aa(2)) = median(data(65,:)) + sqrt(mad(data(65,:))) * randn(1,aa(2)-aa(1)+1);
data(65,aa(3):aa(4)) = median(data(65,:)) + sqrt(mad(data(65,:))) * randn(1,aa(4)-aa(3)+1);
hold on
plot(data(65,:))

figure;plot(data(66,:))
[aa bb]=ginput;
aa=round(aa);

data(66,aa(1):aa(2)) = median(data(66,:)) + sqrt(mad(data(66,:))) * randn(1,aa(2)-aa(1)+1);
%data(66,aa(3):aa(4)) = median(data(66,:)) + sqrt(mad(data(66,:))) * randn(1,aa(4)-aa(3)+1);
hold on
plot(data(66,:))




EEG.data=data;
eeglab redraw


%% EOG GUI BASED AC

% perform EOG artifact correct using channels 65 and 66
EEG = pop_scrls_regression( EEG, [65  66], 2, 0.9999, 0.01, 32, []);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off');

%% EOG NON GUI BASED AC

X=EEG.data;Xbkup=EEG.data;
indata = X(1:64,:);
opt.refdata = X(65:66,:);
opt.M = 2; % filter order
opt.lambda = 0.9999;
opt.sigmal = 0.01;
opt.prec = 32;

[Y,theta,Hh] = scrls_regression(indata, opt);

X(1:64,:)=Y;
EEG.data=X;
eeglab redraw
%[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

data_ac=EEG.data;
figure;plot(data_ac(34,:))
hold on
plot(Xbkup(34,:))
legend({'AC','Original'})

%% correct the stim signal

data_ac=EEG.data;
%stim_ac = data_ac(83,:);
stim_ac = data_ac(79,:);
figure;stem(stim_ac)
hold on
stim_ac(stim_ac>0.9)=1;
stim_ac(stim_ac<0.9)=0;
stem(stim_ac);
stem(tmp)
disp(sum(stim_ac==tmp)/length(tmp))
data_ac(83,:) = stim_ac;
if sum(stim_ac==tmp)/length(tmp) ==1

    EEG.data=data_ac;
    eeglab redraw
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off');

    % examine bad channels
    figure;
    for i=32:36
        plot(data_ac(i,:));
        title(num2str(i))
        ylim([-500 500])
        waitforbuttonpress;
    end

else

    disp('error')
end


%% remove bad channel and interpolate missing values, if necessary
%bad_ch=12;
bad_ch=[36]-2; % since we are removing M1 and M2
if ~isempty(bad_ch)
    EEG.data(bad_ch+2,:)=0;
    eeglab redraw
end
EEG = pop_select( EEG, 'rmchannel',{'M1','M2','SpO2','PLETH','HRate','Stat','TRIGGERS1',...
    'TRIGGERS2','TRIGGERS3','TRIGGERS4','TRIGGERS5','TRIGGERS6','TRIGGERS7','TRIGGERS8',...
    'TRIGGERS9','TRIGGERS10','TRIGGERS12','TRIGGERS13','TRIGGERS14','TRIGGERS15',...
    'TRIGGERS16','STATUS','Counter 2power24','BIP 01','BIP 02','BIP 03','BIP 04'});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off');
eeglab redraw

if ~isempty(bad_ch)

    EEG = pop_interp(EEG, [bad_ch], 'invdist'); % interploate the bad electrodes
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off');
    eeglab redraw
end

% epoching the data. Stim onset tell when the two pulses came on at the end
% of the fixation cross. So go back 1s and go forward 6s.
% A) Pre stim: fixation cross
% B) Stim: this corresponds to the two pulses, lasts for 350.35ms or 350ms
% C) Then the actual snippet video comes on, lasts for 4004ms
% so -1s to 6s captures all the events
EEG = pop_chanevent(EEG, 63,'edge','leading','edgelen',0,'delchan','off');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  'chan63'  }, [-1 6], 'newname', 'Data epoched', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off');
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [-1000 0] ,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off');
eeglab redraw
EEG = pop_select( EEG, 'rmchannel',{'TRIGGERS11'});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off');
eeglab redraw


% epoch the data based on the events -> use the EEGLAB GUI
%file->import event info -> from data channel
% tools -> extract epochs -> use stim channel
% prune channels -> edit -> select data -> remove m1,m2,stim and other
% channels


%% RUN ICA
% to remove artifacts


%% plot ERPs over occipital

EEG = pop_autorej(EEG, 'nogui','on','threshold',200,'eegplot','on'); % automatic bad trial rejection
bad_trials = find(EEG.reject.rejauto);% bad trials
close


data=EEG.data;
data(:,:,bad_trials) = NaN;
% bad_ep =[10 22];
% data =cat(3,data(:,:,1:9), NaN(62,7000,1), data(:,:,10:20),NaN(62,7000,1),...
%     data(:,:,21:end));

chdata=[];
chdata = squeeze(cat(3,chdata,data(28,:,:),data(29,:,:),data(30,:,:),data(27,:,:),...
    data(54,:,:),data(55,:,:)));
figure;
set(gcf,'Color','w')
plot(nanmean(chdata,2),'k','LineWidth',1)
hold on
mb = sort(bootstrp(1000,@nanmean,chdata'));
plot(mb(25,:),'Color',[.5 .5 .5 .5])
plot(mb(975,:),'Color',[.5 .5 .5 .5])
vline([1000 1300 1300+4004],'r')
hline(0,'r')
set(gca,'FontSize',14)
xlabel('Time (ms)')
ylabel('uV')
title('Occipital channels')

% plot ERPs after splitting into neutral and trauma snippets
neutral_epochs= [1 3 5 7 8 11 14 17 19 21 23 24 27 29 30 33 34 37 39 40];
idx_neutral = zeros(40,1);
idx_neutral(neutral_epochs)=1;
trauma_epochs = find(idx_neutral==0);
idx_trauma = zeros(40,1);
idx_trauma(trauma_epochs)=1;
idx_trauma = logical(idx_trauma);
idx_neutral = logical(idx_neutral);

% plot ERPs
data=EEG.data;
chdata1=[];
chdata1 = squeeze(cat(3,chdata1,data(28,:,idx_trauma),data(29,:,idx_trauma),data(30,:,idx_trauma),...
    data(27,:,idx_trauma),data(54,:,idx_trauma),data(55,:,idx_trauma)));
chdata2=[];
chdata2 = squeeze(cat(3,chdata2,data(28,:,idx_neutral),data(29,:,idx_neutral),data(30,:,idx_neutral),...
    data(27,:,idx_neutral),data(54,:,idx_neutral),data(55,:,idx_neutral)));
% ac2 = sum(abs(chdata2)>100);
% ac2 = find(ac2~=0);
% chdata2(:,ac2)=NaN;

figure;
set(gcf,'Color','w')

plot(nanmean(chdata1,2),'k','LineWidth',1)
hold on
mb = sort(bootstrp(1000,@nanmean,chdata1'));
plot(mb(25,:),'Color',[.5 .5 .5 .5])
plot(mb(975,:),'Color',[.5 .5 .5 .5])


plot(nanmean(chdata2,2),'b','LineWidth',1)
hold on
mb = sort(bootstrp(1000,@nanmean,chdata2'));
plot(mb(25,:),'Color',[.2 .5 .8 .5])
plot(mb(975,:),'Color',[.2 .5 .8 .5])


vline([1000 1300 1300+4004],'r')
hline(0,'r')
set(gca,'FontSize',14)
xlabel('Time (ms)')
ylabel('uV')
title('Occipital channels')
legend({'Trauma','','','Neutral','',''})
title('Subject 204022 Occipital Channels')


data=EEG.data;
chdata=squeeze(data(28,:,:));

%% ERPs over cz

% plot ERPs
data=EEG.data;
data(:,:,bad_trials) = NaN;
chdata1=[];
chdata1 = squeeze(cat(3,chdata1,data(15,:,idx_trauma),data(43,:,idx_trauma),data(44,:,idx_trauma),...
    data(47,:,idx_trauma),data(40,:,idx_trauma)));
chdata2=[];
chdata2 = squeeze(cat(3,chdata2,data(15,:,idx_neutral),data(43,:,idx_neutral),data(44,:,idx_neutral),...
    data(47,:,idx_neutral),data(40,:,idx_neutral)));
figure;
set(gcf,'Color','w')

plot(nanmean(chdata1,2),'k','LineWidth',1)
hold on
mb = sort(bootstrp(1000,@nanmean,chdata1'));
plot(mb(25,:),'Color',[.5 .5 .5 .5])
plot(mb(975,:),'Color',[.5 .5 .5 .5])


plot(nanmean(chdata2,2),'b','LineWidth',1)
hold on
mb = sort(bootstrp(1000,@nanmean,chdata2'));
plot(mb(25,:),'Color',[.2 .5 .8 .5])
plot(mb(975,:),'Color',[.2 .5 .8 .5])


vline([1000 1300 1300+4004],'r')
hline(0,'r')
set(gca,'FontSize',14)
xlabel('Time (ms)')
ylabel('uV')
title('Occipital channels')
legend({'Trauma','','','Neutral','',''})
title('Subject 204022 Central Channels')

%% ERPs over frontal

% plot ERPs
data=EEG.data;
data(:,:,bad_trials) = NaN;
chdata1=[];
chdata1 = squeeze(cat(3,chdata1,data(6,:,idx_trauma),data(36,:,idx_trauma),data(37,:,idx_trauma)));
chdata2=[];
chdata2 = squeeze(cat(3,chdata2,data(6,:,idx_neutral),data(36,:,idx_neutral),data(37,:,idx_neutral)));
figure;
set(gcf,'Color','w')

plot(nanmean(chdata1,2),'k','LineWidth',1)
hold on
mb = sort(bootstrp(1000,@nanmean,chdata1'));
plot(mb(25,:),'Color',[.5 .5 .5 .5])
plot(mb(975,:),'Color',[.5 .5 .5 .5])


plot(nanmean(chdata2,2),'b','LineWidth',1)
hold on
mb = sort(bootstrp(1000,@nanmean,chdata2'));
plot(mb(25,:),'Color',[.2 .5 .8 .5])
plot(mb(975,:),'Color',[.2 .5 .8 .5])


vline([1000 1300 1300+4004],'r')
hline(0,'r')
set(gca,'FontSize',14)
xlabel('Time (ms)')
ylabel('uV')
title('Occipital channels')
legend({'Trauma','','','Neutral','',''})
title('Subject 204022 Frontal Channels')


%
%
% data_preAAR = EEG.data;
%
% data1=EEG.data;
% for i=1:size(data1,3)
%     figure;plot(squeeze(data1(27,:,i)))
%     hold on
%     plot(squeeze(data_preAAR(27,:,i)))
%     plot(squeeze(data_preAAR(64,:,i)))
%     plot(squeeze(data_preAAR(63,:,i)))
%     legend('Artifact Corrected','Original','BP2','BP1')
%     title(num2str(i))
% end
%
%
% % using my RLS code
% data=EEG.data;
% x =  squeeze(data(2,:,4));
% d = squeeze(data(64,:,4));
% figure;
% plot(x)
% hold on
% plot(d)
%
%
% lam=0.9999;
% p=5;
% [w1,e,dhat,P,alp] = nik_rls(x,d,lam,p);
%
% plot(e);
% legend('Sig','BP1','AC')


%% ERP SNIPPET ANALYSES (MAIN) GROUP LEVEL SPATIOTEMPORAL CONTRASTS
%  looping over subjects loading the data and looking at grand averaged ERPs

%
% EEG = pop_autorej(EEG, 'nogui','on','threshold',300,'eegplot','on'); % automatic bad trial rejection
%  bad_trials = find(EEG.reject.rejauto);% bad trials

clc;clear
close all

% adding paths
%addpath(genpath('C:\Users\nikic\Documents\MATLAB\Ana EEG\hdEEG_Trauma_PTSD\SAGA_Matlab\SAGA_interface'))
%root_path = 'F:\DATA\EEG Data Anne Richards\STAR\STAR_Pilot_data-selected\Pilot18';
%root_path='F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Snippet_Datasets_Processed';
addpath(genpath('/home/user/Documents/Repositories/STAR_Study_EEG/SAGA_Matlab/SAGA_interface'))
addpath('/home/user/Documents/MATLAB/eeglab2023.1')
%addpath(genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114'))
addpath('/home/user/Documents/MATLAB')
root_path='/media/user/Data/Ana EEG/STAR/Phase 2/Snippet_Datasets_Processed';
cd(root_path)
eeglab

files=findfiles('.set',root_path)';

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
for i=1:length(files)


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


        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        eeglab redraw
        k=k+1;

    end
end


% spatiotemporal cluster correction
%addpath('C/home/user/Documents/Repositories/limo_tools/limo_tools')
addpath(genpath('/home/user/Documents/Repositories/limo_tools/'))
tfce_flag = true;
loop_iter=1000;
t_scores=[];tboot=zeros(62,7000,loop_iter);
p_scores=[];pboot=zeros(62,7000,loop_iter);

parfor ch=1:size(chdata_neutral,1)
    tmp_neutral = squeeze(chdata_neutral(ch,:,:));
    tmp_trauma = squeeze(chdata_trauma(ch,:,:));

    a = zscore(tmp_neutral')';
    I = abs(a)>4;
    tmp_neutral(I) = NaN;

    a = zscore(tmp_trauma')';
    I = abs(a)>4;
    tmp_trauma(I) = NaN;


    [h,p,ci,stats] = ttest(tmp_trauma',tmp_neutral');
    t = stats.tstat;
    if tfce_flag
        t(p>0.05)=0; % only if TFCE
    end
    t_scores(ch,:) = t;
    p_scores(ch,:) = p;


    % null hypothesis testing using permutation tests        
    for loop=1:loop_iter

        d = tmp_trauma - tmp_neutral;
        flip = sign(randn(size(d)));
        d = d.*flip;
        [h0 p0 ci stats0] = ttest(d');
        t0=stats0.tstat;
        if tfce_flag
            t0(p0>0.05)=0; % only if TFCE
        end
        tboot(ch,:,loop)=t0;
        pboot(ch,:,loop)=p0;
    end
end

% get neighborhood distance matrix
neighb = limo_neighbourdist(EEG, 0.40);


if tfce_flag
    E=1;H=2;dh=0.2;
    [tfce_score,~] = limo_tfce(2,t_scores,neighb,1,E,H,dh);
    [tfce_score_boot,~] = limo_tfce(2,tboot,neighb,1,E,H,dh);

    % get the null distribution of tfce score from each bootstrap
    tfce_boot=[];
    for loop=1:size(tfce_score_boot,3)
        a=squeeze(tfce_score_boot(:,:,loop));
        tfce_boot(loop) = max(a(:));
    end

    % threshold the true tfce scores with the null distribution
    tfce_boot=sort(tfce_boot);tfce_score1=tfce_score;
    thresh = tfce_boot(round(0.95*length(tfce_boot)));
    tfce_score1(tfce_score1<thresh)=0;
    figure;
    subplot(3,1,1)
    tt=linspace(-1,6,size(t_scores,2));
    imagesc(tt,1:62,t_scores);
    title('Uncorrected for multiple comparisons')
    subplot(3,1,2)
    imagesc(tt,1:62,tfce_score1)
    ylabel('Channels')
    xlabel('Time')
    title('Spatiotemp. multiple comparison corrected')

    % plot the significant channels
    a=tfce_score1;
    aa=sum(a(:,1000:6000),2);
    aa=aa./max(aa);% normalize
    subplot(3,1,3);    
    topoplot((aa),EEG.chanlocs,'maplimits', [-1 1])
    axis tight    
    title('Relative Channel importance')    
    axis off
    set(gcf,'Color','w')

else
    % 2D spatiotemporal cluster correction
    LIMO.data.chanlocs=[];
    LIMO.data.neighbouring_matrix=neighb;
    [mask,cluster_p,max_th] = ...
        limo_clustering((t_scores.^2),p_scores,...
        (tboot.^2),pboot,LIMO,2,0.05,0);
    figure;subplot(3,1,1)
    tt=linspace(-1,6,size(t_scores,2));
    imagesc(tt,1:62,t_scores);
    title('Uncorrected for multiple comparisons')
    ylabel('Channels')
    xlabel('Time')
    subplot(3,1,2)
    t_scores1=t_scores;
    t_scores1(mask==0)=0;
    imagesc(tt,1:62,t_scores1);
    title('Corrected for multiple comparisons')
    ylabel('Channels')
    xlabel('Time')

    a=mask;
    aa=sum(a(:,3000:6000),2);
    sig_ch_idx = find(aa>0);
    sig_ch = zeros(numel(chMap),1);
    sig_ch(sig_ch_idx)=1;
    sig_ch_all(i,:) = sig_ch;
    subplot(3,1,3);
    imagesc(sig_ch(chMap))
    title('Sig channels 0 to 3s')
    sgtitle(ImaginedMvmt{i})
    axis off
    set(gcf,'Color','w')
    colorbar
end

% make a movie
v=VideoWriter('ERP_EvsN_TFCE');
t_scores1=tfce_score1;
t_scores1 = t_scores1./max(t_scores1(:));
v.FrameRate=16;
open(v)
figure;
M={};
vmin = min(t_scores1(:));
vmax = max(t_scores1(:));
for i=750:2:size(t_scores1,2)
    disp(i)
    tmp = t_scores1(:,i);
    % topoplot(tmp,EEG.chanlocs, 'maplimits', [-6 6],...
    %     'numcontour', 0);    
    topoplot(tmp,EEG.chanlocs, 'maplimits', [-1 1],...
        'numcontour', 0);  
    colorbar
    axis tight
    sgtitle([num2str(tt(i)) ' Sec'])
    frame=getframe(gcf);
    img = frame2im(frame); 
    %M{i} = (img);
    writeVideo(v,img)
    clf
end
close(v)

% plot ERPs
channels={EEG.chanlocs.labels};
ch=[6 15 27];
yoff=[-9 -9.5 -3];
figure;
for i=1:length(ch)
    subplot(3,1,i)
    hold on
    plot(tt,squeeze(nanmean(chdata_neutral(ch(i),:,:),3)),'LineWidth',2);
    plot(tt,squeeze(nanmean(chdata_trauma(ch(i),:,:),3)),'LineWidth',2)
    vline(tt([1000 1300 1300+4004]),'k')
    hline(0,'r')
    set(gca,'FontSize',14)
    title(channels{ch(i)});
    sig_time = t_scores1(ch(i),:);    
    idx = find(sig_time~=0);
    % idx = tt(idx);
    % sigMask= idx(idx>=0);
    sigDiff = diff([0 idx]);
    segStart = find(sigDiff >1);
    %vline(idx(segStart),'k')
    segEnd = [segStart(2:end)-1];
    segStart = idx(segStart);
    segEnd = [idx(segEnd) idx(end)];
    

    for ii = 1:length(segStart)
        if tt(segStart(ii)) >0
            x1 = tt(segStart(ii));
            x2 = tt(segEnd(ii));
            plot([x1 x2], [yoff(i) yoff(i)], 'g', 'LineWidth', 5);
        end
    end

    if i==2
        ylim([-11 8])
    end
    %axis tight
end

%
% figure;
% set(gcf,'Color','w')
% chdata_trauma = chdata_trauma';
% chdata_neutral = chdata_neutral';
%
% plot(nanmean(chdata_trauma,2),'k','LineWidth',1)
% hold on
% mb = sort(bootstrp(1000,@nanmean,chdata_trauma'));
% plot(mb(25,:),'Color',[.5 .5 .5 .5])
% plot(mb(975,:),'Color',[.5 .5 .5 .5])
%
%
% plot(nanmean(chdata_neutral,2),'b','LineWidth',1)
% hold on
% mb = sort(bootstrp(1000,@nanmean,chdata_neutral'));
% plot(mb(25,:),'Color',[.2 .5 .8 .5])
% plot(mb(975,:),'Color',[.2 .5 .8 .5])
%
%
% vline([1000 1300 1300+4004],'r')
% hline(0,'r')
% set(gca,'FontSize',14)
% xlabel('Time (ms)')
% ylabel('uV')
% title('Occipital channels')
% legend({'Trauma','','','Neutral','',''})
% title('Frontal Channels 13 Subjects')
%
% % running t-tests
% res=[];tval=[];
% for i=1:size(chdata_trauma,1)
%     a = chdata_trauma(i,:);
%     b = chdata_neutral(i,:);
%     [h p tb st]=ttest(a,b);
%     tval(i) = st.tstat;
%     if p<0.05
%         res(i)=1;
%     else
%         res(i)=0;
%     end
% end
% figure;stem(res)
%
% x = tinv(0.95,length(files)-1);
% figure;plot(tval,'LineWidth',1)
% axis tight
% vline([1000 1300 1300+4004],'r')
% hline(0,'r')
% set(gca,'FontSize',14)
% xlabel('Time (ms)')
% ylabel('uV')
% title('Occipital channels')
% title('T-value Frontal channels')
% set(gcf,'Color','w')
% hline(x,'--r')



%%

cd('F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Videos')
vr = VideoReader('4S_NOV21-24.mp4');

frame = readFrame(vr);  % Read the next frame
frames = read(vr, [1 60*7]); % Read frames 10 to 20

fr=vr.frameRate;
figure;
for i=100:size(frames,4)
    imshow(frames(:,:,:,i));
    t = i*fr;
    title(num2str(i))
    waitforbuttonpress
end

%%
% going through a loop and now detecting the average intensity around the
% area where the white marker square appears

ylimits= [983:998];
xlimits= [1888:1904];
pixel_val=[];
for i=1:vr.NumFrames
    disp(i)
    tmp = read(vr,i);
    tmp=tmp(ylimits,xlimits,:);
    pixel_val(i) = mean(tmp(:));
end

figure;stem(pixel_val)
pv_bkup = pixel_val;

% now going through and identifying the frame duration of each event
pixel_val(pixel_val<100)=0;
pixel_val(pixel_val>100)=1;

frame_rate = vr.FrameRate;
dt = 1/frame_rate;
frame_counter=[];
done=false;
frameid_past = pixel_val(1);
counter=1;
i=2;
while ~done
    frameid = pixel_val(i);
    if frameid_past == frameid
        counter=counter+1;
        frameid_past = frameid;
        i=i+1;

    else
        tmp = [frameid_past counter counter*dt];
        frame_counter =[frame_counter;tmp];
        frameid_past = frameid;
        i=i+1;
        counter=1;
    end

    if i==length(pixel_val)
        done=true;
    end
end

stim_pulse = frame_counter(:,1);
num_frames = frame_counter(:,2);
duration = frame_counter(:,3);

frame_counter =table(stim_pulse, num_frames, duration);




%% LOADING SAVING NAP DATA TO EDF FILE FORMAT (PROCESSED)

clc;clear;
close all


% adding paths
root_path='F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Participant 24015';
cd(root_path)
addpath(genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114'))


% Define Poly5 file path
filename = '24015_nap-20250312T132123.DATA.poly5'; % Change this to your actual Poly5 file

% Read the Poly5 file
cfg = [];
cfg.dataset = filename;
data = ft_preprocessing(cfg);

% filter the data
d1 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.3,'HalfPowerFrequency2',35, ...
    'SampleRate',1000);

d2 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.15,'HalfPowerFrequency2',15, ...
    'SampleRate',1000);

% export: left and right eog, F3, F4, C3, C4, O1, O2
% common reference to Cz
% eog: .15/15 hz
% others:  .3-35.

% Extract signal data
signal = data.trial{1}; % The actual EEG/EMG/etc. data
time = data.time{1};    % Corresponding time vector
fs = data.fsample;      % Sampling rate
nChannels = size(signal,1); % Number of channels

% Define channel labels (FieldTrip stores them in .label)
channelNames = data.label;

% get channels
% eog1 eog2 f3 f4 c3 c4 01 02 cz (ref)
ch = [65 66 5 7 15 17 30 32 16 ];
data1 = signal(ch,:);

% filter
data1(1:2,:) = filtfilt(d2,data1(1:2,:)')';
data1(3:end,:) = filtfilt(d1,data1(3:end,:)')';

% re-reference to Cz
data1 = data1 - data1(end,:);

% remove Cz from dataset
data1 = data1(1:end-1,:);
channelNames = channelNames(ch(1:end-1));
signal = data1;

% Convert to EDF format using ft_write_data
cfg = [];
cfg.filename = '24015_Nap_Proc.edf';  % Output EDF filename
cfg.dataformat = 'edf';       % Specify EDF format
cfg.Fs = fs;             % Sampling rate
cfg.label = channelNames;      % Channel names
%cfg.refchan = 'all';          % Use all channels
ft_write_data(cfg.filename, signal, 'header',cfg,'dataformat', 'edf');



%% LOADING SAVING NAP DATA TO EDF FILE FORMAT (UNPROCESSED)

clc;clear;
close all


% adding paths
root_path='F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Participant 24015';
cd(root_path)
addpath(genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114'))


% Define Poly5 file path
filename = '24015_nap-20250312T132123.DATA.poly5'; % Change this to your actual Poly5 file

% Read the Poly5 file
cfg = [];
cfg.dataset = filename;
data = ft_preprocessing(cfg);

% filter the data
% d1 = designfilt('bandpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',0.3,'HalfPowerFrequency2',35, ...
%     'SampleRate',1000);

d1 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',35, ...
    'SampleRate',1000);

% export: left and right eog, F3, F4, C3, C4, O1, O2
% common reference to Cz
% eog: .15/15 hz
% others:  .3-35.

% Extract signal data
signal = data.trial{1}; % The actual EEG/EMG/etc. data
time = data.time{1};    % Corresponding time vector
fs = data.fsample;      % Sampling rate
nChannels = size(signal,1); % Number of channels

% Define channel labels (FieldTrip stores them in .label)
channelNames = data.label;

% get channels
% eog1 eog2 f3 f4 c3 c4 01 02 cz (ref)
% ch = [65 66 5 7 15 17 30 32 16 ];
% data1 = signal(ch,:);

% filter
rmpath((genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114')))
rehash toolbox
data1 = filtfilt(d1,signal')';
%data1(3:end,:) = filtfilt(d1,data1(3:end,:)')';

% reference to M1 and M2, channels 1 to 68
ref = data1([13 19],:);
ref = mean(ref,1);
data1(1:68,:) = data1(1:68,:) - ref;




% Convert to EDF format using ft_write_data
addpath(genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114'))
signal=data1;
cfg = [];
cfg.filename = '24015_Nap_Proc_AllChan.edf';  % Output EDF filename
cfg.dataformat = 'edf';       % Specify EDF format
cfg.Fs = fs;             % Sampling rate
cfg.label = channelNames;      % Channel names
%cfg.refchan = 'all';          % Use all channels
ft_write_data(cfg.filename, signal, 'header',cfg,'dataformat', 'edf');





%% ANALYZING SLEEP DATA TO SEE IF LENGTH MATCHES UP

% the first 30s of sleep epoch corresponds to the very first 30k samples of
% the sleep data itself.
% the data has been contralateral referenced to the mastoids
% band pass b/w .3 to 35z
% the last <30s samples is discarded.


clc;clear

filepath='F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Participant 24020';
cd(filepath)
addpath('C:\Users\nikic\Documents\MATLAB\eeglab2023.1')
%addpath(genpath('C:\Users\nikic\Documents\MATLAB\fieldtrip-20250114'))
addpath('C:\Users\nikic\Documents\MATLAB')
eeglab

% find the number of hours and minutes
total_samples = size(EEG.data,2)/1000; %seconds
num_hours = floor(total_samples/3600)
num_sec =  rem(total_samples,3600)
num_30s_seg = total_samples/30

