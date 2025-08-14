%% ANALYZING THE EEG DATA FROM PILOT OF PTSD TRAUMA DATA

clc;clear
close all

% adding paths
addpath(genpath('C:\Users\nikic\Documents\MATLAB\Ana EEG\hdEEG_Trauma_PTSD\SAGA_Matlab\SAGA_interface'))
root_path = 'F:\DATA\EEG Data Anne Richards\STAR\STAR_Pilot_data-selected\Pilot18';
cd(root_path)
addpath('C:\Users\nikic\Documents\MATLAB\eeglab2023.1')

% starting eeglab
eeglab

% call the functions to load the SAGA datafile 
%Poly5toEEGlab
eeglab redraw

%load into matlab workspace
data=EEG.data;

% plot the stim trigger -> the channel differs by recording? 
figure;
stim = data(83,:);
plot(stim)

% extract the stim onsets from the triggers -> gets the 2 stim triggers
d = diff(stim);
d1 = (find(d==-1));
dd = [0.5e3 diff(d1)];
idx = dd<500;
stim_onsets = d1(idx);
vline(stim_onsets,'r')

% extract 1 stim triggers

% rewrite the stim channel 
tmp = zeros(size(stim));
tmp(stim_onsets)=1;
figure;stem(tmp)
data(83,:) = tmp;

% filter the data ...ZERO PADD TO ACCOUNT FOR FILTERING ARTIFACTS
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

% epoch the data based on the events -> use the EEGLAB GUI
%file->import event info -> from data channel
% tools -> extract epochs -> use stim channel
% prune channels -> edit -> select data -> remove m1,m2,stim and other
% channels
EEG = pop_chanevent(EEG, 83,'edge','leading','edgelen',0,'delchan','off');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  'chan83'  }, [-0.5         3.5], 'newname', 'Continuous Data TMSi Name: DATA (start 6/27/2024 4:17 PM) epochs', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [-500 0] ,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
EEG = pop_select( EEG, 'rmchannel',{'M1','M2','SpO2','PLETH','HRate','Stat','TRIGGERS1','TRIGGERS2','TRIGGERS3','TRIGGERS4','TRIGGERS5','TRIGGERS6','TRIGGERS7','TRIGGERS8','TRIGGERS9','TRIGGERS10','TRIGGERS11','TRIGGERS12','TRIGGERS13','TRIGGERS14','TRIGGERS15','TRIGGERS16','STATUS','Counter 2power24'});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
eeglab redraw

% RUN ICA
% to remove artifacts


% plot ERPs
data=EEG.data;
chdata=[];
chdata = squeeze(cat(3,chdata,data(28,:,:),data(29,:,:),data(30,:,:),data(27,:,:),...
    data(54,:,:),data(55,:,:)));
figure;
set(gcf,'Color','w')
plot(mean(chdata,2),'k','LineWidth',1)
hold on
mb = sort(bootstrp(1000,@mean,chdata'));
plot(mb(25,:),'Color',[.5 .5 .5 .5])
plot(mb(975,:),'Color',[.5 .5 .5 .5])
vline(500,'r')
hline(0,'r')
set(gca,'FontSize',14)
xlabel('Time (ms)')
ylabel('uV')
title('Occipital channels')

data_preAAR = EEG.data;

data1=EEG.data;
for i=1:size(data1,3)
    figure;plot(squeeze(data1(27,:,i)))
    hold on
    plot(squeeze(data_preAAR(27,:,i)))
    plot(squeeze(data_preAAR(64,:,i)))
    plot(squeeze(data_preAAR(63,:,i)))
    legend('Artifact Corrected','Original','BP2','BP1')
    title(num2str(i))
end


% using my RLS code
data=EEG.data;
x =  squeeze(data(2,:,4));
d = squeeze(data(64,:,4));
figure;
plot(x)
hold on
plot(d)


lam=0.9999;
p=5;
[w1,e,dhat,P,alp] = nik_rls(x,d,lam,p);

plot(e);
legend('Sig','BP1','AC')


%% EXTRACTING TIMING INFORMATION FROM THE TRIGGERS


clc;clear
close all

% adding paths
addpath(genpath('C:\Users\nikic\Documents\MATLAB\Ana EEG\hdEEG_Trauma_PTSD\SAGA_Matlab\SAGA_interface'))
%root_path = 'F:\DATA\EEG Data Anne Richards\STAR\STAR_Pilot_data-selected\Pilot18';
root_path = 'F:\DATA\EEG Data Anne Richards\STAR\STAR_Pilot_data-selected\Marker_Snippet_Testing_12_11_24';
cd(root_path)
addpath('C:\Users\nikic\Documents\MATLAB\eeglab2023.1')

% starting eeglab
eeglab

% call the functions to load the SAGA datafile 
%Poly5toEEGlab
eeglab redraw

%load into matlab workspace
data=EEG.data;

% plot the stim trigger -> the channel differs by recording? 
figure;hold on
stim = data(83,:);
plot(stim,'Color',[.5 .5 .5 .5])
ylim([-.5 1])

% extract the stim onsets from the triggers -> gets the 2 stim triggers
d = diff(stim);
d1 = (find(d==-1));
dd = [0.5e3 diff(d1)];
idx = dd<500;
stim_onsets = d1(idx);
vline(stim_onsets,'r')

% extract 1 stim triggers
d = diff(stim);
d1 = (find(d==-1)); 
single_stim_idx=[];
for i=2:length(d1)-1
    if abs(d1(i-1) - d1(i)) > 500 && abs(d1(i+1) - d1(i)) > 500
        single_stim_idx = [single_stim_idx i];
    end
end
stim_onsets1 = d1(single_stim_idx);

vline(stim_onsets1,'g')

