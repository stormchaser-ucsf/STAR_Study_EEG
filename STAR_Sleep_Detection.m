
%% PRELIMS


clc;clear
%root_folder = 'F:\DATA\EEG Data Anne Richards\STAR\Phase 2\';

root_folder = '/media/user/Data/Ana EEG/STAR/Phase 2/';
cd(root_folder)

%cd('C:\Users\nikic\Documents\GitHub\STAR_Study_EEG')
% addpath('C:\Users\nikic\Documents\MATLAB\')
% addpath('C:\Users\nikic\Documents\MATLAB\eeglab2023.1')
% addpath(genpath('C:\Users\nikic\Documents\GitHub\STAR_Study_EEG\SAGA_Matlab'))
% addpath('C:\Users\nikic\Documents\GitHub\STAR_Study_EEG\helpers')

addpath('/home/user/Documents/MATLAB/')
addpath('/home/user/Documents/MATLAB/eeglab2023.1')
addpath(genpath('/home/user/Documents/Repositories/STAR_Study_EEG'))


% starting eeglab
eeglab


%% LOAD THE PLOY5 AND EDF FILE

% now load the sleep data for this particular subject and save it as .set
%Poly5toEEGlab

data=EEG.data;
data=double(data);

disp('data loading done')

%% SLEEP STAGING AND FILTERS

% low pass filters
lpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',30, ...
    'SampleRate',1e3);

bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',65,'HalfPowerFrequency2',100, ...
    'SampleRate',1e3);

spFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',8,'HalfPowerFrequency2',16, ...
    'SampleRate',1e3);

spFilt1 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',9,'HalfPowerFrequency2',13, ...
    'SampleRate',1e3);

spFilt2 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',16, ...
    'SampleRate',1e3);

soFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.16,'HalfPowerFrequency2',1.25, ...
    'SampleRate',1e3);

deltaFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',4, ...
    'SampleRate',1e3);


% ANALYSIS FOR ONE SUBJECT , EXAMPLE
% detect SOs and Spindles
% sort by sleep stage
% detect the peak spindle frequency by investigating the power spectrum
% peak

filename=  '24011_SleepStages';
%filepath = 'F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Participant 24012\Sleep Stage Report';
filepath = '/media/user/Data/Ana EEG/STAR/Phase 2/Participant 24011/Sleep Stage Report/';
filename = fullfile(filepath,filename);
sleep_profile = table2array(importfile_sleep(filename,'Data',[3, Inf]));
figure;plot(sleep_profile)

% convert to full 1khz 30s arrays
Fs=1e3;
sleep_staging = zeros(Fs * 30 * length(sleep_profile),1);
k=1;
for i=1:(30*Fs):length(sleep_staging)
    t = i : i+30*Fs-1;
    sleep_staging(t) = sleep_profile(k);
    k=k+1;
end
figure;plot(sleep_staging)

if sleep_staging(end)<size(EEG.data,2)
    sleep_staging(end+1:size(EEG.data,2)) = sleep_staging(end);
end

disp('Filter design and sleep staging loading done')

%% BASIC PRE PROCESSING

% band pass filter the data
d1 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',30, ...
    'SampleRate',Fs);
data = filtfilt(d1,data')';

% re-reference the data to M1/M2
ref = data([13 19],:);
ref = mean(ref,1);
data(1:68,:) = data(1:68,:) - ref;

% restrict to the 62 neural channels, 2 reference and 2 EOG
data = data(1:66,:);


disp('Re referenced and hpf done')

%% BAD CHANNEL IDENTIFICATION

figure;
for i=1:size(data,1)
    plot(data(i,:))
    title(['ch ' num2str(i)])
    waitforbuttonpress
end


%% EYE BLINK AND ARTIFACT REMOVAL

opt.refdata = data(65:66,:);
opt.M=2;
opt.lambda = 0.9999;
opt.sigma = 0.01;
opt.prec=32;
[Y,H,Hh] = scrls_regression( data, opt);
data=Y;
clear Y H Hh

% remove the EOG channels
neural_ch = [1:64];
data = data(neural_ch,:);

disp('Eye blink and artifact correction done')

%% GET PREFERRED SPINDLE FREQUENCY PER CHANNEL

% get center frequency of spindle per channel using fft
topo_spindle_freq=[];
win_length=2048;
parfor ch=1:64
    disp(['Processing Channel ' num2str(ch)]);
    chdata = data(ch,:);
    seg = 1:win_length:length(chdata);
    ch_spindle_freq=NaN(1,length(seg));
    for i=1:length(seg)
        %disp(i/length(seg)*100)
        if i==length(seg)
            tmp = chdata(seg(i):end);
            sleep_seg = sleep_staging(seg(i):end);
        else
            tmp = chdata(seg(i):seg(i+1)-1);
            sleep_seg = sleep_staging(seg(i):seg(i+1)-1);
        end
        %[psdx,ffreq,phasex]=fft_compute(tmp,Fs,1);
        if length(tmp) >= 1024
            [Pxx,F]=pwelch(tmp,1024,256,1024,Fs);


            % remove the 1/f component and see which is freq with largest power
            idx = logical((F>0) .* (F<=25));
            F1=F(idx);
            F1=log2(F1);
            power_spect = Pxx(idx);
            power_spect = log2(power_spect);
            %[bhat p wh se ci t_stat]=robust_fit(F1,power_spect,1);
            warning('off','all');  % Turns off all warnings
            %opts = statset('MaxIter',1000);  % set max iterations to 1000
            tb=fitlm(F1,power_spect,'RobustOpts','on');
            warning('on','all');  % Turns off all warnings
            bhat = tb.Coefficients.Estimate;
            x = [ones(length(F1),1) F1];
            yhat = x*bhat;

            %plot
            % figure;
            % plot(F1,power_spect,'LineWidth',1);
            % hold on
            % plot(F1,yhat,'LineWidth',1);

            % remove the 1/f from the power spectrum
            power_spect = zscore(power_spect - yhat);

            % find which frequency has the most power
            ff = 2.^F1;
            idx = logical((ff>=7.9) .* (ff<=16.1));
            ff = ff(idx);
            power_spindle = power_spect(idx);
            [aa bb]=max(power_spindle);
            pref_freq = ff(bb);

            % store if in N2 sleep
            if sum(sleep_seg==2) > 0
                ch_spindle_freq(i) = pref_freq;
            end
        end

    end

    % figure;hist(ch_spindle_freq)
    % figure;ksdensity(ch_spindle_freq)
    % nanmean(ch_spindle_freq)

    topo_spindle_freq(ch) = nanmedian(ch_spindle_freq);

end

figure;plot(topo_spindle_freq)

%% PERFORM SPINDLE AND SO ANALYSES


slow_spindle_ch = [1:12 33:43 59 60];
fast_spindle_ch=[14:18 20:32 44:58 61:64];
ref_ch = [13 19];
slow_fast_idx=zeros(64,1);
slow_fast_idx(fast_spindle_ch)=1;

% spindle analyses in all channels
I = sleep_staging>0;
grid_sp = detect_spindles(data,I,soFilt,spFilt1,spFilt2,sleep_staging,Fs,slow_fast_idx);

% plotting
x=grid_sp.ch24;
t=linspace(-2,2,4001);
figure;
plot(t,nanmean(x.sp_epochs,1))
plot_beautify

% storing file name
grid_sp.filename = filename;
grid_sp.sleep_staging = sleep_staging;

%store data
grid_sp_exp{subj} = grid_sp;


% Examining SO-spindle nesting in stage 3 sleep
grid_so = SO_analyses(data,I,soFilt,spFilt1,spFilt2,sleep_staging,slow_fast_idx,Fs);
grid_so.filename = filename;
grid_so.sleep_staging = sleep_staging;

% plotting
x=grid_so.ch24;
t=linspace(-2.5,2.5,5001);
figure;
plot(t,nanmean(x.ep_raw,1),'LineWidth',1)
plot_beautify
axis tight
xlabel('Time in sec')
ylabel('uV')
title('Detected SO Ch 24')

% store data
grid_so_exp{subj} = grid_so;

%%


grid_sp_exp={};grid_so_exp={};
for subj=1:length(exp_files)
    filename = [ exp_files{subj}]
    EEG = pop_biosig(filename,'importevent','off');
    eeglab redraw
    data= double(EEG.data);

    % filter the data
    data = filtfilt(lpFilt,data')';

    % reference to mastoids (maybe)
    ref = nanmean(data(7:8,:));
    data = data-ref;

    % get good times
    I = data_good_times_exp{subj}';

    % get the appropriate staging file
    suc=[];
    for k=1:length(staging_files)
        if  regexp(staging_files{k},filename(32:36))
            suc = k;
            break
        end
    end


    % get sleep staging data
    if ~isempty(suc)


        disp('processing data')

        % restrict only to N2 sleep stages
        scores = importfile1(staging_files{subj},...
            'Data', [3, 233]);

        % extract the time and day etc.
        time=[];sleep_score=[];
        for ii=1:size(scores,1)
            xx = scores.Score(ii);
            a = [hour(xx) minute(xx) second(xx)];
            time(ii)=a(1)*3600 + a(2)*60 + a(3);
            sleep_score(ii) = scores.VarName2(ii);
        end
        time=time-time(1);


        % get the closest match between the scored data and the EEG data
        sleep_staging = [];
        for ii=1:length(time)
            temp_sleep_scores = sleep_score(ii) * ones(EEG.srate*30,1);
            sleep_staging = [sleep_staging ;temp_sleep_scores];
        end

        if length(sleep_staging) < length(I)
            sleep_staging = [zeros((length(I) - length(sleep_staging)),1) ; sleep_staging];
        end

        %         if length(sleep_staging) < length(I)
        %             sleep_staging(end+1:length(I)) = 0;
        %         end
        %         figure;
        %         plot(sleep_staging)
        %         title(num2str(subj))

        % 26 and 15 are good subjects in terms of sleep staging

        % spindle analyses in all channels
        I = sleep_staging>0;
        grid_sp = detect_spindles(data(1:6,:),I,soFilt,spFilt1,spFilt2,sleep_staging,400);

        % storing file name
        grid_sp.filename = filename;
        grid_sp.sleep_staging = sleep_staging;

        %store data
        grid_sp_exp{subj} = grid_sp;


        % Examining SO-spindle nesting in stage 3 sleep
        grid_so = SO_analyses(data(1:6,:),I,soFilt,spFilt1,spFilt2,sleep_staging);
        grid_so.filename = filename;
        grid_so.sleep_staging = sleep_staging;

        % store data
        grid_so_exp{subj} = grid_so;

    end
end
cd('C:\Users\nikic\OneDrive\Documents\MATLAB\Ana EEG')
save grid_sp_exp grid_sp_exp -v7.3
save grid_so_exp grid_so_exp -v7.3

% extra stuff, plotting delta vs. spindle activity to see if we can do
% sleep staging that way
xx = abs(hilbert(xx));
yy = abs(hilbert(yy));
% compute power in 30s internals
len=EEG.srate*30;
dp=[];
sp=[];
for i=0:len:length(xx)
    if (i+len) < length(xx)
        dp = [dp median(xx(i+1:i+len))];
        sp = [sp median(yy(i+1:i+len))];
    else
        dp = [dp median(xx(i+1:end))];
        sp = [sp median(yy(i+1:end))];
    end
end

figure;hold on
xlim([0 30])
ylim([0 4])
cmap=parula(length(dp));
for i=1:length(dp)
    plot(dp(i),sp(i),'.','MarkerSize',20,'color',cmap(i,:));
    pause(0.005)
end


% plot an example of raw EEG, spindle detected as well as envelope
ch=data(3,:);
chf=filtfilt(spFilt2,ch);
chfa=abs(hilbert(chf));
sp = grid_sp.ch3;
sp_st_idx = sp.sp_st_N2(4);
sp_st = sp.sp_st(sp_st_idx);
sp_end = sp.sp_end(sp_st_idx);
len = sp_end-sp_st;
figure;
subplot(2,1,1)
plot(ch(sp_st-750:sp_end+750))
vline([751 750+len])
subplot(2,1,2)
plot(chf(sp_st-750:sp_end+750));
hold on
plot(chfa(sp_st-750:sp_end+750))
vline([751 750+len])


% plot power spectrum along with the null 1/f spectrum
tmp=grid_sp.ch4.sp_epochs_N2;
P=[];
for i=1:size(tmp,1)
    [Pxx,F]=pwelch(tmp(i,:),[],[],[],400);
    P = [P Pxx];
end
P=log10(abs(P));
Pb =  sort(bootstrp(1000,@mean,P'));
figure;plot(F,(mean(((P)),2)),'LineWidth',2,'Color','k');
hold on
%plot(F,Pb(1,:),'--k','LineWidth',2)
%plot(F,Pb(end,:),'--k','LineWidth',2)
idx=F;
[fillhandle,msg]=jbfill(idx',((Pb(25,:))),((Pb(975,:)))...
    ,[0.3 0.3 0.3],[0.3 0.3 0.3],1,.2);
hold on
xlim([0 20])
xlabel('Frequency')
ylabel('Power')
set(gcf,'Color','w')
set(gca,'FontSize',14)
box off
%plotting the 1/f spectrum after doing a robust regression fit
x=F;
y=(mean(((P)),2));
[aa bb]=find(x<=20);
x=x(aa);
y=y(aa);
[bhat p wh se ci t_stat]=robust_fit(x,y,1);
yhat = bhat(1)+bhat(2).*x;
plot(x,yhat,'r','LineWidth',1)

Pcolor=[];
fidx = logical((F>=10) .* (F<=16));
for i=1:size(P,2)
    tmp=P(fidx,i);
    Pcolor=[Pcolor tmp];
end
figure;
imagesc(F(fidx),1:size(Pcolor,2),Pcolor(:,20:end)')
colormap turbo
set(gcf,'Color','w')

figure;plot(F(fidx),mean(Pcolor,2),'k','LineWidth',1)
box off
set(gcf,'Color','w')
axis tight
xlim([10 16])
ylim([-0.05 0.8])



% plot an example of a nested SO and spindle
so=grid_so.ch3;
sof=filtfilt(soFilt,ch);
t=800;
so_st = so.so_st(so.N3_so_idx);
so_end = so.so_end(so.N3_so_idx);
sp_st = sp.sp_st;
[nested_per, nested_idx] = nested_spindles(so_st,sp_st,t);
idx=19;
len = so_end(idx)-so_st(idx);
figure;
subplot(3,1,1)
plot(ch(so_st(idx)-750:so_end(idx)+750))
vline([751 750+len])
axis tight
box off
subplot(3,1,2)
plot(sof(so_st(idx)-750:so_end(idx)+750))
axis tight
box off
subplot(3,1,3)
plot(chf(so_st(idx)-750:so_end(idx)+750))
hold on
plot(chfa(so_st(idx)-750:so_end(idx)+750))
axis tight
box off



% in ch4, find those SO that have a spindle nested with them
so_st= grid_so.ch4.so_st;
sp_st= grid_sp.ch4.sp_st;
[nested_per1,idx] = nested_spindles(so_st,sp_st,800);

%%%% some extra stuff to do sleep staging and processing of the eeg data
%%%% for the figures for the paper %%%%%%
% spectrogram for channel 3
data=data(1:6,:)';
[c,s]=pca(data);
data = (s(:,2:end)*c(:,2:end)')';

ch3 = data(3,:);

figure;plot(ch3)
[aa bb]=ginput;
aa=round(aa);
m = median(ch3);
s = std(ch3,1)*1;
for i=1:2:length(aa)
    tmp = randn(length(aa(i):aa(i+1)),1)*s;
    ch3(aa(i):aa(i+1)) = m +tmp;
end
figure;plot(ch3)
ylim([-200 200])

Fs=400;
[S,F,T,P,Fc,Tc] = spectrogram(ch3,2048,0,[],Fs,'power');
figure;
idx=find(F>=30);idx=idx(1);
imagesc(T,F(1:idx),log10(abs(P(1:idx,:))));
set(gca,'YDir','normal')

figure;spectrogram(ch3,8192,4096,0.2:0.5:30,Fs,'psd','yaxis');

[tf, freqs, times]= timefreq(ch3', Fs,'cycles',[0],'wletmethod','dftfilt2',...
    'winsize',16384*1,'freqs',[.5 30],'ntimesout',500);
figure;imagesc(times,freqs,log10(abs(tf)))
set(gca,'Ydir','normal')
caxis([1 4.5])

% using the filterbank approach
freqs = [1.5:2:30];
filtered_data=[];
for i=1:length(freqs)
    disp(freqs(i))
    lo = freqs(i)-1;
    up = freqs(i)+1;
    bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
        'HalfPowerFrequency1',lo,'HalfPowerFrequency2',up, ...
        'SampleRate',400);
    tmp = filtfilt(bpFilt,ch3);
    tmp = abs(hilbert(tmp));
    tmp = (tmp);
    filtered_data = [filtered_data;tmp];
end
tt = (0:length(sleep_staging)-1)*(1/Fs)*(1/3600);
%filtered_data = (zscore(filtered_data'))';
figure;imagesc(tt,freqs,abs((filtered_data)))
set(gca,'Ydir','normal')
caxis([0 10])

B = imgaussfilt(filtered_data,'FilterSize',15);
figure;imagesc(tt,freqs,((B)))
set(gca,'Ydir','normal')
caxis([0 10])

% using eegfilt
[smoothdata] = eegfilt(ch3,Fs,0.5,30);

idx = find(sleep_staging==5);
sleep_staging(idx)=4;

% plotting the data
figure;
ha=tight_subplot(2,1);
axes(ha(1))
tt = (0:length(sleep_staging)-1)*(1/Fs)*(1/3600);
plot(tt,sleep_staging,'k','LineWidth',1)
axis tight
axes(ha(2))
plot(tt,ch3,'Color',[.15 .15 .15 .55],'LineWidth',1)
axis tight
ylim([-200 200])
axes(ha(3))
tmp = filtfilt(deltaFilt,ch3);
plot(tt,tmp,'Color',[.15 .15 .15 .55],'LineWidth',1)
axis tight
ylim([-140 140])
axes(ha(4))
tmp = smooth(zscore(abs(hilbert(filtfilt(spFilt1,ch3)))));
plot(tt,tmp,'Color',[.15 .15 .15 .55],'LineWidth',1)
axis tight
ylim([-20 20])
xticklabels ''
box off
set(gcf,'Color','w')
set(gca,'FontSize',14)

% plotting spindle and SO detection stuff - ch3
% spindle detection
ch4=ch3;
figure;
sp_tmp = filtfilt(spFilt1,ch4);
so_tmp = filtfilt(soFilt,ch4);
figure;plot(sp_tmp)
hold on
plot(abs(hilbert(sp_tmp)),'r')
ylim([-6.25 6.25])
idx=grid_sp.ch4.sp_st_N3;
vline(grid_sp.ch4.sp_st(idx),'r')
vline(grid_so.ch4.so_st,'k')
figure;plot(so_tmp)
xlim([7.02 7.06]*1e5)
%xlim([2.6715 2.6747]*1e6)
ylim([-60 40])
aa=hline(0);
set(aa,'LineWidth',1)

% TOPOPLOT
figure;topoplot([],EEG.chanlocs,'style','blank','emarker',{'.','k',40,1});
set(gcf,'Color','w')
axis tight

% PLOTTING SOME PHASE BASED ANALYSES OF COUPLING IN SPINDLE BAND ACTIVITY
% BETWEEN CHANNELS
ch1=grid_sp.ch3.sp_epochs_N2;
ch2=grid_sp.ch4.sp_epochs_N2;
plv=[];
plotting=1;
for i=1:100
    tmp=ch1(i,:);
    tt=linspace(-2,2,length(tmp));
    if plotting==1
        figure;plot(tt,tmp,'k','LineWidth',1)
        axis tight
        set(gcf,'Color','w')
        xlim([-1.5 1.5])
        axis off
    end
    tmp=filtfilt(spFilt2,tmp);
    if plotting==1
        figure;plot(tt,tmp,'k','LineWidth',1)
        axis tight
        set(gcf,'Color','w')
        xlim([-1.5 1.5])
        axis off
    end
    tmp_phase=angle(hilbert(tmp));
    if plotting==1
        figure;plot(tt,tmp_phase,'k','LineWidth',1)
        axis tight
        set(gcf,'Color','w')
        xlim([-1.5 1.5])
        axis off
    end

    tmp2=ch2(i,:);
    if plotting==1
        figure;plot(tt,tmp2,'b','LineWidth',1)
        axis tight
        set(gcf,'Color','w')
        xlim([-1.5 1.5])
        axis off
    end
    tmp2=filtfilt(spFilt2,tmp2);
    if plotting==1
        figure;plot(tt,tmp2,'b','LineWidth',1)
        axis tight
        set(gcf,'Color','w')
        xlim([-1.5 1.5])
        axis off
    end
    tmp2_phase=angle(hilbert(tmp2));
    if plotting==1
        figure;plot(tt,tmp2_phase,'b','LineWidth',1)
        axis tight
        set(gcf,'Color','w')
        xlim([-1.5 1.5])
        axis off
    end

    % plotting phase on top of each other
    if plotting==1
        figure;plot(tt,tmp_phase,'k','LineWidth',1)
        hold on
        plot(tt,tmp2_phase,'b','LineWidth',1)
        axis tight
        xlim([-0.6 0.6])
        axis off
        set(gcf,'Color','w')
    end

    phase_diff = exp(1i*(tmp_phase-tmp2_phase));
    plv(i)= mean(phase_diff);


end
figure;stem(abs(plv))
figure;rose(angle(plv))
% plotting an example phase coupled channel
alpha=angle(phase_diff);
[t,r]=rose(alpha,20);
figure
polarplot(t,r,'LineWidth',1,'Color','k');
pax=gca;
%pax.RLim = [0 20];
thetaticks(0:45:315);
pax.ThetaAxisUnits='radians';
pax.FontSize=16;
set(gcf,'Color','w')
pax.RTick = [5 10 15 20 ];
pax.GridAlpha = 0.25;
pax.MinorGridAlpha = 0.25;
pax.ThetaMinorGrid = 'off';
%pax.ThetaTickLabel = {'0', ' ', '\pi/2 ', ' ','\pi',' ','3\pi/2',' '};
pax.ThetaTickLabel = ''
%pax.ThetaTickLabel = {'0', ' ', ' ', ' ','\pi',' ',' ',' '};
pax.RTickLabel = {' ',' '};
pax.RAxisLocation=1;
pax.RAxis.LineWidth=1;
pax.ThetaAxis.LineWidth=1;
pax.LineWidth=1;
temp = exp(1i*alpha);
r1 = abs(mean(temp))*max(2*r);
phi = angle(mean(temp));
hold on;
polarplot([phi-0.01 phi],[0 r1],'LineWidth',1.5,'Color','r')
set(gcf,'PaperPositionMode','auto')
set(gcf,'Position',[680.0,865,120.0,113.0])



% BASELINE BLOCK
clc;clear
cd('C:\Users\nikic\OneDrive\Documents\MATLAB\Ana EEG')
load baseline_analyses
cd('C:\Users\nikic\OneDrive\Documents\MATLAB\eeglab2019_0');
eeglab;
close

% low pass filters
lpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',30, ...
    'SampleRate',400);

bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',65,'HalfPowerFrequency2',100, ...
    'SampleRate',400);

% RE-RUN WITH NEW SPINDLE PARAMETERS FROM 11-13 and 13-16
spFilt1 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',11,'HalfPowerFrequency2',13, ...
    'SampleRate',400);

spFilt2 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',16, ...
    'SampleRate',400);

soFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.16,'HalfPowerFrequency2',1.25, ...
    'SampleRate',400);




sleep_staging_folder = 'F:\DATA\EEG Data Anne Richards\Sleep Stage Reports\baseline';
staging_files = findfiles('',sleep_staging_folder)';

grid_sp_baseline={};grid_so_baseline={};
for subj=46:length(baseline_files)
    filename = [ baseline_files{subj}];
    EEG = pop_biosig(filename,'importevent','off');
    eeglab redraw
    data= double(EEG.data);

    % filter the data
    %data = filtfilt(lpFilt,data')';

    % reference to mastoids
    ref = nanmean(data(7:8,:));
    data = data-ref;

    % get good times
    I = data_good_times_baseline{subj}';

    % get the appropriate staging file
    suc=[];
    for k=1:length(staging_files)
        if  regexp(staging_files{k},filename(32:36))
            suc = k;
            break
        end
    end

    % get sleep staging data
    if ~isempty(suc)



        disp('processing data')


        % restrict only to N2 sleep stages
        scores = importfile1(staging_files{suc},...
            'Data', [3, 233]);

        % extract the time and day etc.
        time=[];sleep_score=[];
        for ii=1:size(scores,1)
            xx = scores.Score(ii);
            a = [hour(xx) minute(xx) second(xx)];
            time(ii)=a(1)*3600 + a(2)*60 + a(3);
            sleep_score(ii) = scores.VarName2(ii);
        end
        time=time-time(1);


        % get the closest match between the scored data and the EEG data
        sleep_staging = [];
        for ii=1:length(time)
            temp_sleep_scores = sleep_score(ii) * ones(EEG.srate*30,1);
            sleep_staging = [sleep_staging ;temp_sleep_scores];
        end

        if length(sleep_staging) < length(I)
            sleep_staging(end+1:length(I)) = 0;
        end

        % spindle analyses in all channels
        I = (sleep_staging>0) .* (sleep_staging~=5);
        %I(1:1000)=0;
        grid_sp = detect_spindles(data(1:6,:),I,soFilt,spFilt1,spFilt2,sleep_staging,400);

        % storing file name
        grid_sp.filename = filename;
        grid_sp.sleep_staging = sleep_staging;

        %store data
        grid_sp_baseline{subj} = grid_sp;

        % Examining SO-spindle nesting in stage 3 sleep
        grid_so = SO_analyses(data(1:6,:),I,soFilt,spFilt1,spFilt2,sleep_staging);
        grid_so.filename = filename;
        grid_so.sleep_staging = sleep_staging;

        % store data
        grid_so_baseline{subj} = grid_so;
    end
end

cd('C:\Users\nikic\OneDrive\Documents\MATLAB\Ana EEG')
save grid_sp_baseline grid_sp_baseline -v7.3
save grid_so_baseline grid_so_baseline -v7.3

