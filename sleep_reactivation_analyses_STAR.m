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

%% SAVING ALL THE SLEEP DATA AS SET FROM POLY5

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
cd(root_path)

eeglab
err_files=[];
for subj=28:42
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    eeglab redraw
    subj_name = ['Participant 240' num2str(subj)];
    filepath = fullfile(root_path,subj_name);
    if isfolder(filepath)
        cd(filepath)
        files = dir('*.Poly5');
        nap_file = files(contains({files.name}, 'nap', 'IgnoreCase', true));
        if ~isempty(nap_file)
            nap_file = fullfile(nap_file.folder,nap_file.name);
            Poly5toEEGlab(nap_file);
            %load into eeglab
            [pathname,filename,extension] = fileparts(nap_file);
            filename = [filename(1:end-5) '.set'];
            EEG = pop_loadset('filename',filename,'filepath',...
                pathname);
            eeglab redraw

            setname = ['240' num2str(subj) '_SleepData'];
            setname1 = [setname '.set'];
            EEG = pop_editset(EEG, 'setname', setname);
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            eeglab redraw
            EEG = pop_saveset( EEG, 'filename',...
                setname1,'filepath',pathname);
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            eeglab redraw
        else
            err_files=[err_files;nap_file];
        end
    end
end

%% REACTIVATION SUBJECT SPECIFIC ENGRAM AND AVERAGE ENGRAM (MAIN MAIN)

%clc;clear
%close all
%eeglab
cd('/media/user/Data/Ana EEG/STAR/Phase 2')
root_path=pwd;
load Imag_Snippet_SubjMatched

% get subjects
subj_names={};
for i=1:length(files_imag)
    tmp=files_imag{i};
    subj_names{i}=tmp(67:71);
end

% % reactivation of subject-specific engram (not the average)
% % loop over subjects
% for i=1:length(subj)
%     % get imag data
%     snippet_data = squeeze(chdata_trauma(:,:,i) - chdata_neutral(:,:,i));
%     imag_data = squeeze(chdata_trauma_imag(:,:,i) - chdata_neutral_imag(:,:,i));
% end


% reactivation of the engram derived from average data
d1 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',35, ...
    'SampleRate',1e3);
load dPCA_Shared_AllSubj
load dPCA_Shared_AllSubj_baseline
res_react=[];
res_react_stage=[];
subj_proc={};
subj_not_proc={};
for i=1:length(subj_names)
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    eeglab redraw

    % load the sleep data
    subj = ['Participant ' subj_names{i}];
    nap_file = [subj_names{i} '_SleepData.set'];
    filepath = fullfile(root_path,subj);
    napfile_full = fullfile(filepath,nap_file);
    if isfolder(filepath) && i~=3 && exist(napfile_full)
        EEG = pop_loadset('filename',nap_file,...
            'filepath',filepath);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        eeglab redraw

        % process sleep data
        data = EEG.data;
        data1 = filtfilt(d1,data')';
        ref = data1([13 19],:); % reference to M1 and M2, channels 1 to 68
        ref = mean(ref,1);
        data1(1:68,:) = data1(1:68,:) - ref;
        EEG.data = data1;
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        eeglab redraw

        EEG = pop_select( EEG, 'rmchannel',{'M1','M2','SpO2','PLETH','HRate','Stat','TRIGGERS1',...
            'TRIGGERS2','TRIGGERS3','TRIGGERS4','TRIGGERS5','TRIGGERS6','TRIGGERS7','TRIGGERS8',...
            'TRIGGERS9','TRIGGERS10','TRIGGERS11','TRIGGERS12','TRIGGERS13','TRIGGERS14','TRIGGERS15',...
            'TRIGGERS16','STATUS','Counter 2power24','BIP 01','BIP 02','BIP 03','BIP 04'});% remove bad channels
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off');
        %[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        eeglab redraw
        data=EEG.data;


        % load the sleep staging
        sleep_report_path = fullfile(filepath,'Sleep Stage Report');
        sleep_file = [subj_names{i} '_SleepStages.xls'];
        sleep_filename = fullfile(sleep_report_path,sleep_file);
        if isfolder(sleep_report_path) && exist(sleep_filename)

            sleep_score = import_nap_sleep_scores(sleep_filename, "Data", [2, Inf]);
            sleep_score = sleep_score(2:end,:);
            sleep_staging= [];
            tidx = 1;
            for ii=1:height(sleep_score)
                tmp = sleep_score(ii,:).Var2;
                sleep_staging(tidx:(tidx+30*EEG.srate-1)) = tmp;
                tidx = size(sleep_staging,2) + 1;
            end

            % align the data
            data = data(:,1:length(sleep_staging));

            % get the reactivation
            [n2,n3,vaf_task_shared,vaf_base_shared] = get_reactivation_sleep_stage_subject(data,sleep_staging,...
                whichMarg_total,whichMarg_total_baseline,Wtotal,Vtotal,Wtotal_baseline,Vtotal_baseline);

            % store
            res_react =  [res_react; [nanmean([vaf_task_shared vaf_base_shared])]];
            res_react_stage = [res_react_stage; [nanmean(n2) nanmean(n3)]];
            subj_proc = [subj_proc  subj_names{i}];
        else
            subj_not_proc = [subj_not_proc  subj_names{i}];
        end
    else
        subj_not_proc = [subj_not_proc  subj_names{i}];
    end
end

save reactivation_analyses_NREM res_react res_react_stage ...
    subj_proc subj_not_proc  -v7.3


fprintf('\nProcessed %d number of subjects\n', length(subj_proc));

figure;boxplot(res_react)
xticks(1:2)
ylim([-3 1])
xticklabels({'Engram','Baseline'})
plot_beautify
ylabel('VAF')
[p,h,stats]=signrank(res_react(:,1),res_react(:,2));
title(['Reactivation pval ' num2str(p)])

figure;boxplot(res_react_stage)
xticks(1:2)
ylim([-1 0.8])
xticklabels({'N2','N3'})
plot_beautify
ylabel('VAF')
[p,h,stats]=signrank(res_react_stage(:,1),res_react_stage(:,2));
title(['Reactivation pval ' num2str(p)])

%% DETECTING IF SHARED SIGNAL IS REPLAYED IN SLEEP (MAIN)

[
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

cd('/media/user/Data/Ana EEG/STAR/Phase 2/Participant 24016/')
eeglab
Poly5toEEGlab

% load the dataset
% EEG = pop_loadset('filename','24014_snippet-20250226T125238.set','filepath',...
%     'F:\DATA\EEG Data Anne Richards\STAR\Phase 2\Participant 24014\');
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );


% call the functions to load the SAGA datafile
%Poly5toEEGlab
%eeglab redraw

%%%% perform sleep analyses filter and clean data
data = EEG.data;
d1 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',35, ...
    'SampleRate',EEG.srate);
data1 = filtfilt(d1,data')';

% reference to M1 and M2, channels 1 to 68
ref = data1([13 19],:);
ref = mean(ref,1);
data1(1:68,:) = data1(1:68,:) - ref;
EEG.data = data1;
eeglab redraw

% remove bad channels 
EEG = pop_select( EEG, 'rmchannel',{'M1','M2','SpO2','PLETH','HRate','Stat','TRIGGERS1',...
    'TRIGGERS2','TRIGGERS3','TRIGGERS4','TRIGGERS5','TRIGGERS6','TRIGGERS7','TRIGGERS8',...
    'TRIGGERS9','TRIGGERS10','TRIGGERS11','TRIGGERS12','TRIGGERS13','TRIGGERS14','TRIGGERS15',...
    'TRIGGERS16','STATUS','Counter 2power24','BIP 01','BIP 02','BIP 03','BIP 04'});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off');
eeglab redraw

% get the eeg data
data=EEG.data;

% now getting the sleep staging files and aligning it
filename = fullfile(pwd,'Sleep Stage Report/24016_SleepStages.xls');
sleep_score = import_nap_sleep_scores(filename, "Data", [2, Inf]);
sleep_score = sleep_score(2:end,:);

sleep_staging= [];
tidx = 1;
for i=1:height(sleep_score)
    tmp = sleep_score(i,:).Var2;
    sleep_staging(tidx:(tidx+30*EEG.srate-1)) = tmp;
    tidx = size(sleep_staging,2) + 1;
end
figure;
plot(sleep_staging)

data = data(:,1:length(sleep_staging));



%% Sleep reactivation: task shared dPCs vs baseline shared dPCs (MAIN)
%(MAIN)
% Sleep data:
%   data = channels x samples, e.g. 62 x 7205598
%
% Task dPCA:
%   Wtotal, Vtotal, whichMarg_total
%
% Baseline dPCA:
%   Wtotal_baseline, Vtotal_baseline, whichMarg_total_baseline
%
% Assumption:
%   whichMarg == 1 : condition-specific dPCs
%   whichMarg == 2 : shared time dPCs

clear vaf_task_shared vaf_base_shared

p = gcp('nocreate');
if isempty(p)
    parpool('threads')
end

%%%%%

% Parameters
Fs = 1000;
winSec = 3;
winSamp = winSec * Fs;
stepSamp = winSamp / 2;      % 50% overlap
nDim = 11;                   % use first 11 dPCs / 95% VAF dimensionality
Xsleep = double(data);       % channels x samples

[nChan, nSamples] = size(Xsleep);
startIdx = 1:stepSamp:(nSamples - winSamp + 1);
nWins = numel(startIdx);
fprintf('Sleep duration = %.2f min\n', nSamples/Fs/60);
fprintf('Number of 3 s windows with 50%% overlap = %d\n', nWins);

% Select shared/time dPCs from task and baseline
taskShared_all = find(whichMarg_total == 2);
baseShared_all = find(whichMarg_total_baseline == 2);

% Restrict to first nDim dPCs
taskShared = taskShared_all(taskShared_all <= nDim);
baseShared = baseShared_all(baseShared_all <= nDim);

% Match number of components for fair comparison
k = min(numel(taskShared), numel(baseShared));
taskShared = taskShared(1:k);
baseShared = baseShared(1:k);

fprintf('Using %d shared/time dPCs.\n', k);
fprintf('Task shared dPCs: ');
disp(taskShared);
fprintf('Baseline shared dPCs: ');
disp(baseShared);

% Preallocate
vaf_task_shared = nan(nWins,1);
vaf_base_shared = nan(nWins,1);
segment_sleep_stages = nan(nWins,1);

% Loop through sleep snippets
parfor w = 1:nWins

    idx1 = startIdx(w);
    idx2 = idx1 + winSamp - 1;

    Xwin = Xsleep(:, idx1:idx2);      % channels x 3000 samples

    % Center each sleep window channel-wise
    Xwin = Xwin - mean(Xwin,2);

     % get sleep stage
    tmp_sleep_stage = sleep_staging(idx1:idx2);
    segment_sleep_stages(w) = mode(tmp_sleep_stage);

    denom = norm(Xwin,'fro')^2;

    % % error correction
    % if denom == 0 || any(isnan(Xwin(:)))
    %     continue
    % end

    % Reconstruction using engram
    Z_task = Wtotal(:,taskShared)' * Xwin;
    Xhat_task = Vtotal(:,taskShared) * Z_task;
    vaf_task_shared(w) = 1 - norm(Xwin - Xhat_task,'fro')^2 / denom;

    % Reconstruction using baseline  dPCs
    Z_base = Wtotal_baseline(:,baseShared)' * Xwin;
    Xhat_base = Vtotal_baseline(:,baseShared) * Z_base;
    vaf_base_shared(w) = 1 - norm(Xwin - Xhat_base,'fro')^2 / denom;
end


%%%% SORTING RESULTS BY SLEEP STAGE 
n2=[];n3=[];
idx2 = find(segment_sleep_stages == 2);
idx3 = find(segment_sleep_stages == 3);
n2 = vaf_task_shared(idx2); %- vaf_base_shared(idx2);
n3 = vaf_task_shared(idx3);% - vaf_base_shared(idx3);
if length(n2) > length(n3)
    n3(end+1:length(n2)) = NaN;
elseif length(n2) < length(n3)
    n2(end+1:length(n3)) = NaN;
end

figure;boxplot([n2 n3])
ylim([-0.5 1.5])
xticks(1:2)
xticklabels({'N2','N3'})
ylabel('Sleep reactivation')

%%% PLOTTING

% Time axis
tMin = (startIdx + winSamp/2 - 1) / Fs / 60;
% Plot VAF over sleep
figure; hold on
plot(tMin, smooth(vaf_task_shared,10), 'LineWidth', 1.5);
plot(tMin, smooth(vaf_base_shared,10), 'LineWidth', 1.5);
plot(tMin, segment_sleep_stages/3, 'k','LineWidth', 1.5);
ylim([-2 2])
xlabel('Sleep time (min)');
ylabel('VAF');
legend({'Engram', 'Baseline',...
    'Sleep Stage'}, ...
    'Location','best');
title('Sleep reactivation of Engram');
box off;

% Smooth plot for visualization only
smoothWin = 25; % windows; 25 windows = ~37.5 s with 50% overlap
vaf_task_smooth = movmean(vaf_task_shared, smoothWin, 'omitnan');
vaf_base_smooth = movmean(vaf_base_shared, smoothWin, 'omitnan');
figure; hold on
plot(tMin, vaf_task_smooth, 'LineWidth', 2);
plot(tMin, vaf_base_smooth, 'LineWidth', 2);
xlabel('Sleep time (min)');
ylabel('Smoothed VAF');
legend({'Engram', 'Non-task baseline structure'}, ...
    'Location','best');
xlim([0 120])
title('Sleep reactivation');
box off;
plot_beautify

% Summary statistics
fprintf('\nTask shared/time dPC sleep VAF:\n');
fprintf('Mean ± SEM = %.4f ± %.4f\n', ...
    mean(vaf_task_shared,'omitnan'), ...
    std(vaf_task_shared,'omitnan') / sqrt(sum(~isnan(vaf_task_shared))));
fprintf('Median = %.4f\n', median(vaf_task_shared,'omitnan'));

fprintf('\nBaseline shared/time dPC sleep VAF:\n');
fprintf('Mean ± SEM = %.4f ± %.4f\n', ...
    mean(vaf_base_shared,'omitnan'), ...
    std(vaf_base_shared,'omitnan') / sqrt(sum(~isnan(vaf_base_shared))));
fprintf('Median = %.4f\n', median(vaf_base_shared,'omitnan'));

% Direct comparison: task shared vs baseline shared
diff_vaf = vaf_task_shared - vaf_base_shared;
fprintf('\nTask shared - baseline shared VAF:\n');
fprintf('Mean ± SEM = %.4f ± %.4f\n', ...
    mean(diff_vaf,'omitnan'), ...
    std(diff_vaf,'omitnan') / sqrt(sum(~isnan(diff_vaf))));
fprintf('Median = %.4f\n', median(diff_vaf,'omitnan'));

% Note: windows overlap, so this p-value treats windows as observations
% and is mainly descriptive.
[p_win,~,stats_win] = signrank(vaf_task_shared, vaf_base_shared);
fprintf('\nWindow-wise signed-rank test, task shared vs baseline shared:\n');
fprintf('p = %.4g, signedrank = %.3f\n', p_win, stats_win.signedrank);







% Better summary using non-overlapping blocks for less inflated statistics
% Each block is 1 min of sleep.
blockMin = 1;
blockDurSamp = blockMin * 60 * Fs;
nBlocks = floor(nSamples / blockDurSamp);
block_task = nan(nBlocks,1);
block_base = nan(nBlocks,1);
for b = 1:nBlocks
    blockStartSamp = (b-1)*blockDurSamp + 1;
    blockEndSamp   = b*blockDurSamp;
    winInBlock = startIdx >= blockStartSamp & startIdx <= (blockEndSamp - winSamp + 1);
    block_task(b) = mean(vaf_task_shared(winInBlock), 'omitnan');
    block_base(b) = mean(vaf_base_shared(winInBlock), 'omitnan');
end

[p_block,~,stats_block] = signrank(block_task, block_base);
fprintf('\n1-min block signed-rank test:\n');
fprintf('Task shared mean block VAF = %.4f\n', mean(block_task,'omitnan'));
fprintf('Baseline shared mean block VAF = %.4f\n', mean(block_base,'omitnan'));
fprintf('p = %.4g, signedrank = %.3f\n', p_block, stats_block.signedrank);

% Boxplot of 1-min block means
figure; hold on
boxplot([block_task, block_base], ...
    'Labels', {'Engram', 'Baseline'});
ylabel('Mean VAF per 1-min sleep block');
title('Sleep reactivation: task shared vs baseline shared dPCs');
box off;
xlim([0.5 2.5])
plot_beautify


figure; hold on
data1 = [block_task, block_base];
figure; hold on
boxplot(data1, ...
    'Labels', {'Engram','Baseline'}, ...
    'Symbol', '', ...        % hides outliers
    'Widths', 0.5);

% Get box handles
h = findobj(gca, 'Tag', 'Box');
h = flipud(h);   % makes order match columns

% Colors
colors = [
    0.2 0.45 0.85;   % box 1
    0.85 0.35 0.35   % box 2
];

% Shade boxes
for i = 1:2
    patch(get(h(i),'XData'), get(h(i),'YData'), colors(i,:), ...
        'FaceAlpha', 0.35, ...
        'EdgeColor', colors(i,:), ...
        'LineWidth', 1.5);
end

% Put median lines on top
med = findobj(gca, 'Tag', 'Median');
set(med, 'Color', 'k', 'LineWidth', 1.5);

ylabel('Mean VAF per 1-min sleep block');
title('Sleep reactivation');

box off
set(gca, 'FontSize', 12);
xlim([.5 2.5])
plot_beautify



%% Sleep replay with 50% EEG segment overlap: (comparing condn vs. time vs. baseline)
% condition-specific dPCs vs time dPCs vs random manifolds
%
% data            = sleep EEG, channels x samples
% Fs              = 1000
% Wtotal          = dPCA decoder axes, channels x components
% Vtotal          = dPCA encoder axes, channels x components
% whichMarg_total = marginalization labels
%
% Assumes:
% whichMarg_total == 1 : condition-specific dPCs
% whichMarg_total == 2 : time dPCs

p = gcp('nocreate');
if isempty(p)
    parpool('threads')
end

% Parameters
Fs = 1000;
winSec = 3;
winSamp = winSec * Fs;
stepSamp = winSamp/2;       % 50% overlap = 1.5 s step

nDim = 11;
nNull = 5;

Xsleep = double(data);      % channels x samples

[nChan, nSamples] = size(Xsleep);

startIdx = 1:stepSamp:(nSamples - winSamp + 1);
nWins = numel(startIdx);
fprintf('Number of 3s windows with 50%% overlap = %d\n', nWins);

% Select dPCs
condComps_all = find(whichMarg_total == 1);
timeComps_all = find(whichMarg_total == 2);

% Restrict to first 11 dPCs
condComps = condComps_all(condComps_all <= nDim); % differentiates viewing and imagined
timeComps = timeComps_all(timeComps_all <= nDim); % engram
nCondComps = numel(condComps);
fprintf('Condition-specific comps:\n');
disp(condComps);
fprintf('Time comps:\n');
disp(timeComps);

% Preallocate
vaf_cond = nan(nWins,1);
vaf_time = nan(nWins,1);
vaf_null = nan(nWins,nNull);
segment_sleep_stages = nan(nWins,1);

% Main loop
parfor w = 1:nWins

    idx1 = startIdx(w);
    idx2 = idx1 + winSamp - 1;

    Xwin = Xsleep(:,idx1:idx2);       % channels x time
    Xwin = Xwin - mean(Xwin,2);       % center each channel

    % get sleep stage
    tmp_sleep_stage = sleep_staging(idx1:idx2);
    segment_sleep_stages(w) = mode(tmp_sleep_stage);

    % ovrall variance
    denom = norm(Xwin,'fro')^2;

    % % check if error
    % if denom == 0 || any(isnan(Xwin(:)))
    %     continue
    % end

    % Condition-specific dPC VAF
    Zcond = Wtotal(:,condComps)' * Xwin;
    Xhat_cond = Vtotal(:,condComps) * Zcond;
    vaf_cond(w) = 1 - norm(Xwin - Xhat_cond,'fro')^2 / denom;

    % Time dPC VAF (engram)
    Ztime = Wtotal(:,timeComps)' * Xwin;
    Xhat_time = Vtotal(:,timeComps) * Ztime;
    vaf_time(w) = 1 - norm(Xwin - Xhat_time,'fro')^2 / denom;

    % baseline manifold
    Z_base = Wtotal_baseline(:,baseShared)' * Xwin;
    Xhat_base = Vtotal_baseline(:,baseShared) * Z_base;
    vaf_base_shared(w) = 1 - norm(Xwin - Xhat_base,'fro')^2 / denom;

    % Random manifold null
    % Same dimensionality as condition-specific dPC subspace
    for p = 1:nNull
        A = randn(nChan, nCondComps);
        [Q,~] = qr(A,0);              % random orthonormal manifold
        Xhat_null = Q * Q' * Xwin;
        vaf_null(w,p) = 1 - norm(Xwin - Xhat_null,'fro')^2 / denom;
    end
end

% Null confidence intervals per sleep window
null_mean = mean(vaf_null,2,'omitnan');
null_lo   = prctile(vaf_null,2.5,2);
null_hi   = prctile(vaf_null,97.5,2);

% Time axis
tMin = (startIdx + winSamp/2 - 1) / Fs / 60;

% Plot
figure; hold on

% Null 95% interval
fill([tMin fliplr(tMin)], ...
     [null_lo' fliplr(null_hi')], ...
     [0.8 0.8 0.8], ...
     'FaceAlpha',0.4, ...
     'EdgeColor','none');

plot(tMin, null_mean, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
plot(tMin, smooth(vaf_cond,10), 'LineWidth', 1.5);
plot(tMin, smooth(vaf_time,10), 'LineWidth', 1.5);
% sleep_stage
plot(tMin,segment_sleep_stages/3, 'LineWidth', 1.5)

xlabel('Sleep time (min)');
ylabel('VAF');
legend({'Random manifold 95% CI', ...
        'Random mean', ...
        'Condition-distinguishing submanifold', ...
        'Engram submanifold',...
        'Sleep stages'}, ...
        'Location','best');

title('Sleep reactivation - Single Participant');
box off;
xlim([0 110])
plot_beautify
ylim([-1 0.75])


%%%% boxplots
data1 = [vaf_cond vaf_time];
figure; hold on
boxplot(data1, ...
    'Labels', {'Condition-specific','Engram'}, ...
    'Symbol', '', ...        % hides outliers
    'Widths', 0.5);

% Get box handles
h = findobj(gca, 'Tag', 'Box');
h = flipud(h);   % makes order match columns

% Colors
colors = [
    0.2 0.45 0.85;   % box 1
    0.85 0.35 0.35   % box 2
];

% Shade boxes
for i = 1:2
    patch(get(h(i),'XData'), get(h(i),'YData'), colors(i,:), ...
        'FaceAlpha', 0.35, ...
        'EdgeColor', colors(i,:), ...
        'LineWidth', 1.5);
end

% Put median lines on top
med = findobj(gca, 'Tag', 'Median');
set(med, 'Color', 'k', 'LineWidth', 1.5);

ylabel('VAF');
box off
set(gca, 'FontSize', 12);
xlim([.5 2.5])
plot_beautify

%%%% reactivation sleep stage
n2=[];n3=[];
idx2 = find(segment_sleep_stages == 2);
idx3 = find(segment_sleep_stages == 3);
n2 = -vaf_cond(idx2) + vaf_time(idx2);
n3 = -vaf_cond(idx3) + vaf_time(idx3);
if length(n2) > length(n3)
    n3(end+1:length(n2)) = NaN;
elseif length(n2) < length(n3)
    n2(end+1:length(n3)) = NaN;
end
figure;boxplot([n2 n3])
ylim([0 2])
xticks(1:2)
xticklabels({'N2','N3'})
ylabel('Sleep reactivation')

%%% Summary statistics
fprintf('\nCondition-specific dPC sleep VAF:\n');
fprintf('Mean ± SEM = %.4f ± %.4f\n', ...
    mean(vaf_cond,'omitnan'), ...
    std(vaf_cond,'omitnan')/sqrt(sum(~isnan(vaf_cond))));
fprintf('Median = %.4f\n', median(vaf_cond,'omitnan'));

fprintf('\nEngram dPC sleep VAF:\n');
fprintf('Mean ± SEM = %.4f ± %.4f\n', ...
    mean(vaf_time,'omitnan'), ...
    std(vaf_time,'omitnan')/sqrt(sum(~isnan(vaf_time))));
fprintf('Median = %.4f\n', median(vaf_time,'omitnan'));

fprintf('\nRandom manifold null VAF:\n');
fprintf('Mean null = %.4f\n', mean(null_mean,'omitnan'));

%%% Window-wise significance versus random null
p_cond = nan(nWins,1);
p_time = nan(nWins,1);

for w = 1:nWins
    p_cond(w) = (1 + sum(vaf_null(w,:) >= vaf_cond(w))) / (1 + nNull);
    p_time(w) = (1 + sum(vaf_null(w,:) >= vaf_time(w))) / (1 + nNull);
end

sig_cond = p_cond < 0.05;
sig_time = p_time < 0.05;

fprintf('\nFraction significant windows, condition-specific vs null: %.3f\n', ...
    mean(sig_cond,'omitnan'));

fprintf('Fraction significant windows, time dPCs vs null: %.3f\n', ...
    mean(sig_time,'omitnan'));

% Compare condition-specific vs time dPC VAF across windows
[p_cond_vs_time,~,stats_cond_vs_time] = signrank(vaf_cond, vaf_time);
fprintf('\nCondition-specific vs Time dPC VAF across windows:\n');
fprintf('p = %.4g, signedrank = %.3f\n', ...
    p_cond_vs_time, stats_cond_vs_time.signedrank);


