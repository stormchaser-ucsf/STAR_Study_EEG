function [n2,n3] = get_reactivation_sleep_stage(data,sleep_staging,...
                    whichMarg_total,whichMarg_total_baseline)



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

% fprintf('Using %d shared/time dPCs.\n', k);
% fprintf('Task shared dPCs: ');
% disp(taskShared);
% fprintf('Baseline shared dPCs: ');
% disp(baseShared);

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
n2 = vaf_task_shared(idx2) - vaf_base_shared(idx2);
n3 = vaf_task_shared(idx3) - vaf_base_shared(idx3);
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


