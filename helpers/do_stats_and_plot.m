function out = do_stats_and_plot(chdata_neutral,chdata_trauma,EEG)



% spatiotemporal cluster correction
%addpath('C/home/user/Documents/Repositories/limo_tools/limo_tools')
addpath(genpath('/home/user/Documents/Repositories/limo_tools/'))
tfce_flag = false;
loop_iter=1000;
t_scores=[];tboot=zeros(62,7000,loop_iter);
p_scores=[];pboot=zeros(62,7000,loop_iter);
parpool('threads')
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
delete(gcp)
chMap=1:62;

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
    aa=aa./max(aa);% normalize
    subplot(3,1,3);
    topoplot((aa),EEG.chanlocs,'maplimits', [-1 1])
    axis tight
    title('Relative Channel importance')
    axis off
    set(gcf,'Color','w')



    %%% no idea what this is for
    % sig_ch_idx = find(aa>0);
    % sig_ch = zeros(numel(chMap),1);
    % sig_ch(sig_ch_idx)=1;
    % sig_ch_all(i,:) = sig_ch;
    % subplot(3,1,3);
    % imagesc(sig_ch(chMap))
    % title('Sig channels 0 to 3s')
    % %sgtitle(ImaginedMvmt{i})
    % axis off
    % set(gcf,'Color','w')
    % colorbar
end



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


end