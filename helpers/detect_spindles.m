function grid_sp = detect_spindles(data_main,I,soFilt,spFilt1,spFilt2,sleep_staging,Fs,slow_fast_idx)
%function grid_sp = detect_spindles(data_main,I,bpFilt,spFilt1,spFilt2,,sleep_staging)

% identify the ied times large high freq spikes
I=logical(I);
II=[];
win = Fs*1; % discarding 1s of data either side of an artifact
for i=1:size(data_main,1)
    % threshold
    II(i,:) = (I);
    xx=find(abs(data_main(i,:))>250);
    for j=1:length(xx)
        if xx(j)+win < length(I) && xx(j)>win+1
            II(i,[-win:win]+xx(j))=0;
        elseif xx(j)+win>length(I)
            II(i,xx(j):end)=0;
        end
    end
end


grid_sp=struct;
for i=1:size(data_main,1)
    disp(['Processing Channel ' num2str(i)])
    I=II(i,:);

    %if i<=2

    if slow_fast_idx(i) == 0
        temp1 = filtfilt(spFilt1,data_main(i,:));
    else
        temp1 = filtfilt(spFilt2,data_main(i,:));
    end

    % removing artifiacts in the spindle amplitudes
    xx=find(abs(temp1)>50);
    win=1*Fs;
    for j=1:length(xx)
        if xx(j)+win < length(I) && xx(j)>win+1
            II(i,[-win:win]+xx(j))=0;
        elseif xx(j)+win>length(I)
            II(i,xx(j):end)=0;
        end
    end

    I=II(i,:);

    % spindle detection
    slo = filtfilt(soFilt,data_main(i,:));
    temp = abs(hilbert(temp1));
    sp_amp = temp(logical(I));



    sp_thresh1 = mean(sp_amp) + 3*std(sp_amp);
    sp_thresh2 = mean(sp_amp) + 1.0*std(sp_amp);

    % find the peaks and peak locations in NREM sleep, given by index I
    [peak_val peak_loc] = findpeaks(temp);

    %indices when peaks above 3 S.D.
    amp_times = find(peak_val>sp_thresh1);

    % timings of these peaks in the original datavector
    peak_I = zeros(1,length(temp));
    peak_I(peak_loc(amp_times)) = 1;
    peak_I = peak_I .* I;
    aa= find(peak_I==1);

    % check if the duration criteria are met (300ms to 3500ms) with at
    % least 1 S.D
    sp_peak=[];sp_st=[];sp_end=[];
    for j=1:length(aa)
        % get the time when it went above 2nd thresh pre peak
        k = aa(j);
        while temp(k) > sp_thresh2 && k>1
            k=k-1;
        end
        st = k+1;

        % get the time when it went above 2nd thresh post peak
        k = aa(j);
        while temp(k) > sp_thresh2 && k < length(I)
            k=k+1;
        end
        stp = k-1;
        dur = stp-st;

        % getting number of peaks in ongoing spindle freq event
        [peak_val1 peak_loc1] = findpeaks(temp1(st:stp));
        no_of_peaks = length(peak_val1);

        if dur*(1000/Fs) > 300 && dur*(1000/Fs)<2000 && no_of_peaks>=3 && st > 1000
            %sp_peak = [sp_peak aa(j)];
            %sp_epochs = [sp_epochs; temp1(aa(j)-220:aa(j)+220)];
            sp_st=[sp_st st];
            sp_end = [sp_end stp];
        end
    end

    [sp_st sp_end] = combine_events(sp_st,sp_end,300,2000,Fs);
    name = ['ch' num2str(i)];
    grid_sp.(name).sp_st=sp_st;
    grid_sp.(name).sp_end=sp_end;

    % z-score temp1
    %m=mean(temp1(logical(I)));
    %s=std(temp1(logical(I)));
    %temp1=(temp1-m)./s;

    sp_epochs=[];sp_epochs_tf=[];x=[];so_epochs=[];fa_epochs=[];
    sp_sp_epochs=[];
    win = 2*Fs; %2s of data
    for j=1:length(sp_st)
        data_seg = temp1(sp_st(j):sp_end(j));
        [a b]= max(data_seg);
        b=sp_st(j)+b-1;
        %sp_epochs = [sp_epochs; data_main(i,b-win:b+win)];
        if b+win <= length(I)
            sp_epochs = [sp_epochs; data_main(i,b-win:b+win)];
            sp_sp_epochs = [sp_sp_epochs;temp1(b-win:b+win)];
        end
        data_seg = temp(sp_st(j):sp_end(j));
        if b+win <= length(I)
            [a b]= max(data_seg);
            b=sp_st(j)+b-1;
            sp_epochs_tf = [sp_epochs_tf;data_main(i,b-win:b+win)];
            x = [ x;temp(b-win:b+win)];
            so_epochs=[so_epochs;slo(b-win:b+win)];
        end
    end

    grid_sp.(name).sp_epochs=sp_epochs;
    grid_sp.(name).so_epochs = so_epochs;
    grid_sp.(name).fa_epochs = fa_epochs;
    grid_sp.(name).I = I;
    grid_sp.(name).sp_sp_epochs = sp_sp_epochs;
    if slow_fast_idx(i) == 0
        grid_sp.(name).type = 'Slow spindle';
    else
        grid_sp.(name).type = 'Fast spindle';
    end

    if i==13 || i==19
        grid_sp.(name).type = 'Reference';
    end

    %grid_sp.(name).x = x;

    %grid_sp.(name).sp_epochs_tf=sp_epochs_tf;

    %     [tf, freqs, times] = timefreq(sp_epochs_tf', Fs,'tlimits',[-2 2],...
    %         'freqs',[0 30],'cycles',[.1 .1]);
    %
    %     %grid_sp.(name).tf=tf;
    %     x=log10(abs(tf));
    %     x=mean(x,3);
    %     %x=mean(tf,3);
    %     %x=abs(x);
    %     fi = (freqs>7) .* (freqs<17);
    %     ti = (times>-0.5) .* (times<0.5);
    %     x=x(logical(fi),logical(ti));
    %     % figure;plot(freqs(fi==1),mean(log10(x),2));title(name)
    %     xx=mean(x,2);
    %     ff=freqs(fi==1);
    %     [aa bb]=max(xx);
    %     grid_sp.(name).peak_freq = ff(bb);
    %     grid_sp.(name).ff=ff;
    %     grid_sp.(name).fdata=xx;

    % restrict to N2 sleep alone
    I_sleep = sleep_staging==2;
    sp_st_N2=[];sp_end_N2=[];sp_epochs_N2=[];so_epochs_N2=[];
    sp_sp_epochs_N2=[];
    for ii=1:length(sp_st)
        if I_sleep(sp_st(ii)) == 1
            sp_st_N2 = [sp_st_N2 ii];
            sp_end_N2 = [sp_end_N2 ii];
            sp_epochs_N2 = [sp_epochs_N2; sp_epochs(ii,:)];
            so_epochs_N2 = [so_epochs_N2;so_epochs(ii,:)];
            sp_sp_epochs_N2 = [sp_sp_epochs_N2; sp_sp_epochs(ii,:)];
        end
    end
    grid_sp.(name).sp_st_N2 = sp_st_N2;
    grid_sp.(name).sp_end_N2 = sp_end_N2;
    grid_sp.(name).sp_epochs_N2 = sp_epochs_N2;
    grid_sp.(name).so_epochs_N2 = so_epochs_N2;
    grid_sp.(name).sp_sp_epochs_N2 = sp_sp_epochs_N2;

    % look at SO amplitude in Stage 3 sleep
    I_sleep = sleep_staging==3;
    sp_st_N3=[];sp_end_N3=[];sp_epochs_N3=[];so_epochs_N3=[];
    sp_sp_epochs_N3=[];
    for ii=1:length(sp_st)
        if I_sleep(sp_st(ii)) == 1
            sp_st_N3 = [sp_st_N3 ii];
            sp_end_N3 = [sp_end_N3 ii];
            sp_epochs_N3 = [sp_epochs_N3; sp_epochs(ii,:)];
            so_epochs_N3 = [so_epochs_N3;so_epochs(ii,:)];
            sp_sp_epochs_N3 = [sp_sp_epochs_N3; sp_sp_epochs(ii,:)];
        end
    end
    grid_sp.(name).sp_st_N3 = sp_st_N3;
    grid_sp.(name).sp_end_N3 = sp_end_N3;
    grid_sp.(name).sp_epochs_N3 = sp_epochs_N3;
    grid_sp.(name).so_epochs_N3 = so_epochs_N3;
    grid_sp.(name).sp_sp_epochs_N3 = sp_sp_epochs_N3;



end


end