function [topo_spindle_freq,topo_spindle_freq_N3] = ...
    topo_spindle_frequency(data,sleep_staging,onebyf_flag)
%function topo_spindle_freq = topo_spindle_frequency(data,onebyf_flag)
% INPUT: 
% data
% onebyf_flag is 1 or 0



% get center frequency of spindle per channel using fft
topo_spindle_freq=[];
topo_spindle_freq_N3=[];
win_length=4096;
if onebyf_flag==0
    onebyf_flag = false;
else
    onebyf_flag = true;
end
Fs=1e3;

parfor ch=1:64
    disp(['Processing Channel ' num2str(ch)]);
    chdata = data(ch,:);
    seg = 1:win_length:length(chdata);
    ch_spindle_freq=NaN(1,length(seg));
    ch_spindle_freq_N3=NaN(1,length(seg));
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

            if onebyf_flag
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
            end

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

             % store if in N3 sleep
            if sum(sleep_seg==3) > 0
                ch_spindle_freq_N3(i) = pref_freq;
            end

        end

    end

    % figure;hist(ch_spindle_freq)
    % figure;ksdensity(ch_spindle_freq)
    % nanmean(ch_spindle_freq)

    topo_spindle_freq(ch) = nanmedian(ch_spindle_freq);
    topo_spindle_freq_N3(ch) = nanmedian(ch_spindle_freq_N3);

end

%figure;plot(topo_spindle_freq)

% plot on brain
% fast_spindle_idx = topo_spindle_freq>=median(topo_spindle_freq);
% slow_spindle_idx = topo_spindle_freq<median(topo_spindle_freq);
% figure;
% topoplot(topo_spindle_freq, EEG.chanlocs); 
% axis tight
% clim([min(topo_spindle_freq) max(topo_spindle_freq)])
% disp([min(topo_spindle_freq) max(topo_spindle_freq)])
% colorbar
% if onebyf_flag
%     title("Preferred Spindle Freq 1/f")
% else
%     title("Preferred Spindle Freq no 1/f yes eye blink")
% end
