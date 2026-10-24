function [grid_so] = ...
    SO_analyses_STAR(data_main,I,soFilt,spFilt1,spFilt2,sleep_staging,Fs,frontal_elec,posterior_elec)
%function [grid_so] = SO_analyses(data,I,bpFilt,spFilt1,spFilt2,sleep_staging)

% identify the ied times large high freq spikes
I=logical(I);
II=NaN(size(data_main));
win = Fs*2; % discarding 2s of data either side of an artifact
parfor i=1:size(data_main,1)
    % artifact threshold
    tmp = (I);
    xx=find(abs(data_main(i,:))>250);
    for j=1:length(xx)
        if xx(j)+win < length(I) && xx(j)>win+1
            tmp([-win:win]+xx(j))=0;
        elseif xx(j)+win>length(I)
            tmp(xx(j):end)=0;
        end
    end
    II(i,:) = tmp;
end



grid_so = struct;
for i=1:size(data_main,1)

    disp(['Processing SOs at Ch ' num2str(i)]);
    I=II(i,:);
    name = ['ch' num2str(i)];

    if sum(abs(data_main(i,:))) > 0
        so_data = data_main(i,:);
        so_data = filtfilt(soFilt,so_data);

        if ismember(i,frontal_elec)
            sp_data = abs(hilbert(filtfilt(spFilt1,data_main(i,:))));
        elseif ismember(i,posterior_elec)
            sp_data = abs(hilbert(filtfilt(spFilt2,data_main(i,:))));
        end

        [so_st_new so_end_new] = SO_detect_STAR(-so_data,[],[],I,Fs);
        %[so_st_new so_end_new] = combine_events(so_st_new,so_end_new,500,4000);

        % get the epochs
        ep=[];ep_raw=[];
        for ii=1:length(so_st_new)
            temp=so_data(so_st_new(ii):so_end_new(ii));
            [aa bb]=max(temp);
            bb = bb + so_st_new(ii)-1;
            ep(ii,:) = so_data(bb+[-600:600]);
            ep_raw(ii,:) = data_main(i,bb+[-600:600]);
        end
        ep=(ep')';

        % storing
        grid_so.(name).so_st = so_st_new;
        grid_so.(name).so_end = so_end_new;
        grid_so.(name).so_ep = ep;
        grid_so.(name).ep_raw = ep_raw;

        % get the stage 3 background Spindle power
        % plot the phase-locked SP
        I_sleep = sleep_staging==3;
        ep_sp=[];ep_so=[];N3_so_idx=[];
        for ii=1:length(so_st_new)
            if I_sleep(so_st_new(ii)) == 1
                N3_so_idx = [N3_so_idx ii];
                temp=sp_data(so_st_new(ii):so_end_new(ii));
                [aa bb]=max(temp);
                bb = bb + so_st_new(ii)-1;
                ep_sp = [ep_sp; sp_data(bb+[-600:600])];

                temp=so_data(so_st_new(ii):so_end_new(ii));
                [aa bb]=max(temp);
                bb = bb + so_st_new(ii)-1;
                ep_so = [ep_so; so_data(bb+[-600:600])];
            end
        end

        grid_so.(name).ep_sp_N3 = ep_sp;
        grid_so.(name).ep_so_N3 = ep_so;
        grid_so.(name).N3_so_idx = N3_so_idx;

    else
        grid_so.(name).so=[];
    end
end
