
data=EEG.data;

% plot the stim trigger -> the channel differs by recording?
figure;
stim = data(83,:);
%stim = data(79,:);
stem(stim)

% remove unncesary triggers at the end
% [aa bb]=ginput;
% aa=round(aa);
% vline(aa,'r')
% stim(aa:end)=0;

% extract the stim onsets from the triggers
d = [0 diff(stim)];
d1 = (find(d==1)); % these are all the transitions from 0 to 1
d1 = [d1 d1(end)+600];
% identifying stim onsets
stim_onsets=[];
done=false;k=1;
while ~done

    % now go forward till you have atleast 500 consecutive zeros.
    done1=false;
    init_k = k;
    while ~done1
        if (d1(k+1)-d1(k) < 500) 
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
stim_onsets = stim_onsets(2:3:end); % taking into account just the two pulse flashes
