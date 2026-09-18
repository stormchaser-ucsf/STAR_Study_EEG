function [start_filtered, end_filtered] = ...
    combine_events(start_times,end_times,separation_threshold,duration_threshold,fs)

% separation_threshold and duration_threshold are in milliseconds




start_filtered = [];
end_filtered = [];

i = 0;

Fss=1e3; %have to normalize to get in milliseconds conversion

while i < length(start_times)
    i = i + 1;  
    j = i;
    tempStart = start_times(i);
    
    if (j+1) < length(start_times)
        separation = start_times(j+1) - end_times(j);
        duration = end_times(j+1) - tempStart;
        
        while separation*(Fss/fs) <= separation_threshold && duration*(Fss/fs) <= duration_threshold
            if (j+1) < length(start_times)
            j = j + 1;
            else 
                break;
            end
            separation = start_times(j+1) - end_times(j);
            duration = end_times(j+1) - tempStart;
            if separation*(Fss/fs) > separation_threshold || duration*(Fss/fs) > duration_threshold
                i = j; break;
            end
        end
        
        tempEnd = end_times(j);
        start_filtered  = [start_filtered tempStart];
        end_filtered = [end_filtered tempEnd];
        
    else
        start_filtered = [start_filtered start_times(i)];
        end_filtered = [end_filtered end_times(j)];
    end
    
end

end






