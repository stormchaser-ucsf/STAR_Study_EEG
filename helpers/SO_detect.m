function [so_st_new, so_end_new, k] = SO_detect(temp11,sleep_st,sleep_end,II,Fs)
%
%  [so_st_new so_end_new] = SO_detect(temp11,sleep_st,sleep_end,II)
%SO_detect




ds=temp11;

% find all positive to negative zero crossing.
zero_cross=[];
for i=2:length(ds)
   % disp(i/length(ds))
    if ds(i) <0 && ds(i-1) >=0 
        zero_cross = [zero_cross; i];
    end
end

%I=(zero_cross > sleep_st) .* (zero_cross < sleep_end);
III = zeros(1,length(ds));
III(zero_cross)=1;
zero_cross1 = III .* II;
%zero_cross = zero_cross(logical(II));
zero_cross = find(zero_cross1==1);

% subsequent positive to negative zero crossing between 0.9 and 3s
so_st=[];
so_end=[];
for i=2:length(zero_cross)
    d= zero_cross(i) - zero_cross(i-1);
    if (d*(1000/Fs)) > 900 && (d*(1000/Fs)) < 3000
        so_st = [so_st; zero_cross(i-1)];
        so_end = [so_end; zero_cross(i)];
    end
end

% amplitude difference between the peak and trough
so_amp=[];
for i=1:length(so_st)
    temp=ds(so_st(i):so_end(i));
    % find the zero crossing points
    for j=2:length(temp)
        if temp(j)>0 && temp(j-1) <=0
            break
        end
    end
    aa1 = max(abs(temp(1:j-1)));
    aa2 = max(abs(temp(j:end)));   
  % IMPORTANT POINT: IF THE BELOW IF CODE IS GOING TO BE UNCOMMENTED, THEN
  % PLEAST ALSO UNCOMMENT THE LINK FOLLOWING THE ELSE STATEMENT ASSIGNING
  % ZERO IF CONDITION NOT SATISFIED
    %  if aa1 > 5
        so_amp=[so_amp;aa1+aa2]; 
   % else
    %    so_amp=[so_amp;0]; 
    %end
end


% Ampliude thresholds to get atleast 100 SO
so_amp_hist=sort(so_amp(so_amp>0),'ascend');
so_st1=[];so_end1=[];
k=0.96;
so_st_new=[];
so_end_new=[];
tracking=[];
while k>=0.90 %% && length(so_st_new)<100
    tracking= [tracking k];
    k=k-0.005;
    so_st1=[];so_end1=[];
    so_amp_thresh = so_amp_hist(ceil(k*length(so_amp_hist)));
    IL= so_amp >= so_amp_thresh;
    so_st1=so_st(IL);
    so_end1=so_end(IL);
    
    so_st_new=[];
    so_end_new=[];
    for i=1:length(so_st1)
        temp=(ds(so_st1(i):so_end1(i)));
        j=1;
        while temp(j)<=0
            j=j+1;
        end
        temp1=temp(1:j-1);
        temp2=diff(temp1);
        zc=[];
        for j=2:length(temp2)
            if temp2(j) >0 && temp2(j-1)<0
                zc = [zc;j-1];
            end
        end
        if length(zc)==1
            so_st_new=[so_st_new;so_st1(i)];
            so_end_new=[so_end_new;so_end1(i)];
        end
    end
    
end

