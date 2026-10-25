function so_synchrony = build_so_synchrony(grid_so,a,b)

num_ch = a;
len = b;

so_synchrony = zeros(num_ch,len);

for j=1:num_ch
    name = ['ch' num2str(j)];
    tmp = grid_so.(name).so_st;
    tmp1 = grid_so.(name).so_end;
    for k=1:length(tmp)
        so_synchrony(j,tmp(k):tmp1(k)) = 1;
    end
end


end