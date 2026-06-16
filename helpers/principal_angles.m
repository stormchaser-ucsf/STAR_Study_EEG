function [angles,Va,Vb,pcap] = principal_angles(Wa,Wb)
%  [angles,Va,Vb,pcap] = principal_angles(Wa,Wb)
%
% COMPUTES PRINCIPAL ANGLES BETWEEN TWO MANIFOLDS WA AND WB
% INPUT
% Wa:   N X p matrix, with N channels/neurons/measurements and p
%       dimensions
% Wb:   N X q matrix, with N channels/neurons/measurements and q
%       dimensions. q can be different from p
% OUTPUT
% angles: Principal angles between Wa and Wb
% Va:   Directions in the basis of Wa that give smallest angle with Wb
% Vb:   Directions in the basis of Wb that give smallest angles with Wa
% pcap: Percent variance captured between the two manifolds
% N.N. Nov 2020


% check if orthonormal basis
Wa_norm = Wa'*Wa;
check = Wa_norm-eye(size(Wa_norm));
Wb_norm = Wb'*Wb;
check1 = Wb_norm-eye(size(Wb_norm));

if norm(check,'fro')<1e-10 && norm(check1,'fro')<1e-10 % they are orthogonal
    [angles,Va,Vb,pcap] = find_angles(Wa,Wb);
else
    % do QR decomposition if not orthogonal
    [qa1,ra1]=qr(Wa,0);
    [qa2,ra2]=qr(Wb,0);
    [angles,Va,Vb,pcap] = find_angles(qa1,qa2);
end

function [out,out_Va,out_Vb,out_pcap] = find_angles(in1,in2)
    prod_mat = in1'*in2;
    [out_Va,s,out_Vb]=svd(prod_mat); % use the SVD method
    out=acosd(diag(s));

    % now get the VAF captured betwen the two
    tmp=zeros(2,1);
    num = (in1*in1')*(in2*in2')*(in1*in1');
    den = in2*in2';
    tmp(1) = trace(num)/trace(den);

    num = (in2*in2')*(in1*in1')*(in2*in2');
    den = in1*in1';
    tmp(2) = trace(num)/trace(den);
    out_pcap = mean(tmp);
    
end

end