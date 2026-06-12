function [prin_angles] = compute_prin_angles_PC(dataTensor,k)
%function [prin_angles] = compute_prin_angles_PC(dataTensor)
%dataTensor -> time X channels X conditions
prin_angles=[];
for i=1:size(dataTensor,3)
    Xa = (squeeze(dataTensor(:,:,i)));    
    [Wai,s1,l1]=pca(Xa);    
    for j=i+1:size(dataTensor,3)
        Xa = (squeeze(dataTensor(:,:,j)));
        [Waj,s2,l2]=pca(Xa);    
        temp = subspacea(Wai(:,1:k),Waj(:,1:k))*180/pi;
        prin_angles = [prin_angles temp];        
    end
end
end