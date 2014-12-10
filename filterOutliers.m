function filterOutliers

[idx,vgene,mutation,coverage] = textread('Vties.txt','%d %s %d %d');
fid = fopen('outliers.txt');
for i = 1:size(idx,1)
    if mutation(i,1)/coverage(i,1)>0.6
        fprintf(fid,'%d %s %d %d',idx(i,1),char(vgene{i,1}),mutation(i,1),coverage(i,1));
    end
end
fclose(fid);