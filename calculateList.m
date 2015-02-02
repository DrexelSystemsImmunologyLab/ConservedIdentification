function mix_list = calculateList(germ_db, vlength, mut_fre,prob_threshold)
%germ_db: file name for germline database
%vlength: length of aligned V region
%mut_fre: mutation frequency
%prob_threshold: gene and allele pairs with mixing probability greater than
%this will be labeled as mixing pairs

%% read germlines
mat = fastaread(germ_db);

%% calulate length of V genes
l = zeros(size(mat));
for i = 1:size(mat,1)
    l(i,1) = length(mat(i).Sequence);    
end

%% trim vlength
vlength = round(vlength);
if vlength>min(l)
    vlength = min(l);
end

%% create difference matrix
d_matrix = zeros(size(mat,1),size(mat,1));
for i = 1:size(mat,1)
    for j = i+1:size(mat,1)
        string1 = mat(i).Sequence(end-vlength+1:end);
        string2 = mat(j).Sequence(end-vlength+1:end);        
        d_matrix(i,j) = sum(string1~=string2);
    end
end

%% calculate prob. of mixing
M = vlength;
N = ceil(M*mut_fre);
p_matrix = zeros(size(d_matrix));
for i = 1:size(p_matrix,1)
    for j = i+1:size(p_matrix)
        K = d_matrix(i,j);        
        p = 0;
        for k = ceil(K/2):K
            p = p + hygepdf(k,M,K,N)*0.33^k;
        end
        p_matrix(i,j) = p;
    end
end

%% find pairs with prob. greater than threshold
[row,col] = find(p_matrix>prob_threshold);

%% output list
mix_list = cell(length(row),2);
for i = 1:size(row,1)
    mix_list{i,1} = mat(row(i)).Header;
    mix_list{i,2} = mat(col(i)).Header;
end

%% print output
outfilename = ['Vties_',strrep(germ_db,'.fasta',''),'_',num2str(vlength),'_',num2str(mut_fre),'_',num2str(prob_threshold),'.txt'];
fid = fopen(outfilename,'w');
for i = 1:size(mix_list,1)
    fprintf(fid,'%s|%s\r',char(mix_list{i,1}),char(mix_list{i,2}));
end
fclose(fid);
