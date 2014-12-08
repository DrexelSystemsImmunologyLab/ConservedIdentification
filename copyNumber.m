function copy_number = copyNumber(file_name)

inFileName = [strrep(file_name,'.fasta',''),'/tempSeq.fasta'];
vlength = textread([strrep(file_name,'.fasta',''),'/Vties.txt'],'%*d %*s %*d %d');
seq = fastaread(file_name);
mat = fastaread(inFileName);

order = zeros(size(mat));
sequence = cell(size(mat));
for i = 1:size(mat,1)
    order(i) = str2double(mat(i).Header);
    sequence{i,1} = mat(i).Sequence;
end

[~,ia,ic] = unique(sequence,'stable');
copy_number = zeros(size(seq,1),2);

for i = 1:size(ia,1)
    p = find(ic==i);
    temp = vlength(p);
    max_p = find(temp==max(temp)); 
    max_p = max_p(1);
    copy_number(order(p(max_p)),1) = size(p,1);
    copy_number(order(p(max_p)),2) = order(p(max_p));
    for j = 1:size(p,1)
        if j~=max_p
            copy_number(order(p(j)),2) = order(p(max_p));
        end
    end    
end

