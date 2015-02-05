function [vgene, mismatch, aligned_length] = compareV(sequence,vposition)
% compare germlines, output germline(s) with least mismatches

load('germlines.mat');      % germline sequence data
sequence = sequence(1:vposition);
sbinary = zeros(size(sequence));
sbinary(sequence=='A') = 1;     % change string to doubles
sbinary(sequence=='C') = 4;
sbinary(sequence=='G') = 9;
sbinary(sequence=='T') = 16;
aligned_length = length(sbinary);       % number of aligned nucleotides
m1 = zeros(size(Vbinary,1),aligned_length);
m2 = zeros(size(Vbinary,1),aligned_length);
for i = 1:size(Vbinary,1)
    m1(i,:) = Vbinary{i,1}(end-aligned_length+1:end);    % germline
    m2(i,:) = sbinary;      % sequence
end
mutation = m1 - m2;
mcount = (mutation~=0);
msum = sum(mcount,2);
mismatch = min(msum);       % least number of mismatches
p = find(msum==mismatch);       % alleles with least mismatches
vgene = Vname(p);
str = strjoin(vgene','|');      % converting cell to string
vgene = str;