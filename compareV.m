function [vgene, mismatch, aligned_length, indel] = compareV(sequence,vposition)
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

indels = zeros(size(p));      % checking for potential indels
if mismatch>=18
    for i = 1:size(p,1)
        for j = 1:size(mutation,2)-29
            if sum(mutation(p(i),j:j+29)~=0)>=18      % 18 mutations in a sliding window of 30
                indels(i,1) = 1;      % flag as potential indels
            end
        end
    end
end

if sum(indels==0)/length(indels)>0.5
    vgene = Vname(p(indels==0));      % remove V genes that are potential indels
    indel = 0;
else
    vgene = Vname(p(indels==0));      % keep all the V genes
    indel = 1;      % label the sequence as indel
end
str = strjoin(vgene','|');      % converting cell to string
vgene = str;