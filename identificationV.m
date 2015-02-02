function [vgene, vposition, mismatch, aligned_length] = identificationV(sequence)
% finds V anchors in sequence, outputs J group and position of J in sequence

vgene = '';
vposition = 0;
mismatch = 0;
aligned_length = 0;
sequence = upper(sequence);

%% Look for DXXXYYC
pattern = 'GAC.........TATTACTGT';
p = regexp(sequence,pattern);
if ~isempty(p)
    vposition = p + 12;     % change vposition to 'TATTACTGT'
end

%% Look for YYC
if vposition==0
    pattern = 'TATTACTGT';
    p = strfind(sequence,pattern);
    if ~isempty(p)
        if length(p)>1      % Multiple anchors
            d = zeros(size(p));
            for i = 1:length(p)
                if length(sequence)>p+21
                    d(i) = sum(sequence(p-12:p-10)~='GAC');
                end
            end
            [~,idx] = min(d);     % the one with least mutation in D
            vposition = p(idx);
        else
            vposition = p;
        end    
    end
end

%% Look for DXXXXXC
if vposition==0
    pattern = 'GAC...............TGT';
    p = regexp(sequence,pattern);
    if ~isempty(p)
        if length(p)>1      % Multiple anchors 
            d = zeros(size(p));
            for i = 1:length(p)
                if length(sequence)>p+17
                    d(i) = sum(sequence(p+12:p+17)~='TATTAC');
                end
            end
            [~,idx] = min(d);     % the one with least mutation in YY
            vposition = p(idx) + 12;     % change vposition to 'TATTACTGT'
        else
            vposition = p + 12;     % change vposition to 'TATTACTGT'
        end
    end
end

%% compare germlines, output germline(s) with least mismatches
if vposition~=0
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
end

