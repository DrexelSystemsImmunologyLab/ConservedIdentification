function [vgene, vposition, mismatch, aligned_length, indel] = identificationV(sequence)
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
        vposition = p + 12;     % change vposition to 'TATTACTGT'        
    end
end

%% Compare with germlines
if vposition(1)~=0
    vgene = cell(size(vposition));
    mismatch = zeros(size(vposition));
    aligned_length = zeros(size(vposition));
    indel = zeros(size(vposition));
    for i = 1:length(vposition)
        [vgene{i}, mismatch(i), aligned_length(i), indel(i)] = compareV(sequence,vposition(i));
    end

    if length(vgene)>1      % if multiple anchors
        freq = mismatch./aligned_length;
        [~,idx] = min(freq);     % choose the one with least mutation frequency
        if length(idx)>1
            idx = idx(aligned_length==max(aligned_length));     % if equal frequency, choose the one with longest alignment
        end
        idx = idx(1);
        vgene = char(vgene{idx});
        vposition = vposition(idx);
        mismatch =mismatch(idx);
        aligned_length = aligned_length(idx);
        indel = indel(idx);
    else
        vgene = char(vgene);
    end
end
