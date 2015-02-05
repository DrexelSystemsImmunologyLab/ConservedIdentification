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
        vposition = p + 12;     % change vposition to 'TATTACTGT'        
    end
end

%% Compare with germlines
vgene = cell(size(vposition));
mismatch = zeros(size(vposition));
aligned_length = zeros(size(vposition));
for i = 1:length(vposition)
    [vgene{i}, mismatch(i), aligned_length(i)] = compareV(sequence,vposition(i));
end

if length(vgene)>1      % if multiple anchors
    freq = mismatch./aligned_length;
    [~,idx] = min(freq);     % choose the one with least mutation frequency
    vgene = char(vgene{idx});
    mismatch =mismatch(idx);
    aligned_length = aligned_length(idx);
else
    vgene = char(vgene);
end
