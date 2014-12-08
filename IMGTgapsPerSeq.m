function sequence = IMGTgapsPerSeq(sequence,vname)



mat = fastaread('HHV-functional-gapposition.fasta');

pp = strfind(vname,'|');
if ~isempty(pp)
    vname = vname(1:pp(1)-1);
end
pp = strfind(vname,'/');
if ~isempty(pp)
    vname = vname(1:pp(1)-1);
end

found = 0;
for j = 1:size(mat,1)
    if strcmp(mat(j).Header,vname)==1            
        found = j;         
    end
end

if found~=0
    gaps_position = strfind(mat(found).Sequence, '.');%-(start(2,1)-start(1,1));
    for k = 1:length(gaps_position)
        if gaps_position(k)>0
            sequence = [sequence(1:gaps_position(k)-1), '-', sequence(gaps_position(k):end)];
        end
    end    
end

%gap_sequence = sequence;