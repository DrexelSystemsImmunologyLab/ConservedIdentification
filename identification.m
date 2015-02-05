function identification(file_name)
% v gene and j gene identification

tic;

%% make new folder
folder_name = strrep(file_name,'.fasta','');
mkdir(folder_name);
disp('making new directory done');

%% read fasta file
seq = fastaread(file_name);
s = [num2str(size(seq,1)),' sequence(s) loading complete'];
disp(s);

%% make output files
outFileName1 = [folder_name,'/noFoundJ.txt'];     % order number of sequences without J
fid1 = fopen(outFileName1,'w');
outFileName2 = [folder_name,'/noFoundV.txt'];     % order number of sequences with J but without V
fid2 = fopen(outFileName2,'w');
outFileName3 = [folder_name,'/intactCDR3.txt'];    % order number of sequences with both J and V
fid3 = fopen(outFileName3,'w');
outFileName4 = [folder_name,'/Jgene.txt'];     % order number and identified J gene groups
fid4 = fopen(outFileName4,'w');
outFileName5 = [folder_name,'/Vgene.txt'];     % order number and identified V genes, number of mismatches and number of aligned nucleotides
fid5 = fopen(outFileName5,'w');
outFileName6 = [folder_name,'/CDR3.txt'];     % order number and CDR3 nucleotides and amino acids
fid6 = fopen(outFileName6,'w');
outFileName7 = [folder_name,'/numbers.txt'];     % order number and positions of J and V anchor
fid7 = fopen(outFileName7,'w');

%% identify V and J gene(s)
for i = 1:size(seq,1)
    sequence = seq(i).Sequence;
    [jgroup, jposition] = identificationJ(sequence);      % identify J
    if jgroup==0        % not finding J
        sequence = seqrcomplement(sequence);        % reverse compliment
        [jgroup, jposition] = identificationJ(sequence);      % identify J in reverse compliment
    end
    if jgroup==0        % not finding J in reverse compliment
        fprintf(fid1,'%d\r',i);     % put order number into 'noFoundJ.txt'
    end
    if jgroup~=0 && jposition~=0     % J is found
        [vgene, vposition, mismatch, aligned_length] = identificationV(sequence);    % identify V to the left of J
        if isempty(vgene)      % not finding V
            fprintf(fid2,'%d\r',i);     % put order number into 'noFoundV.txt'
        elseif ~isempty(vgene);     % V is found, thus CDR3 is intact
            fprintf(fid3,'%d\r',i);     % put order number into 'intactCDR3.txt'
            fprintf(fid4,'%d\t%d\r',i,jgroup);     % put order number and J gene group into 'Jgene.txt'
            fprintf(fid5,'%d\t%s\t%d\t%d\r',i,vgene,mismatch,aligned_length);     % put order number and V genes, number of mismatches and number of aligned nucleotides into 'Vgene.txt'
            if mod(jposition - vposition - 19,3)==0
                fprintf(fid6,'%d\t%s\t%s\r',i,sequence(vposition+6:jposition-14),nt2aa(sequence(vposition+6:jposition-14),'ACGTonly',false));     % put order number and in-frame CDR3 nucleotides and amino acids into 'CDR3.txt'
            else
                fprintf(fid6,'%d\t%s\t\r',i,sequence(vposition+6:jposition-14));     % put order number and out-of-frame CDR3 nucleotides into 'CDR3.txt'
            end  
            fprintf(fid7,'%d\t%d\t%d\r',i,vposition,jposition);     % order number and positions of J and V anchor in 'number.txt'
        end
    end
end
fclose('all');
disp('identification done');

%% calculate v ties
inFileName = [folder_name,'/Vgene.txt'];     % read V identification stats
fid = fopen(inFileName,'r');
dataArray = textscan(fid,'%d%s%d%d','delimiter','\t');
fclose(fid);
vlength = mean(dataArray{4});
mut_freq = mean(dataArray{3})/vlength;
mix_list = calculateList('HHV.fasta', vlength, mut_freq,0.01);

%% apply v ties
outFileName = [folder_name,'/Vties.txt'];
fid = fopen(outFileName,'w');
vgenes = dataArray{2};
for i = 1:size(vgenes,1)
    str = char(vgenes{i});
    c = strsplit(str,'|')';    
    idx = ismember(mix_list(:,1),c);
    c = union(c,mix_list(idx,2));       % apply v ties
    for j = 1:length(c)
        p = strsplit(char(c{j}),'*');
        c{j} = char(p{1});     % only keep gene information, discard alleles 
    end
    c = unique(c);
    if size(c,1)==1
        str = strjoin(c,'|');
    elseif size(c,2)==1
        str = strjoin(c','|');
    end
    fprintf(fid,'%d\t%s\t%d\t%d\r',dataArray{1}(i),str,dataArray{3}(i),dataArray{4}(i));     % output information into 'Vties.txt'
end
fclose(fid);
disp('applying v ties done');

toc;
