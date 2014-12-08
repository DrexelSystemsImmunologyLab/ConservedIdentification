function [germlines,vseqs] = makeTempSeq(file_name)

seq = fastaread(file_name);
mat1 = fastaread('HHV.fasta');
mat2 = fastaread('HHJ.fasta');

germlines = cell(size(seq));
vseqs = cell(size(seq));
%jseqs = cell(size(seq));

[rev,vp,~] = textread([strrep(file_name,'.fasta',''),'/number.txt'],'%*d %d %d %d');
inFileName = [strrep(file_name,'.fasta',''),'/Vties.txt'];
fid = fopen(inFileName,'r');
dataArray = textscan(fid,'%d%s%d%d','delimiter',' ');
order = dataArray{1,1};
vgene = dataArray{1,2};
clear dataArray;
fclose(fid);
jgroup = textread([strrep(file_name,'.fasta',''),'/Jgroup.txt'],'%*d %d');
jgene = cell(size(jgroup));
for i = 1:size(jgene,1)
    if jgroup(i,1)==1
        jgene{i,1} = 'IGHJ1/4/5';
    elseif jgroup(i,1)==2
        jgene{i,1} = 'IGHJ2';
    elseif jgroup(i,1)==3
        jgene{i,1} = 'IGHJ3';
    elseif jgroup(i,1)==4
        jgene{i,1} = 'IGHJ6';
    end
end
fid = fopen([strrep(file_name,'.fasta',''),'/CDR3.txt'],'r');
formatSpec = '%*d%s%*s';
CDR3 = textscan(fid,formatSpec,'delimiter',' ');
CDR3 = CDR3{1};
fclose(fid);

outFileName = [strrep(file_name,'.fasta',''),'/tempSeq.fasta'];
fid = fopen(outFileName,'w');

for i = 1:size(vgene,1)
    vname = char(vgene{i,1});
    bar = strfind(vname,'|');
    if ~isempty(bar)
        temp_vname = vname(1:bar(1)-1);
    else
        temp_vname = vname;
    end
    for j = 1:size(mat1,1)
        if strcmp(mat1(j).Header,temp_vname)==1
            if j==17
                vseq = mat1(j).Sequence(1:strfind(mat1(j).Sequence,'TATCACTGT')+5);
            elseif j==42
                vseq = mat1(j).Sequence(1:strfind(mat1(j).Sequence,'TATTGCTGT')+5);
            else
                vseq = mat1(j).Sequence(1:strfind(mat1(j).Sequence,'TATTACTGT')+5);
            end
        end
    end
            
    jname = char(jgene{i,1});
    if strcmp(jname,'IGHJ1/4/5')==1
        jname = 'IGHJ1';
    end
    for j = 1:size(mat2,1)
        if strcmp(mat2(j).Header,jname)==1
            jseq = mat2(j).Sequence;
        end
    end    
    
    germline = [vseq,char(CDR3{i,1}),jseq];
    
    junction_length = length(char(CDR3{i,1}));
    length_germline = length(vseq) + junction_length + length(jseq);
    
    if rev(i)==0
        sequence = seq(order(i)).Sequence;
    else
        sequence = seqrcomplement(seq(order(i)).Sequence);                    
    end           
    
    temp = length(vseq)-5-vp(i);
    if temp>0
        nSequence = [vseq(1:temp),sequence];
    else
        nSequence = sequence(-1*temp+1:end);
    end
    if length(nSequence)>=length_germline
        nSequence = nSequence(1:end-(length(nSequence)-length_germline));
    else
        temp = length_germline-length(nSequence);
        nSequence = [nSequence,jseq(end-temp+1:end)];
    end

    p = strfind(nSequence,'N');
    nSequence(p) = germline(p);

    fprintf(fid,'>%d\r\n',order(i));
    fprintf(fid,'%s\r\n',nSequence);
    germlines{order(i),1} = germline;
    vseqs{order(i),1} = vseq;
    %jseqs{order(i),1} = jseq;
end