function conservedIdentification(file_name,seq)

seq = fastaread(file_name);
s = [num2str(size(seq,1)),' sequence loading complete'];
disp(s);

new_file_name1 = strrep(file_name,'.fasta','/number.txt');
fid1 = fopen(new_file_name1,'w');
new_file_name2 = strrep(file_name,'.fasta','/CDR3.txt');
fid2 = fopen(new_file_name2,'w');
new_file_name3 = strrep(file_name,'.fasta','/Jgroup.txt');
fid3 = fopen(new_file_name3,'w');
new_file_name4 = strrep(file_name,'.fasta','/intactCDR.fasta');
fid4 = fopen(new_file_name4,'w');
new_file_name5 = strrep(file_name,'.fasta','/Vgene.txt');
fid5 = fopen(new_file_name5,'w');
new_file_name6 = strrep(file_name,'.fasta','/Vallele.txt');
fid6 = fopen(new_file_name6,'w');
new_file_name7 = strrep(file_name,'.fasta','/MultipleYYC.fasta');
fid7 = fopen(new_file_name7,'w');
new_file_name8 = strrep(file_name,'.fasta','/Jonly.fasta');
fid8 = fopen(new_file_name8,'w');
new_file_name9 = strrep(file_name,'.fasta','/Nofound.fasta');
fid9 = fopen(new_file_name9,'w');

load('germlines.mat');

for i = 1:size(seq,1)
    %if mod(i,floor(size(seq,1)/20))==0        
    %    a = floor(i/floor(size(seq,1)/20));
    %    if a>1
    %        fprintf(repmat('\b',1,21));
    %    end
    %    s = [repmat('=',1,a),repmat('-',1,20-a)];
    %    fprintf('%s\r',s);       
    %end
    Vfound = 0;
    Jfound = 0;
    rev = 0;
    %%%%%%%%%%%%%%%%%%%%%%%
    sequence = strrep(seq(i).Sequence,'@','A');
    %%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(strfind(sequence,'TGGTCACCGTCTCCTCAG'))
        Jfound = 1;
        Jposition = strfind(sequence,'TGGTCACCGTCTCCTCAG');
        Jgene = 1;    %IGHJ1/4/5
    end
    if ~isempty(strfind(sequence,'TGGTCACTGTCTCCTCAG'))
        Jfound = 1;
        Jposition = strfind(sequence,'TGGTCACTGTCTCCTCAG');
        Jgene = 2;    %IGHJ2
    end
    if ~isempty(strfind(sequence,'TGGTCACCGTCTCTTCAG'))
        Jfound = 1;
        Jposition = strfind(sequence,'TGGTCACCGTCTCTTCAG');
        Jgene = 3;    %IGHJ3
    end
    if ~isempty(strfind(sequence,'CGGTCACCGTCTCCTCAG'))
        Jfound = 1;
        Jposition = strfind(sequence,'CGGTCACCGTCTCCTCAG');
        Jgene = 4;    %IGHJ6
    end
    if Jfound==0
        if ~isempty(strfind(sequence,'TGGTCACCGTCTCCT'))
            Jfound = 1;
            Jposition = strfind(sequence,'TGGTCACCGTCTCCT');
            Jgene = 1;
        end
        if ~isempty(strfind(sequence,'TGGTCACTGTCTCCT'))
            Jfound = 1;
            Jposition = strfind(sequence,'TGGTCACTGTCTCCT');
            Jgene = 2;
        end
        if ~isempty(strfind(sequence,'TGGTCACCGTCTCTT'))
            Jfound = 1;
            Jposition = strfind(sequence,'TGGTCACCGTCTCTT');
            Jgene = 3;
        end
        if ~isempty(strfind(sequence,'CGGTCACCGTCTCCT'))
            Jfound = 1;
            Jposition = strfind(sequence,'CGGTCACCGTCTCCT');
            Jgene = 4;
        end
    end
    if Jfound==0
        if ~isempty(strfind(sequence,'TGGTCACCGTCT'))
            Jfound = 1;
            Jposition = strfind(sequence,'TGGTCACCGTCT');
            Jgene = 1;
        end
        if ~isempty(strfind(sequence,'TGGTCACTGTCT'))
            Jfound = 1;
            Jposition = strfind(sequence,'TGGTCACTGTCT');
            Jgene = 2;
        end
        if ~isempty(strfind(sequence,'TGGTCACCGTCT'))
            Jfound = 1;
            Jposition = strfind(sequence,'TGGTCACCGTCT');
            Jgene = 3;
        end
        if ~isempty(strfind(sequence,'CGGTCACCGTCT'))
            Jfound = 1;
            Jposition = strfind(sequence,'CGGTCACCGTCT');
            Jgene = 4;
        end
    end
    if Jfound==1
        if ~isempty(strfind(sequence(1:Jposition),'TATTACTGT'))
            Vfound = 1;
            Vposition = strfind(sequence(1:Jposition),'TATTACTGT');
        end
    end
    if Jfound==0 && Vfound==0
        sequence = seqrcomplement(sequence);
        rev = 1;
        if ~isempty(strfind(sequence,'TGGTCACCGTCTCCTCAG'))
            Jfound = 1;
            Jposition = strfind(sequence,'TGGTCACCGTCTCCTCAG');
            Jgene = 1;
        end
        if ~isempty(strfind(sequence,'TGGTCACTGTCTCCTCAG'))
            Jfound = 1;
            Jposition = strfind(sequence,'TGGTCACTGTCTCCTCAG');
            Jgene = 2;
        end
        if ~isempty(strfind(sequence,'TGGTCACCGTCTCTTCAG'))
            Jfound = 1;
            Jposition = strfind(sequence,'TGGTCACCGTCTCTTCAG');
            Jgene = 3;
        end
        if ~isempty(strfind(sequence,'CGGTCACCGTCTCCTCAG'))
            Jfound = 1;
            Jposition = strfind(sequence,'CGGTCACCGTCTCCTCAG');
            Jgene = 4;
        end
        if Jfound==0
            if ~isempty(strfind(sequence,'TGGTCACCGTCTCCT'))
                Jfound = 1;
                Jposition = strfind(sequence,'TGGTCACCGTCTCCT');
                Jgene = 1;
            end
            if ~isempty(strfind(sequence,'TGGTCACTGTCTCCT'))
                Jfound = 1;
                Jposition = strfind(sequence,'TGGTCACTGTCTCCT');
                Jgene = 2;
            end
            if ~isempty(strfind(sequence,'TGGTCACCGTCTCTT'))
                Jfound = 1;
                Jposition = strfind(sequence,'TGGTCACCGTCTCTT');
                Jgene = 3;
            end
            if ~isempty(strfind(sequence,'CGGTCACCGTCTCCT'))
                Jfound = 1;
                Jposition = strfind(sequence,'CGGTCACCGTCTCCT');
                Jgene = 4;
            end
        end
        if Jfound==0
            if ~isempty(strfind(sequence,'TGGTCACCGTCT'))
                Jfound = 1;
                Jposition = strfind(sequence,'TGGTCACCGTCT');
                Jgene = 1;
            end
            if ~isempty(strfind(sequence,'TGGTCACTGTCT'))
                Jfound = 1;
                Jposition = strfind(sequence,'TGGTCACTGTCT');
                Jgene = 2;
            end
            if ~isempty(strfind(sequence,'TGGTCACCGTCT'))
                Jfound = 1;
                Jposition = strfind(sequence,'TGGTCACCGTCT');
                Jgene = 3;
            end
            if ~isempty(strfind(sequence,'CGGTCACCGTCT'))
                Jfound = 1;
                Jposition = strfind(sequence,'CGGTCACCGTCT');
                Jgene = 4;
            end
        end
        if Jfound==1
            if ~isempty(strfind(sequence(1:Jposition),'TATTACTGT'))
                Vfound = 1;
                Vposition = strfind(sequence(1:Jposition),'TATTACTGT');
            end
        end        
    end
    if Vfound~=0
        if length(Vposition)>1
            Vfound = 2;
        end        
    end
    if Vfound==1 && Jfound==1
        Jposition = Jposition(end);
        fprintf(fid1,'%d %d %d %d\r\n',i,rev,Vposition,Jposition);
        if mod(Jposition - Vposition - 19,3)==0
            fprintf(fid2,'%d %s %s\r\n',i,sequence(Vposition+6:Jposition-14),nt2aa(sequence(Vposition+6:Jposition-14),'ACGTonly',false));
        else
            fprintf(fid2,'%d %s outofframe\r\n',i,sequence(Vposition+6:Jposition-14));
        end  
        fprintf(fid3,'%d %d\r\n',i,Jgene);
        fprintf(fid4,'>%d\r\n',i);
        fprintf(fid4,'%s\r\n',sequence);
        
        sequence = sequence(1:Vposition);
        sbinary = zeros(size(sequence));
        sbinary(sequence=='A') = 1;
        sbinary(sequence=='C') = 4;
        sbinary(sequence=='G') = 9;
        sbinary(sequence=='T') = 16;
        l = length(sbinary);
        m1 = zeros(size(Vbinary,1),l);
        m2 = zeros(size(Vbinary,1),l);
        for j = 1:size(Vbinary,1)
            m1(j,:) = Vbinary{j,1}(end-l+1:end);
            m2(j,:) = sbinary;
        end
        mutation = m1 - m2;
        mcount = (mutation~=0);
        msum = sum(mcount,2);
        [mmin,p] = minCI(msum);
        Vallele = char(Vname{p(1)});
        star = strfind(char(Vname{p(1)}),'*');
        Vgene = char(Vname{p(1)}(1:star-1));
        if length(p)>1
            for j = 2:length(p)
                Vallele = [Vallele,'|',char(Vname{p(j)})];
                star = strfind(char(Vname{p(j)}),'*');
                if isempty(strfind(Vgene,char(Vname{p(j)}(1:star-1))))
                    Vgene = [Vgene,'|',char(Vname{p(j)}(1:star-1))];
                end
            end
        end
        fprintf(fid5,'%d %s %d %d\r\n',i,Vgene,mmin,l);
        fprintf(fid6,'%d %s %d %d\r\n',i,Vallele,mmin,l);  
    elseif Vfound==2 && Jfound==1   %Multiple YYC
        fprintf(fid7,'>%d\r\n',i);
        fprintf(fid7,'%s\r\n',seq(i).Sequence);
    elseif Vfound==0 && Jfound==1
        fprintf(fid8,'>%d\r\n',i);
        fprintf(fid8,'%s\r\n',seq(i).Sequence);
    elseif Jfound==0
        fprintf(fid9,'>%d\r\n',i);
        fprintf(fid9,'%s\r\n',seq(i).Sequence);
    end
end
fclose('all');