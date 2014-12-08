function conservedIdentificationMultipleYYC(file_name)

file_name1 = strrep(file_name,'.fasta','/MultipleYYC.fasta');
seq = fastaread(file_name1);
load('germlines.mat');
s = [num2str(size(seq,1)),' sequence loading complete'];
disp(s);
new_file_name1 = strrep(file_name,'.fasta','/number.txt');
fid1 = fopen(new_file_name1,'a');
new_file_name2 = strrep(file_name,'.fasta','/CDR3.txt');
fid2 = fopen(new_file_name2,'a');
new_file_name3 = strrep(file_name,'.fasta','/Jgroup.txt');
fid3 = fopen(new_file_name3,'a');
new_file_name4 = strrep(file_name,'.fasta','/Vgene.txt');
fid4 = fopen(new_file_name4,'a');
new_file_name5 = strrep(file_name,'.fasta','/Vallele.txt');
fid5 = fopen(new_file_name5,'a');
count = 0;

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
    
    Jposition = Jposition(end);
    mleast = length(sequence);
    cleast = 0;
    vpleast = 0;
    for j = 1:size(Vposition,2)
        if Vposition(j)>=60
            sequence2 = sequence(1:Vposition(j));
            sbinary = zeros(size(sequence2));
            sbinary(sequence2=='A') = 1;
            sbinary(sequence2=='C') = 4;
            sbinary(sequence2=='G') = 9;
            sbinary(sequence2=='T') = 16;
            l = length(sbinary);
            m1 = zeros(size(Vbinary,1),l);
            m2 = zeros(size(Vbinary,1),l);
            for k = 1:size(Vbinary,1)
                m1(k,:) = Vbinary{k,1}(end-l+1:end);
                m2(k,:) = sbinary;
            end
            mutation = m1 - m2;
            mcount = (mutation~=0);
            msum = sum(mcount,2);
            [mmin,p] = minCI(msum);
            if mmin<=mleast
                mleast = mmin;
                cleast = l;
                vpleast = Vposition(j);
                Vallele = char(Vname{p(1)});
                star = strfind(char(Vname{p(1)}),'*');
                Vgene = char(Vname{p(1)}(1:star-1));
                if length(p)>1
                    for k = 2:length(p)
                        Vallele = [Vallele,'|',char(Vname{p(k)})];
                        star = strfind(char(Vname{p(k)}),'*');
                        if isempty(strfind(Vgene,char(Vname{p(k)}(1:star-1))))
                            Vgene = [Vgene,'|',char(Vname{p(k)}(1:star-1))];
                        end
                    end
                end
            end
        end
    end
    if mleast/cleast<0.15    %may need to adjust threshold here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fprintf(fid1,'%d %d %d %d\r\n',str2double(seq(i).Header),rev,vpleast,Jposition);
        if mod(Jposition - vpleast - 19,3)==0
            fprintf(fid2,'%d %s %s\r\n',str2double(seq(i).Header),sequence(vpleast+6:Jposition-14),nt2aa(sequence(vpleast+6:Jposition-14),'ACGTonly',false));
        else
            fprintf(fid2,'%d %s outofframe\r\n',str2double(seq(i).Header),sequence(vpleast+6:Jposition-14));
        end  
        fprintf(fid3,'%d %d\r\n',str2double(seq(i).Header),Jgene);  
        fprintf(fid4,'%d %s %d %d\r\n',str2double(seq(i).Header),Vgene,mleast,cleast);
        fprintf(fid5,'%d %s %d %d\r\n',str2double(seq(i).Header),Vallele,mleast,cleast);
        count = count + 1;
    end   
end
fclose('all');
s = [num2str(count),' sequence(s) added'];
disp(s);

