function [vlabel,geneGroup] = Vties(file_name)
%Make vties based on average V length and mutation

dirString = strrep(file_name,'.fasta','');
inFileName = [dirString,'/Vgene.txt'];
fid = fopen(inFileName,'r');
formatSpec = '%d%s%d%d';
dataArray = textscan(fid,formatSpec,'Delimiter',' ','ReturnOnError',false);
order = dataArray{1,1};
Vgene = dataArray{1,2};
mutation = dataArray{1,3};
vlength = dataArray{1,4};
clear dataArray;

average_mutation = mean(mutation);
average_vlength = mean(vlength);
mutation_frequency = average_mutation/average_vlength;

if mutation_frequency<0.15
    mutation_frequency = 0.05;
elseif mutation_frequency<0.3 && mutation_frequency>=0.15
    mutation_frequency = 0.15;
elseif mutation_frequency>0.3
    mutation_frequency = 0.3;
end

if average_vlength<150
    average_vlength = 100;
elseif average_vlength<200 && average_vlength>=150
    average_vlength = 150;
elseif average_vlength>=200
    average_vlength = 200;
end

dbName = ['mix_sm_',num2str(average_vlength),'_',num2str(mutation_frequency),'.mat'];
load(dbName);
load('HHVdatabase.mat');

geneGroup = cell(49,2);
for i = 1:size(genename,1)
    geneGroup{i,1} = char(genename{i,1});
    string = char(genename{i,1});
    for j = 1:size(M,1)
        if ~isempty(strfind(char(M{j,1}),char(genename{i,1})))
            string = [string,'|',char(M{j,1})];
        end
    end
    c = strsplit(string,'|');
    c = unique(c);
    string = '';
    for j = 1:length(c);
        string = [string,'|',char(c{j})];
    end
    geneGroup{i,2} = string(2:end);
end
vlabel = geneGroup(:,2);
vlabel = unique(vlabel,'stable');

Vties = cell(size(Vgene));
for i = 1:size(Vgene,1)
    c = strsplit(char(Vgene{i,1}),'|');
    temp_string = '';
    for j = 1:size(c,2)
        p = find(strcmp(genename,char(c{j})),1);
        temp_string = [temp_string,'|',char(geneGroup{p,2})];
    end
    c = strsplit(temp_string(2:end),'|');
    c = unique(c);
    temp_string = '';
    for j = 1:size(c,2)            
        temp_string = [temp_string,'|',char(c{j})];
    end
    temp_string = temp_string(2:end);
    Vties{i,1} = temp_string;
end

outFileName = [dirString,'/Vties.txt'];
fid = fopen(outFileName,'w');
for i = 1:size(Vties,1)
    fprintf(fid,'%d %s %d %d\r\n',order(i,1),char(Vties{i,1}),mutation(i,1),vlength(i,1));
end

fclose('all');