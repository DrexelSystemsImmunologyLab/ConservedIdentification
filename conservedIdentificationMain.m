function conservedIdentificationMain(file_name,clone_threshold,cp_threshold,gap_method,subject,subset,tissue,disease,date,lab,person)

tic;

string = strrep(file_name,'.fasta','');
mkdir(string);
disp('making new directory done');

conservedIdentification(file_name);

disp('intact identification done');
FileInfo = dir([string,'/Jonly.fasta']);
if FileInfo.bytes~=0
    conservedIdentificationNoYYC(file_name);
end
disp('NoYYC identification done');

FileInfo = dir([string,'/MultipleYYC.fasta']);
if FileInfo.bytes~=0
    conservedIdentificationMultipleYYC(file_name);
end
disp('MultipleYYC identification done');

[vlabel,geneGroup] = Vties(file_name);

[germlines,vseqs] = makeTempSeq(file_name);

copy_number = copyNumber(file_name);
saved_name = [string,'/copy_number.mat'];
save(saved_name,'copy_number');
disp('Collapsing complete');

makeTable(file_name,gap_method,subject,subset,tissue,disease,date,lab,person,germlines,vseqs);
disp('Making master table complete');

% makeClones(file_name,clone_threshold,cp_threshold);
% disp('Making clones complete');
% 
% outputResult(file_name,vlabel,geneGroup)
% 
% collapseClones(file_name)
% disp('Collapsing clones complete');

% modifyClones(file_name)
% disp('Making clones with single copy complete');

% countMutationClone(file_name)
% mutationAnalysis(file_name)
% disp('Analyzing mutation complete');

s = [file_name,' done'];
disp(s);
toc;

clear;
fclose('all');