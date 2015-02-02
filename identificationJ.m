function [jgroup, jposition] = identificationJ(sequence)
% finds J string in sequence, outputs J group and position of J in sequence

J1 = 'TGGTCACCGTCTCCTCAG'; %IGHJ1/4/5
J2 = 'TGGTCACTGTCTCCTCAG'; %IGHJ2
J3 = 'TGGTCACCGTCTCTTCAG'; %IGHJ3
J6 = 'CGGTCACCGTCTCCTCAG'; %IGHJ6

jgroup = 0;
jposition = 0;
sequence = upper(sequence);

p1 = strfind(sequence, J1);
p2 = strfind(sequence, J2);
p3 = strfind(sequence, J3);
p4 = strfind(sequence, J6);
if ~isempty(p1)
    jgroup = 1;
    jposition = p1(end);
elseif ~isempty(p2)
    jgroup = 2;
    jposition = p2(end);
elseif ~isempty(p3)
    jgroup = 3;
    jposition = p3(end);
elseif ~isempty(p4)
    jgroup = 6;
    jposition = p4(end);
end

%% shorten by 3 nts
if jgroup==0 && jposition==0    
    p1 = strfind(sequence, J1(1:end-3));
    p2 = strfind(sequence, J2(1:end-3));
    p3 = strfind(sequence, J3(1:end-3));
    p4 = strfind(sequence, J6(1:end-3));    
    if ~isempty(p1)
        jgroup = 1;
        jposition = p1(end);
    elseif ~isempty(p2)
        jgroup = 2;
        jposition = p2(end);
    elseif ~isempty(p3)
        jgroup = 3;
        jposition = p3(end);
    elseif ~isempty(p4)
        jgroup = 6;
        jposition = p4(end);
    end    
end

%% shorten by 6 nts
if jgroup==0 && jposition==0    
    p1 = strfind(sequence, J1(1:end-6));
    p2 = strfind(sequence, J2(1:end-6));
    p3 = strfind(sequence, J3(1:end-6));
    p4 = strfind(sequence, J6(1:end-6));    
    if ~isempty(p1)
        jgroup = 1;
        jposition = p1(end);
    elseif ~isempty(p2)
        jgroup = 2;
        jposition = p2(end);
    elseif ~isempty(p3)
        jgroup = 3;
        jposition = p3(end);
    elseif ~isempty(p4)
        jgroup = 6;
        jposition = p4(end);
    end    
end



    