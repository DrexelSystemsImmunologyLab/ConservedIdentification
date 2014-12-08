function n = strcmpn(string1,string2)

n = 0;
for i = 1:length(string1)
    if strcmp(string1(i),string2(i))~=1 && strcmp(string1(i),'N')~=1 && strcmp(string2(i),'N')~=1
        n = n + 1;
    end
end