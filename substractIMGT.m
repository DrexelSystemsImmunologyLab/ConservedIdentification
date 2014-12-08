function name = substractIMGT(header)

p = strfind(header,'|');
name = header(p(1)+1:p(2)-1);

