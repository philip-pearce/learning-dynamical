function replace_string(D,nonSphericalFlag)
%Read file and write to variable
fid = fopen('DNEB_clean/source/fitness.f90','r');
f=textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
txt = f{:};
fclose(fid);
%Edit ORIG_DIM variable
split_string = strsplit(txt{8},'=');
new_string = strcat({split_string{1,1}},{'= '},{num2str(D)});
txt{8} = new_string{1};
%Edit SCALING_METHOD variable
split_string = strsplit(txt{9},'=');
new_string = strcat({split_string{1,1}},{'= '},{num2str(nonSphericalFlag+1)});
txt{9} = new_string{1};
%Write to file
fid = fopen('DNEB_clean/source/fitness.f90','w');
fprintf(fid,'%s\n',txt{1:end})
fclose(fid);

end