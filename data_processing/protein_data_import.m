function [xyz] = protein_data_import(num_files)
%This function imports data from DE Shaw research (.dcd files) to Matlab.
%Requires: matdcd-1.0 (https://www.ks.uiuc.edu/Development/MDTools/matdcd/)

xyz = [];
j=1:num_files;
jj = cellstr(num2str(j(:),'%03d'))';
for i=1:num_files
%Change this line depending on file name
str_temp = strcat('2F4K-0-c-alpha-',jj(i),'.dcd');
%Change second input to readdcd depending on number of molecules
xyz = [xyz;readdcd(str_temp{1},1:35)];
end
end

