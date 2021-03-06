function [barrierMinima,barrierEnergy,barrierCoords] = read_string_barriers
%Read barrier file and write to variable
fid = fopen('../DNEB_clean/output/barriers.txt','r');
f=textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
%Number of barriers
numBarriers = round((size(f{1},1)-1)/7);
%check barrier file of expected size
if mod((size(f{1},1)-1),7) ==0
    for i = 1:numBarriers
        rowID = 7*(i-1)+7;
        trimmedString = strtrim(f{1}{rowID});
        energyString = strsplit(trimmedString);
        trimmedString2 = strtrim(f{1}{rowID-5});
        minimaString = strsplit(trimmedString2);
        trimmedString3 = strtrim(f{1}{rowID-2});
        coordsString = strsplit(trimmedString3);
        barrierMinima(i,1) = str2num(minimaString{4});
        barrierMinima(i,2) = str2num(minimaString{6});
        if strcmp(energyString{3},'Infinity')==1
            barrierEnergy(i)=NaN;
        else
        barrierEnergy(i) = str2num(energyString{3});
        end
        dim = size(coordsString,2) - 3;
        for j = 1:dim
            barrierCoords(i,j) = str2num(coordsString{j+3});
        end
    end
else
    fprintf('ERROR: barrier text file not of expected size. Not reading in olution - check file directly.')
    barrierMinima=0;
    barrierEnergy=0;
    barrierCoords = 0;
end
fclose(fid);
end

