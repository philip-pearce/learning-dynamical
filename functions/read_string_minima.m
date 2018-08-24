function [minimaCoords,minimaEnergy] = read_string_minima
%Read minima file and write to variable
fid = fopen('output/minima.txt','r');
f=textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');

%number of minima
numMinima = round(size(f{1},1)/4);
%check minima file of expected size
if mod(size(f{1},1),4) ==0
    for i = 1:numMinima
        %Read in minima coords and energy using characteristics of output
        %file
        rowID = 4*(i-1)+2;
        trimmedString = strtrim(f{1}{rowID});
        coordsString = strsplit(trimmedString);
        dim = size(coordsString,2) - 1;
        for j = 1:dim
            minimaCoords(i,j) = str2num(coordsString{j+1});
        end
        trimmedString = strtrim(f{1}{rowID+1});
        energyString = strsplit(trimmedString);
        minimaEnergy(i) = str2num(energyString{2});
    end
else
    fprintf('ERROR: minima text file not of expected size. Not reading in olution - check file directly.')
    minimaCoords = 0;
    minimaEnergy = 0;
end
fclose(fid);
end

