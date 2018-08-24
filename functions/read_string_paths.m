function [pathCoords,pathLength,pathEnergy] = read_string_path(gmm)
%Read barrier file and write to variable
fid = fopen('../DNEB_clean/output/path.txt','r');
f=textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');


idx = find(~cellfun(@isempty,f{1}));
numDimensions = find(diff(idx)~=1,1);
numBarriers = size(idx,1)/numDimensions;

    for i = 1:numBarriers
        for j = 1:numDimensions
            pathCoords{i}(j,:) = str2num(f{1}{idx(j+numDimensions*(i-1))});
        end
    end
numImages = size(pathCoords{1},2);
pathLength = zeros(numBarriers,numImages);
for i =1:numBarriers
for j =1:numImages
    pathEnergy(i,j) = -log(pdf(gmm,pathCoords{i}(:,j)'));
        if j>1
        pathLength(i,j) = pathLength(i,j-1) + norm(pathCoords{i}(:,j) - pathCoords{i}(:,j-1));
        end
end
end

fclose(fid);
end

