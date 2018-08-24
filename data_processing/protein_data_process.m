function [xyzTransform,rmsd,newpoints] = protein_data_process(xyz,xyzTailPoints,varargin)
%Input: xyz = matrix: rows are timepoints; columns are coordinates of
%       carbon alpha atoms (x1,y1,z1, ...). Obtain raw data from DE Shaw
%       Research.

%       xyz_tailPoints = vector of residues to remove (e.g. atoms on the tail).
%       e.g. for villin use [1:5,31:35]

%OPTIONAL inputs:
%       xyzNative = matrix: coordinates of native state; obtainable from PDB file.
%       Rows are atoms; columns are (x,y,z). (Use mean of xyz otherwise).
%       For villin, see read_pdb.m to obtain this.

%       cutOffDistance = distance to cut off pairwise distances. (Use 8A
%       otherwise).

%       subSampleFactor = factor to subsample by (whole number > 1).
%       Otherwise don't subsample.

%Output: xyzTransform = new observations after PCA

%Example (for villin):
%protein_data_process(xyz,[1:5,31:35],xyzNative,8,5)

if nargin > 2
    xyzNative = varargin{1};
else
    x_mean=xyz_mean(:,1:3:end);
    y_mean=xyz_mean(:,2:3:end);
    z_mean=xyz_mean(:,3:3:end);
    xyzNative = [x_mean' y_mean' z_mean'];
end
if nargin>3
        cutOffDistance = varargin{2};
else
        cutOffDistance = 8;
end
if nargin>4
        subSampleFactor = varargin{3};
else
        subSampleFactor = 1;
end
%--------------------------------------------------------------------------
%Remove specified atoms
xyzTailColumns = [];
for i = xyzTailPoints
xyzTailColumns = [xyzTailColumns, 3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3];
end
xyz(:,xyzTailColumns) = [];
xyzNative(xyzTailPoints,:) = [];

%Align with native state and calculate RMSD
for i = 1:size(xyz,1)
points{i} = [xyz(i,1:3:end)' xyz(i,2:3:end)' xyz(i,3:3:end)'];
[~, newpoints{i}] = procrustes(xyzNative,points{i},'scaling',false);
rmsd(i,1) = sqrt(mean(sum((xyzNative - newpoints{i}).^2,2)));
end

%Subsample
newpoints = newpoints(1:subSampleFactor:end);
rmsd = rmsd(1:subSampleFactor:end);

%Perform PCA on distances < cutOffDistance
a = cellfun(@(x) pdist(x),newpoints,'UniformOutput',false);
b = cell2mat(a(:));
b(b>cutOffDistance)=cutOffDistance;
[~,xyzTransform,~,~,explained] = pca(b);

end

