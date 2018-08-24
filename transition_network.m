function [MFPT,minimaCoords,barrierCoords,minimaEnergy,barrierEnergy,barrierMinima,minimaHessian,barrierHessian,pathCoords,pathLength,pathEnergy,removeRows] = step_3b_transitions_from_landscape(gmm)
%This code takes a Gaussian mixture model probability density function (PDF) as
%input and finds the minima and minimum energy paths between minima on the
%energy landscape corresponding to the PDF, then calculates the mean first
%passage time between each state using Kramers turnover rates.

%Example:
%load('example_GMMs/example_GMM_villin_5D.mat')
%[MFPT,minimaCoords,barrierCoords,minimaEnergy,barrierEnergy,barrierMinima,minimaHessian,barrierHessian,pathCoords,pathLength,pathEnergy,removeRows] = transition_network(gmm);

%Input: gmm (GMM fit to data)
%Output: MFPT = matrix of mean first passage times between states (minima)
%        minimaCoords = locations of minima
%        barrierCoords = locations of saddle points between minima
%        minimaEnergy = energy of minima
%        barrierEnergy = energy of saddle points
%        barrierMinima = minima connected by each barrier (first row = minima connected by first barrier etc)
%        minimaHessian = hessian at minima
%        barrierHessian = hessian at barriers
%        pathCoords = coordinates along path between minima
%        pathLength = length along each path
%        pathEnergy = enery along each path
%        removeRows = list of barriers removed owing to lack of convergence or passing through an intermediate minimum

addpath(genpath('DNEB_clean'))
addpath(genpath('functions'))

N=gmm.NumComponents;
D=gmm.NumVariables;
%Convert covariance matrix to appropriate format if necessary
if size(gmm.Sigma,1)==1
    sigma_temp = zeros(gmm.NumVariables,gmm.NumVariables,N);
    for i=1:N
    sigma_temp(:,:,i) = diag(gmm.Sigma(:,:,i));
    end
    gmm = gmdistribution(gmm.mu,sigma_temp,gmm.ComponentProportion);
end

%Edit orig_dim and scaling_method variables (for scaling)
nonSphericalFlag = 1;
replace_string(D,nonSphericalFlag);
%input
generate_input_pp(gmm);
%run terminal commands to run dneb (make sure each file has appropriate permissions set)
fprintf('Running minimiser...\n')
cd DNEB_clean
[a,b] = system('./make_parallel','-echo');
%Find minima
[a,b] = system('./neb_minima','-echo');
%Get coordinates and energy of minima
[minimaCoords,minimaEnergy] = read_string_minima;
%Open parallel threads
setenv('OMP_NUM_THREADS','8');
%Find barriers
[a,b] = system('./neb_barriers','-echo');
%Get coordinates and energy of barriers
[barrierMinima,barrierEnergy,barrierCoords] = read_string_barriers;
%Get coordinates and energy along paths
[pathCoords,pathLength,pathEnergy] = read_string_paths(gmm);
cd ..




 removeRows = [];
% %Remove any path that goes through an intermediate minimum
for i = 1:size(pathCoords,2)
A = pdist2(minimaCoords(setdiff([1:size(minimaCoords,1)],barrierMinima(i,:)),:),pathCoords{i}');
[M,I] = min(A(:));
[I_row, I_col] = ind2sub(size(A),I);
if M<0.5
    removeRows = [removeRows,i];
    fprintf(['Path from minimum ',num2str(barrierMinima(i,1)),' to ',num2str(barrierMinima(i,2)),' seems to go through an intermediate minimum and has been removed.\n'])
end
end
% 
% %remove rows with NaN - these have not converged
removeRowsNan = find(isnan(barrierEnergy));
removeRows = union(removeRows,removeRowsNan);


% %Plot individual paths that are retained or removed
delete paths/*.png
delete paths/removed/*.png
delete paths/retained/*.png
for i = 1:size(pathCoords,2)
figure
scatter(minimaCoords(:,1),minimaCoords(:,2))
hold on
scatter(pathCoords{i}(1,:),pathCoords{i}(2,:),'rx')
if ismember(i,removeRows)==1
saveas(gcf,['paths/removed/Path ',num2str(i),'.png'])
else
  saveas(gcf,['paths/retained/Path ',num2str(i),'.png'])
end
end
%Plot all retained paths on one figure
retainedRows = setdiff([1:size(pathCoords,2)],removeRows);
figure
scatter(minimaCoords(:,1),minimaCoords(:,2))
hold on
for i=retainedRows
    scatter(pathCoords{i}(1,:),pathCoords{i}(2,:),'rx')
end
saveas(gcf,['paths/RetainedPaths.png'])
close all

%Here it may be necessary to add any paths in "retained" that go through intermediate minima using the command below:
%removeRows = union(removeRows,[])

%Remove spurious barriers that were identified above
barrierCoords(removeRows,:)=[];
barrierEnergy(removeRows) = [];
barrierMinima(removeRows,:) = [];

%Calculate mean first passage times
minimaHessian = calculateHessian(gmm,minimaCoords);
barrierHessian = calculateHessian(gmm,barrierCoords);
MFPT = MFPTmatrixWithHessian(minimaEnergy',minimaHessian,barrierEnergy',barrierHessian,[barrierMinima(:,1) [1:size(barrierMinima(:,1),1)]' barrierMinima(:,2)],1,1);


