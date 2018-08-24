function Gradient = calculateGradient(gmm,points)
%Input GMM - Gaussian mixture model
%      points - points where to calculate Hessian (minima and/or barriers)
%Output Gradient - gradient at the points. Each element  of cell corresponds to a row in
%the points variable


N=gmm.NumComponents;
D = gmm.NumVariables;
P = gmm.ComponentProportion;
num_points = size(points,1);
%Initialize gradient at each point
G = cell(num_points,1);
[G{:,1}] = deal(zeros(D,1));

for i =1:size(points,1)
    x = points(i,:);
for j = 1:N
    gmmSeparate{j} = gmdistribution(gmm.mu(j,:),gmm.Sigma(:,:,j));
    G{i,1} = G{i,1} - P(j)*pdf(gmmSeparate{j},x) * (inv(gmm.Sigma(:,:,j))*(x-gmm.mu(j,:))');
end
Gradient{i} = -G{i,1}/pdf(gmm,points(i,:));
end