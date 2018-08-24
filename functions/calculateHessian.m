function Hessian = calculateHessian(gmm,points)
%Input GMM - Gaussian mixture model
%      points - points where to calculate Hessian (minima and/or barriers)
%Output H - Hessian at the points. Each element  of cell corresponds to a row in
%points and consists of the full Hessian matrix
N=gmm.NumComponents;
D = gmm.NumVariables;
P = gmm.ComponentProportion;
num_points = size(points,1);
%Initialize Hessian matrices at each point
H = cell(num_points,1);
[H{:,1}] = deal(zeros(D,D));
Hessian = zeros(num_points,D,D);
for i =1:size(points,1)
    x = points(i,:);
for j = 1:N
    gmmSeparate{j} = gmdistribution(gmm.mu(j,:),gmm.Sigma(:,:,j));
    H{i,1} = H{i,1} + P(j)*pdf(gmmSeparate{j},x) * (inv(gmm.Sigma(:,:,j))*((x-gmm.mu(j,:))'*(x-gmm.mu(j,:)))*inv(gmm.Sigma(:,:,j))-inv(gmm.Sigma(:,:,j)));
end
Hessian(i,:,:) = -H{i,1}/pdf(gmm,points(i,:));
end

