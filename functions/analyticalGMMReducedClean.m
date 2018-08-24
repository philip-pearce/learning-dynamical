function [gmmReduced,V] = analyticalGMMReducedClean(gmm,dim)
%INPUTS: gmm = Gaussian mixture model to be reduced in dimension
%        dim = reduced number of dimensions

%OUTPUTS: gmmReduced = Gaussian mixture model in reduced dimensions
%         V = PCA transformation

N=gmm.NumComponents;
%Convert covariance matrix to appropriate format if necessary
if size(gmm.Sigma,1)==1
    sigma_temp = zeros(gmm.NumVariables,gmm.NumVariables,N);
    for i=1:N
    sigma_temp(:,:,i) = diag(gmm.Sigma(:,:,i));
    end
    gmm = gmdistribution(gmm.mu,sigma_temp,gmm.ComponentProportion);
end



%Perform PCA on means ONLY
C = cov(gmm.mu-mean(gmm.mu));
[V,S] = svd(C);
%Apply pca to means
new_mean = (gmm.mu-mean(gmm.mu))*V(:,1:dim);

        % Apply PCA to covariance matrix
        new_sigma = zeros(dim, dim, N);
        for i=1:N
            new_sigma(:,:,i) = V(:,1:dim)' * gmm.Sigma(:,:,i) * V(:,1:dim);
        end

        %% Generate reduced GMM
         gmmReduced = gmdistribution(new_mean, new_sigma, gmm.ComponentProportion);

end

