function gmm = protein_fit_GMM(xyzTransform,D)
%Fits a GMM to processed data in the required number of dimensions (D)
%(clustering is first done in 3D)

%Cluster first by fitting a GMM in 3D (note number of components is found
%by finding the optimal number of kmeans clusters according to evalclusters then adding 4 to increase the fit)
clust = zeros(size(xyzTransform,1),10);
for i=1:10
clust(:,i) = kmeans(xyzTransform(:,1:2),i,'emptyaction','singleton',...
'replicate',5,'MaxIter',2000);
end
va = evalclusters(xyzTransform(:,1:2),clust,'CalinskiHarabasz')
numClusters = va.OptimalK+4;
gmm = fitgmdist(xyzTransform(:,1:3),numClusters,'CovarianceType','diagonal','Options',statset('MaxIter',2000),'Replicates',20);
clusterX = cluster(gmm,xyzTransform(:,1:3));
%Now fit indivdual gaussians to each cluster and combine
N = gmm.NumComponents;
for i = 1:N
gmmIndividual{i} = fitgmdist(xyzTransform(clusterX==i,1:D),1,'CovarianceType','diagonal','Options',statset('MaxIter',2000),'Replicates',20);
end
for i = 1:N
mu(i,:) = gmmIndividual{i}.mu;
Sigma(:,:,i) = gmmIndividual{i}.Sigma;
end
for i = 1:N
weight(i) = size(xyzTransform(clusterX==i,:),1)/size(xyzTransform,1);
end
gmm = gmdistribution(mu,Sigma,weight);
end

