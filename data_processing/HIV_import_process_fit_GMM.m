%This script converts raw data from Zanini et al (p17 sequences) to a GMM.

%Load consensus sequence(s)
load('p17/consensus_p17.mat')
%Number of consensus sequences
C = 21;
s{1} = consensus;
s{2} = consensus_A1;
s{3} = consensus_A2;
s{4} = consensus_B;
s{5} = consensus_C;
s{6} = consensus_D;
s{7} = consensus_F1;
s{8} = consensus_G;
s{9} = consensus_H;
s{10} = consensus_K;
s{11} = consensus_01_AE;
s{12} = consensus_02_AG;
s{13} = consensus_03_AB;
s{14} = consensus_04_CPX;
s{15} = consensus_06_CPX;
s{16} = consensus_07_BC;
s{17} = consensus_08_BC;
s{18} = consensus_10_CD;
s{19} = consensus_11_CPX;
s{20} = consensus_12_BF;
s{21} = consensus_14_BG;
for i=1:C
reads{i,1} = '0';
days{i,1} = '0';
end
%Extract relevant data from files
data=[];
%set ID to zero for consensus sequences
id = zeros(C,1);
fprintf('Extracting data...');
for i=1:11
    if i~=4
        filename = ['p17/haplotypes_p' num2str(i) '_p17.fasta']
        read_temp = fastaread(filename);
        data = [data;read_temp];
        id_temp = i*ones(size(read_temp,1),1);
        %Patient ID for each sample
        id = [id;id_temp];
    end
end
        clear filename read_temp id_temp

%Extract information from data structure
counter = C;
            for j = 1:size(data,1)
                counter = counter+1;
            temp = strsplit(data(j).Header,'_');
            temp2 = strsplit(data(j).Header,'n.reads: ');
            days(counter,1) = temp(2);
            reads(counter,1) = temp2(2);
            s{counter} = data(j).Sequence;
            clear temp temp2
            end
            
%Align sequences
fprintf('Aligning sequences...');
sequences = multialign(s);

%Convert relevant strings to numbers
days = cellfun(@str2num, days);
reads = cellfun(@str2num, reads);

%Sort by sample time
fprintf('Sorting sequences by sample time (from infection)...');
[days , I] = sort(days);
reads = reads(I);
id = id(I);
sequences = sequences(I,:);

%Array - 0 = consensus (original); 1 = non-consensus
for i =1:size(sequences,1)
diff(i,:) = sequences(1,:) ~= sequences(i,:);
end
diff = double(diff);
%Marker for number of days
[~,~,id_days] = unique(days);

%First remove patient 7
diff(id==7,:) = [];
days(id==7) = [];
reads(id==7) = [];
sequences(id==7,:) = [];
id(id==7) = [];
%normalise to reads per thousand
for i = 1:size(reads,1)
reads_scaled(i) = floor(reads(i)/sum(reads(id==id(i)&days==days(i)))*1000);
end
%note the first few (consensus) rows should be populated with zero reads
reads_scaled(isnan(reads_scaled)==1)=0;
diff_reads = [];
sequences_reads =[];
id_reads = [];
%timepoints (days) for each patient (ID)
for i = 1:max(id)
    current_id_days{i} = unique(days(id==i));
end
counter = 0;
idx = [];
%loop over sequences
for i = 1:size(diff,1)
    counter = counter+1;
    %ignore conesnsus sequences
    if id(i)>0
        %Only take final 5 timepoints
    if days(i)>current_id_days{id(i)}(end-5)
        %locations of sequences we keep
        idx = [idx;counter];
        %sequences repeated by the number of times read
        diff_reads = [diff_reads; repmat(diff(i,:),reads_scaled(i),1)];
        sequences_reads = [sequences_reads; repmat(sequences(i,:),reads_scaled(i),1)];
        %associated patient for each repeat
        id_reads = [id_reads;repmat(id(i),reads_scaled(i),1)];
    end
    end
end
%PCA
[coeff,score,latent,tsquared,explained,mu]=pca(diff_reads);
%unique sequences in PC space
score_unique = (diff(idx,:) - mean(diff_reads))*coeff;
%consensus sequence(s) in PC space
score_consensus = (diff(1,:) - mean(diff_reads))*coeff;

for i = 1:11
patientMeans(i,:) = mean(score(id_reads==i,:));
end

%Fit GMM
%First fit gaussian to each patient
clear mu Sigma
clusterX = id_reads;
clusterX(clusterX>4) = clusterX(clusterX>4)-1;
clusterX(clusterX>6) = clusterX(clusterX>6)-1;
id_reads = clusterX;
D=10;
N = max(clusterX);
for i = 1:N
gmmIndividual{i} = fitgmdist(score(id_reads==i,1:D),1,'CovarianceType','diagonal','Options',statset('MaxIter',2000),'Replicates',20);
end
for i = 1:N
mu(i,:) = gmmIndividual{i}.mu;
Sigma(:,:,i) = gmmIndividual{i}.Sigma;
end
for i = 1:N
weight(i) = size(score(id_reads==i,:),1)/size(score,1);
end
gmm = gmdistribution(mu,Sigma,weight);