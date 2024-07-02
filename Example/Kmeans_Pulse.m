%% K-Means Clustering for Latent Features
clc
clear
close all

tic

%% Load Results
load GMIM_Pulse.mat 

Encode = readmatrix(sprintf('Data_Pulse_Finetune/Encode'));
Loss = readmatrix(sprintf('Data_Pulse_Finetune/Loss'));
Loss_Epoch = readmatrix(sprintf('Data_Pulse_Finetune/Min_Loss_Epoch'));    
DecodeSa_Raw = readmatrix(sprintf('Data_Pulse_Finetune/Decode_Sa'));

DecodeSa = exp(DecodeSa_Raw);


%% Cluster Analysis
Cluster = {'Num of Clusters = 1'
           'Num of Clusters = 2'
           'Num of Clusters = 3'
           'Num of Clusters = 4'
           'Num of Clusters = 5'
           'Num of Clusters = 6'
           'Num of Clusters = 7'
           'Num of Clusters = 8'
           'Num of Clusters = 9'
           'Num of Clusters = 10'};

NCluster = [1:10]';
Encode_N = Encode;


%% K-means
fprintf('K-means\n');

Label_KM = [];

for j = 1:length(Cluster)
    [Label_KM{1,j},C{1,j},sumd{1,j}] = kmeans(Encode_N,NCluster(j,1),'MaxIter',1000, ...
                                                                     'Replicates',1000); 
                                                                 
    CH_En_Raw = evalclusters(Encode_N,Label_KM{1,j},'CalinskiHarabasz');
    CH_En_KM(1,j) = CH_En_Raw.CriterionValues;

    DB_En_Raw = evalclusters(Encode_N,Label_KM{1,j},'DaviesBouldin');
    DB_En_KM(1,j) = DB_En_Raw.CriterionValues;
    
    Sil_En_Raw = evalclusters(Encode_N,Label_KM{1,j},'Silhouette','distance','Euclidean');
    Sil_En_KM(1,j) = Sil_En_Raw.CriterionValues;

    fprintf('%s finsihed.\n',Cluster{j,1});                                                    
end
      

for i = 1:size(Encode_N,1)
    for j = 1:length(NCluster)
        Dist(1,j) = sum(sumd{1,j});
    end
end

fprintf('\n');                                                    


%% Save Results
clearvars -except Loss Loss_Epoch ...
                  DecodeSa ...
                  Cluster NCluster ...
                  Encode ...
                  Label_KM    Label_GM     Label_SP ...
                  Dist ...
                  CH_En_KM    DB_En_KM      Sil_En_KM ...
                  CH_En_GM    DB_En_GM      Sil_En_GM ...
                  Sim_KM      Sim_GM
              
              
save('Cluster_Pulse.mat')

toc


