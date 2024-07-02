%% Conditional Spectra
clear
close all
clc

tic


%% Load Clusters
load GMIM_Pulse.mat
LogSa_Pulse = LogSa;

load GMIM_NonPulse.mat
LogSa_NonPulse = LogSa;

load Cluster_Pulse.mat
N_Clu_Pulse = find(Sil_En_KM == max(Sil_En_KM));
Label_Pulse = Label_KM{N_Clu_Pulse};

load Cluster_NonPulse.mat
N_Clu_NonPulse = find(Sil_En_KM == max(Sil_En_KM));
Label_NonPulse = Label_KM{N_Clu_NonPulse};


%% Number of GMs
N_GM = 20;

% From USGS Disaggregation Tool
R = 16.54;
epsilon = 0.90;

% From NIST (2011)
Por = 0.905-0.188*R+1.337*epsilon;
Prob_Pulse = exp(Por)/(1+exp(Por));

N_Pulse = round(N_GM*Prob_Pulse);
N_NonPulse = N_GM-N_Pulse;


%% Pusle-Type GMs in Clusters
for i = 1:N_Clu_Pulse
    N_Clu_PulseGM(i,1) = length(find(Label_Pulse == i));
    LogSa_Pulse_Clu{i,1} = LogSa_Pulse(find(Label_Pulse == i),:);
    LogSa_Pulse_Clu_Mean{i,1} = mean(LogSa_Pulse(find(Label_Pulse == i),:));
end

for i = 1:N_Clu_Pulse
    N_Pulse_Clu(i,1) = round(N_Pulse*N_Clu_PulseGM(i,1)/length(Label_Pulse));
    if i == N_Clu_Pulse
        N_Pulse_Clu(i,1) = N_Pulse-sum(N_Pulse_Clu(1:i-1,1));
    end
end


%% Non-Pusle-Type GMs in Clusters
for i = 1:N_Clu_NonPulse
    N_Clu_NonPulseGM(i,1) = length(find(Label_NonPulse == i));
    LogSa_NonPulse_Clu{i,1} = LogSa_NonPulse(find(Label_NonPulse == i),:);
    LogSa_NonPulse_Clu_Mean{i,1} = mean(LogSa_NonPulse(find(Label_NonPulse == i),:));
end

for i = 1:N_Clu_NonPulse
    N_NonPulse_Clu(i,1) = round(N_NonPulse*N_Clu_NonPulseGM(i,1)/length(Label_NonPulse));
    if i == N_Clu_NonPulse
        N_NonPulse_Clu(i,1) = N_NonPulse-sum(N_NonPulse_Clu(1:i-1,1));
    end
end


%% Load Location-Related Data
load TargetSa.mat
load GMPE.mat


%% Period
T_Raw = TargetSa(:,1);


%% Interpolation
Target = pchip(T_Raw,TargetSa(:,2),T);
TargetSa = [T Target];

GMPE_median = pchip(T_Raw,GMPE(:,2),T);
GMPE_mediansigma = pchip(T_Raw,GMPE(:,3),T);
GMPE = [T GMPE_median GMPE_mediansigma];


%% Target Event
T_star = 1;

T_ID = find(TargetSa(:,1) == T_star);

sigma_star = log(GMPE(T_ID,3))-log(GMPE(T_ID,2));

epsilon = (log(TargetSa(T_ID,2))-log(GMPE(T_ID,2)))/sigma_star;


%% Ground motion prediction model calcs
periods = TargetSa(:,1);
for i=1:length(periods)
    median_sa(i,1) = GMPE(i,2);
    sigma(i,1) = log(GMPE(i,3))-log(GMPE(i,2));
    rho(i,1) = baker_jayaram_correlation(T_star, periods(i,1));
end

CMS = exp( log(median_sa) + epsilon.* rho .* sigma);
CMS_sigma = sigma .* sqrt(1-rho.^2);

for i = 1:length(sigma)
    for j = 1:length(sigma)
        CS_COV(i,j) = baker_jayaram_correlation(periods(i,1), periods(j,1))*sigma(i,1)*sigma(j,1);
    end
end


%% Realization
Mean = log(CMS);
COV = CS_COV-CS_COV(:,T_ID)*CS_COV(T_ID,:)/sigma_star^2;


%% Realizations
N_Real = 100;
fprintf('Run %.0f Sets of Realizations. \n',N_Real);

for m = 1:N_Real

rng(m);
Real = lhsnorm(Mean,COV,20);

T_Raw2 = [logspace(-2,log10(1),15) logspace(log10(1.1),log10(2),10)]';
T_Raw2 = unique(round(T_Raw2*100)/100);

for i = 1:length(T_Raw2)
    ID(i,1) = find(T > 0.999*T_Raw2(i,1) & T < 1.001*T_Raw2(i,1));
end

Real2 = Real(:,ID);

for j = 1:length(Real2(1,:))
    Real2_Std(1,j) = std(Real2(:,j));
end

Real3 = pchip(T_Raw2,Real2,T);

Real_Plot{m,1} = Real3;


%% Realization to Cluster Centriod
LogSa_Clu_Mean = [LogSa_Pulse_Clu_Mean
                  LogSa_NonPulse_Clu_Mean];

for i = 1:length(Real3(:,1))
    for j = 1:(N_Clu_Pulse+N_Clu_NonPulse)
        Error_Mean_Raw(i,j) = sum((Real3(i,:)-LogSa_Clu_Mean{j,1}).^2);
    end
end


%% Different Order Realization to Cluster Centriod 
N_Clu = [N_Pulse_Clu
         N_NonPulse_Clu];
  
Comb = unique(perms([1:length(N_Clu)]), 'rows');

     
for i = 1:length(Comb)
    Error_Mean{i,1} = Error_Mean_Raw;
    for j = 1:length(Comb(i,:))        
        [B1, RealID{i,j}] = mink(Error_Mean{i,1}(:,Comb(i,j)),N_Clu(Comb(i,j),1));
        Error_Mean{i,1}(RealID{i,j},:) = 9999;
    end
end


for i = 1:length(Comb)
    RealID_P{i,1} = {};
    for j = 1:length(N_Pulse_Clu)
        RealID_P{i,1}(j,1) = {RealID{i,find(Comb(i,:) == j)}};
    end

    RealID_NP{i,1} = {};  
    for j = 1:length(N_NonPulse_Clu)
        RealID_NP{i,1}(j,1) = {RealID{i,find(Comb(i,:) == j+length(N_Pulse_Clu))}};
    end
end

          
%% Pulse-Type GM Selection
for p = 1:length(Comb)
    for q = 1:length(LogSa_Pulse_Clu)
        for i = 1:length(LogSa_Pulse_Clu{q,1}(:,1))
            for j = 1:length(RealID_P{p,1}{q,1})
                Error_P{p,1}{q,1}(i,j) = sum((Real3(RealID_P{p,1}{q,1}(j,1),:)-LogSa_Pulse_Clu{q,1}(i,:)).^2);
            end
        end
    end
end


for p = 1:length(Comb)
    for q = 1:length(LogSa_Pulse_Clu)
        for j = 1:length(RealID_P{p,1}{q,1})
            [~, GMID_P{p,1}{q,1}(j,1)] = mink(Error_P{p,1}{q,1}(:,j),1);
            Error_P{p,1}{q,1}(GMID_P{p,1}{q,1}(j,1),:) = 9999;
        end
    end
end


%% Non-Pulse-Type GM Selection
for p = 1:length(Comb)
    for q = 1:length(LogSa_NonPulse_Clu)
        for i = 1:length(LogSa_NonPulse_Clu{q,1}(:,1))
            for j = 1:length(RealID_NP{p,1}{q,1})
                Error_NP{p,1}{q,1}(i,j) = sum((Real3(RealID_NP{p,1}{q,1}(j,1),:)-LogSa_NonPulse_Clu{q,1}(i,:)).^2);
            end
        end
    end
end


for p = 1:length(Comb)
    for q = 1:length(LogSa_NonPulse_Clu)
        for j = 1:length(RealID_NP{p,1}{q,1})
            [~, GMID_NP{p,1}{q,1}(j,1)] = mink(Error_NP{p,1}{q,1}(:,j),1);
            Error_NP{p,1}{q,1}(GMID_NP{p,1}{q,1}(j,1),:) = 9999;
        end
    end
end


%% Selected GMs
LogCS = log(CMS)';
LogCS_Sigma = CMS_sigma';

for p = 1:length(Comb)
    LogSa_Sel{p,1} = [];

    for o = 1:length(N_Pulse_Clu)
        LogSa_Sel{p,1} = [LogSa_Sel{p,1}
                           LogSa_Pulse_Clu{o,1}(GMID_P{p,1}{o,1},:)];
    end
    for o = 1:length(N_NonPulse_Clu)
        LogSa_Sel{p,1} = [LogSa_Sel{p,1}
                           LogSa_NonPulse_Clu{o,1}(GMID_NP{p,1}{o,1},:)];
    end

    % LogSa_Sel{p,1} = [LogSa_Pulse_Clu{1,1}(GMID_P{p,1}{1,1},:)
    %                   LogSa_Pulse_Clu{2,1}(GMID_P{p,1}{2,1},:)
    %                   LogSa_Pulse_Clu{3,1}(GMID_P{p,1}{3,1},:)
    %                   LogSa_NonPulse_Clu{1,1}(GMID_NP{p,1}{1,1},:)
    %                   LogSa_NonPulse_Clu{2,1}(GMID_NP{p,1}{2,1},:)
    %                   LogSa_NonPulse_Clu{2,1}(GMID_NP{p,1}{3,1},:)];

    Mean_LogSa_Sel{p,1} = mean(LogSa_Sel{p,1});
    
    for j = 1:length(LogSa_Sel{p,1}(1,:))
        Std_LogSa_Sel{p,1}(1,j) = std(LogSa_Sel{p,1}(:,j));
    end
    
    Error_CS(p,1) = sum((Mean_LogSa_Sel{p,1}-LogCS).^2)+sum((Std_LogSa_Sel{p,1}-LogCS_Sigma).^2);
end

    
% Best Combination
ID_Comb = find(Error_CS == min(Error_CS),1);
ID_Pulse{m,1} = GMID_P{ID_Comb,1};
ID_NonPulse{m,1} = GMID_NP{ID_Comb,1};
Error_CS_Sim(m,1) = min(Error_CS);


% Print process
if mod(m,N_Real*0.1)==0
    fprintf('Finish %.0f Sets of Realizations. \n',m);
end   
end


%% Selected GM ID
ID_Min = find(Error_CS_Sim == min(Error_CS_Sim));

for j = 1:N_Clu_Pulse
    ID_Pulse_Raw{j,1} = find(Label_Pulse == j);
    ID_Pulse_Final_Cell{j,1} = ID_Pulse_Raw{j,1}(ID_Pulse{ID_Min,1}{j,1});
end
ID_Pulse_Final = cell2mat(ID_Pulse_Final_Cell);

for j = 1:N_Clu_NonPulse
    ID_NonPulse_Raw{j,1} = find(Label_NonPulse == j);
    ID_NonPulse_Final_Cell{j,1} = ID_NonPulse_Raw{j,1}(ID_NonPulse{ID_Min,1}{j,1});
end
ID_NonPulse_Final = cell2mat(ID_NonPulse_Final_Cell);

LogSa_Final = [LogSa_Pulse(ID_Pulse_Final,:)
               LogSa_NonPulse(ID_NonPulse_Final,:)];

Mean_LogSa_Final = mean(LogSa_Final);
for j = 1:length(LogSa_Final(1,:))
    Std_LogSa_Final(1,j) = std(LogSa_Final(:,j));
end


%% Save Data
save('GM_Selection.mat','ID_Pulse_Final','ID_NonPulse_Final',...
                        'LogSa_Final','Mean_LogSa_Final','Std_LogSa_Final', ...
                        'CMS', 'CMS_sigma','T');

toc

