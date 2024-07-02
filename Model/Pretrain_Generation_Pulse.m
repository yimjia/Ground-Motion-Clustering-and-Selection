%% Pretrain Data Generation
clear
close all
clc

tic


%% Load Data
load GMIM_Pulse.mat


%% Number of GMs
N_GM = size(LogSa,1);


%% Number of Pretrained GMs (User-Defined)
N_Pre = 10000;


%% Repetitions
Rep = ceil((N_Pre/N_GM)/10)*10;


%% Sa
% Find scaled Sa
ID = Scale_Period/(T(2)-T(1));
Scale_Sa = LogSa(1,ID);

% Shift to the scaled Sa
LogSa = LogSa-Scale_Sa;


% Include noise
for i = 1:N_GM
    for j = 1:size(LogSa,2)
        T_Ceil(j+1,1) = (T(2)-T(1))*(j);
        if j == ID
            LogSa_Pre_Ceil{i,1}(:,j) = normrnd(LogSa(i,j),abs(LogSa(i,j))*0.00,[Rep,1]);
        else
            rng(1+i+j);
            LogSa_Pre_Ceil{i,1}(:,j) = normrnd(LogSa(i,j),abs(LogSa(i,j))*(abs(j-ID)*0.0015+0.20),[Rep,1]);
        end
    end
end

LogSa_Pre_Raw = cell2mat(LogSa_Pre_Ceil);


% Smooth
for i = 1:size(LogSa_Pre_Raw,1)
    LogSa_Pre(i,:) = smooth(T,LogSa_Pre_Raw (i,:),(T(2)-T(1))*5,'rloess');
end


% Inverse shift
LogSa = LogSa+Scale_Sa;

LogSa_Pre = LogSa_Pre+Scale_Sa;


% Find Scaled Sa
LogSa_Pre(:,ID) = Scale_Sa;


%% Save Data
LogSa = LogSa_Pre;

save(sprintf('GMIM_Pulse_Pretrain.mat'),'LogSa');

toc


