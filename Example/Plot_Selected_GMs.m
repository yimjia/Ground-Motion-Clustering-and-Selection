%% Plot the Selected GMs with Conditional Spectra
clc
close all
clear

tic


%% Load DAta
load GM_Selection.mat


%% Number of GMs
N_GM = length(LogSa_Final(:,1));
N_Pulse = length(ID_Pulse_Final);
N_NonPulse = length(ID_NonPulse_Final);


%% Conditional Spectra
LogCS = log(CMS)';
LogCS_Sigma = CMS_sigma';


%% Selected GMs
Fig = figure('units','normalized','Position',[0.05 0.05 0.88 0.88]);
hold on
box on
grid on

pbaspect([1 1 1])

p8 = plot(T,LogSa_Final(6:end,:), 'c-', 'linewidth', 2);
p7 = plot(T,LogSa_Final(1:5,:), '-.','color',[0 0.8 0],'linewidth', 2);

p1 = plot(T, log(CMS), 'k-', 'linewidth', 6);
p2 = plot(T, log(exp(log(CMS)+CMS_sigma)), 'k--', 'linewidth', 3);
p3 = plot(T, log(exp(log(CMS)-CMS_sigma)), 'k--', 'linewidth', 3);

p4 = plot(T,Mean_LogSa_Final, 'b-.', 'linewidth', 6);
p5 = plot(T, Mean_LogSa_Final+Std_LogSa_Final, 'b:', 'linewidth', 3);
p6 = plot(T, Mean_LogSa_Final-Std_LogSa_Final, 'b:', 'linewidth', 3);

h = legend([p1 p2 p7(1) p8(1) p4 p5], ...
           'Conditional mean spectrum', ...
           'Conditional mean +/- conditional \sigma', ...
           sprintf('%d pulse-type GMs',N_Pulse), ...
           sprintf('%d non-pulse-type GMs',N_NonPulse), ...
           sprintf('Mean of %d GMs',N_GM), ...
           sprintf('Mean of %d GMs +/- \\sigma',N_GM), ...
           'location','northeast');

    
xlabel('T (s)');
ylabel('{ln\it Sa}','Interpreter','latex');

xlim([min(T) max(T)])
ylim([floor(min(min(LogSa_Final)))-1 ceil(max(max(LogSa_Final)))])
xticks([T(1) max(T)*0.25:max(T)*0.25:max(T)])
yticks([-3:1:2])

a1=gca;
set(a1,'FontSize',32)
set(a1,'linewidth',2)
set(h,'FontSize',20)

saveas(Fig, 'lnSa of Selected GMs.tiff');


toc
