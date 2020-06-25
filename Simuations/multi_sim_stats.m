options = statset('MaxIter',1000);
const = load('sim_const.mat');
slope = load('sim_slope.mat');
burst = load('sim_burst.mat');

const_amp_pos = const.amp_pos(const.amp_pos<1000);
slope_amp_pos = slope.amp_pos(slope.amp_pos<1000);
burst_amp_pos = burst.amp_pos(burst.amp_pos<1000);

const_amp_var = const.amp_var(abs(const.amp_var)<1000);
slope_amp_var = slope.amp_var(abs(slope.amp_var)<1000);
burst_amp_var = burst.amp_var(abs(burst.amp_var)<1000);

const_amp_neg = const.amp_neg(abs(const.amp_neg)<1000);
slope_amp_neg = slope.amp_neg(abs(slope.amp_neg)<1000);
burst_amp_neg = burst.amp_neg(abs(burst.amp_neg)<1000);

const_time_pos = const.time_stat_pos;
slope_time_pos = slope.time_stat_pos;
burst_time_pos = burst.time_stat_pos;

const_time_neg = const.time_stat_neg;
slope_time_neg = slope.time_stat_neg;
burst_time_neg = burst.time_stat_neg;

const_time_var = const.time_stat_var;
slope_time_var = slope.time_stat_var;
burst_time_var = burst.time_stat_var;

const_rate_pos = const.amp_pos(const.amp_pos<1000)./const.time_stat_pos(const.amp_pos<1000);
slope_rate_pos = slope.amp_pos(slope.amp_pos<1000)./slope.time_stat_pos(slope.amp_pos<1000);
burst_rate_pos = burst.amp_pos(burst.amp_pos<1000)./burst.time_stat_pos(burst.amp_pos<1000);

const_rate_neg = const.amp_neg(abs(const.amp_neg)<1000)./const.time_stat_neg(abs(const.amp_neg)<1000);
slope_rate_neg = slope.amp_neg(abs(slope.amp_neg)<1000)./slope.time_stat_neg(abs(slope.amp_neg)<1000);
burst_rate_neg = burst.amp_neg(abs(burst.amp_neg)<1000)./burst.time_stat_neg(abs(burst.amp_neg)<1000);

const_rate_var = const.amp_var(abs(const.amp_var)<1000)./const.time_stat_var(abs(const.amp_var)<1000);
slope_rate_var = slope.amp_var(abs(slope.amp_var)<1000)./slope.time_stat_var(abs(slope.amp_var)<1000);
burst_rate_var = burst.amp_var(abs(burst.amp_var)<1000)./burst.time_stat_var(abs(burst.amp_var)<1000);


%% show all in same figure
figure; 
histogram(const_amp_pos,'BinWidth',10,'FaceColor','k');hold on; 
histogram(slope_amp_pos,'BinWidth',10,'FaceColor','m'); 
histogram(burst_amp_pos,'BinWidth',10,'FaceColor','c');
title('Positive amplitudes'); xlabel('Amp [A.U]'); ylabel('Probability');
legend('const','slope', 'bursts');

figure; 
histogram(const_amp_neg,'BinWidth',10,'FaceColor','k'); hold on;
histogram(slope_amp_neg,'BinWidth',10,'FaceColor','m'); 
histogram(burst_amp_neg,'BinWidth',10,'FaceColor','c');
legend('const','slope', 'bursts');
title('Negative amplitudes'); xlabel('Amp [A.U]'); ylabel('Probability');


figure;
subplot(2,2,1);
histogram(const_time_pos,'FaceColor','k','BinWidth',0.1);
xlim([0,5]); ylim([0,200]);
ylabel('Counts');
title('Positive durations'); xlabel('Time [min]');
subplot(2,2,[3,4]);
histogram(const_time_var,'FaceColor','m','BinWidth',0.5);
xlim([0,45]); ylabel('Counts');
title('Quiescent durations'); xlabel('Time [min]'); 
subplot(2,2,2);
histogram(const_time_neg,'FaceColor','c','BinWidth',0.1); 
xlim([0,5]); ylim([0,200]); ylabel('Counts');
title('Negative durations'); xlabel('Time [min]');
suptitle('Constant signal times'); 

figure;
subplot(2,2,1);
histogram(slope_time_pos,'FaceColor','k','BinWidth',0.1);
xlim([0,5]); ylim([0,8]);
ylabel('Counts');
title('Positive durations'); xlabel('Time [min]');
subplot(2,2,[3,4]);
histogram(slope_time_var,'FaceColor','m','BinWidth',0.5);
xlim([0,45]); ylabel('Counts');
title('Quiescent durations'); xlabel('Time [min]'); 
subplot(2,2,2);
histogram(slope_time_neg,'FaceColor','c','BinWidth',0.1); 
xlim([0,5]); ylim([0,8]); ylabel('Counts'); 
title('Negative durations'); xlabel('Time [min]');
suptitle('slope signal times'); 

figure;
subplot(2,2,1);
histogram(burst_time_pos,'FaceColor','k','BinWidth',0.1);
xlim([0,5]); 
ylabel('Counts');
title('Positive durations'); xlabel('Time [min]');
subplot(2,2,[3,4]);
histogram(burst_time_var,'FaceColor','m','BinWidth',0.5);
xlim([0,45]); ylabel('Counts');
title('Quiescent durations'); xlabel('Time [min]'); 
subplot(2,2,2);
histogram(burst_time_neg,'FaceColor','c','BinWidth',0.1); 
xlim([0,5]); ylabel('Counts');
title('Negative durations'); xlabel('Time [min]');
suptitle('burst signal times'); 

figure;
histogram(const_time_neg,'Normalization','Probability','FaceColor','k');
hold on; histogram(slope_time_neg,'Normalization','Probability','FaceColor','m');
histogram(burst_time_neg,'Normalization','Probability','FaceColor','c');
title('Lengths of negative increments'); xlabel('Time [min]');
legend('const','slope', 'burst');ylabel('Probability');

figure; 
histogram(const_rate_pos,'BinWidth',10,'Normalization','Probability','FaceColor','k'); hold on; 
ylabel('Probability'); 
histogram(slope_rate_pos,'BinWidth',10,'Normalization','Probability','FaceColor','m');
histogram(burst_rate_pos,'BinWidth',10,'Normalization','Probability','FaceColor','c');
title('Positive rates'); xlabel('Rate [A.U/min]'); legend('const','slope', 'burst');

figure; 
histogram(const_rate_neg,'BinWidth',10,'Normalization','Probability','FaceColor','k'); hold on; 
ylabel('Probability');
histogram(slope_rate_neg,'BinWidth',10,'Normalization','Probability','FaceColor','m');
histogram(burst_rate_neg,'BinWidth',10,'Normalization','Probability','FaceColor','c');
title('Negative rates'); xlabel('Rate [A.U/min]'); legend('const','slope', 'burst');


%% show burst with slope

figure; 
histogram(slope_amp_pos,'BinWidth',10,'FaceColor','m'); hold on;
histogram(burst_amp_pos,'BinWidth',10,'FaceColor','c');
title('Positive amplitudes'); xlabel('Amp [A.U]'); ylabel('Probability');
legend('slope', 'burst'); 

figure;
gmm_x = 0:0.5:500;
gmm = fitgmdist(slope_amp_pos(abs(slope_amp_pos)<500),5,'Options',options);
gmm_y = pdf(gmm,gmm_x');
plot(gmm_x,gmm_y,'m','lineWidth',2); hold on;
gmm = fitgmdist(burst_amp_pos(abs(burst_amp_pos)<500),5,'Options',options);
gmm_y = pdf(gmm,gmm_x');
plot(gmm_x,gmm_y,'c','lineWidth',2);
title('Multi Gaussian fitting'); xlabel('Intensity [A.U]'); ylabel('Probability');
legend('slope','burst');

figure; 
histogram(slope_amp_neg,'BinWidth',10,'FaceColor','m'); hold on;
histogram(burst_amp_neg,'BinWidth',10,'FaceColor','c');
legend('slope', 'burst'); ylabel('Probability');
title('Negative amplitudes'); xlabel('Amp [A.U]');

figure;
gmm_x = -500:0.5:0;
gmm = fitgmdist(slope_amp_neg(abs(slope_amp_neg)<500),5,'Options',options);
gmm_y = pdf(gmm,gmm_x');
plot(gmm_x,gmm_y,'m','lineWidth',2); hold on;
gmm = fitgmdist(burst_amp_neg(abs(burst_amp_neg)<500),5,'Options',options);
gmm_y = pdf(gmm,gmm_x');
plot(gmm_x,gmm_y,'c','lineWidth',2);
title('Multi Gaussian fitting'); xlabel('Intensity [A.U]'); ylabel('Probability');
legend('slope','burst');

figure; 
histogram(slope_amp_var,'BinWidth',10,'Normalization','Probability','FaceColor','m'); hold on;
histogram(burst_amp_var,'BinWidth',10,'Normalization','Probability','FaceColor','c');
title('Quiescent amplitudes'); xlabel('Amp [A.U]');
legend('slope', 'burst'); ylabel('Probability'); 

figure;
hold on; histogram(slope_time_pos,'BinWidth',0.1,'Normalization','Probability','FaceColor','m');
histogram(burst_time_pos,'BinWidth',0.1,'Normalization','Probability','FaceColor','c');
legend('slope', 'burst'); ylabel('Probability');
title('Lengths of positive increments'); xlabel('Time [min]');

figure;
hold on; histogram(const_time_var,'BinWidth',1,'FaceColor','m');
histogram(burst_time_var,'BinWidth',1,'FaceColor','c');
legend('const', 'burst'); ylabel('Probability');
title('Lengths of quiescent increments'); xlabel('Time [min]');

figure;
histogram(slope_time_neg(slope_time_neg<10),'BinWidth',0.1,'FaceColor','m'); hold on;
histogram(burst_time_neg,'BinWidth',0.1,'FaceColor','c');
title('Lengths of negative increments'); xlabel('Time [min]');
legend('slope', 'burst'); ylabel('Probability');

figure; 
histogram(slope_rate_pos,'BinWidth',10,'Normalization','Probability','FaceColor','m'); hold on;
histogram(burst_rate_pos,'BinWidth',10,'Normalization','Probability','FaceColor','c');
title('Positive rates'); xlabel('Rate [A.U/min]'); legend('slope', 'burst'); ylabel('Probability');

figure; 
histogram(slope_rate_var(abs(slope_rate_var)<150),'BinWidth',10,'Normalization','Probability','FaceColor','m'); hold on;
histogram(burst_rate_var(abs(burst_rate_var)<150),'BinWidth',10,'Normalization','Probability','FaceColor','c');
title('Quiescent rates'); xlabel('Rate [A.U/min]'); 
ylabel('Probability'); legend('slope', 'burst'); xlim([-150,150]);

figure; 
histogram(slope_rate_neg,'BinWidth',10,'Normalization','Probability','FaceColor','m'); hold on;
histogram(burst_rate_neg,'BinWidth',10,'Normalization','Probability','FaceColor','c');
title('Negative rates'); xlabel('Rate [A.U/min]'); legend('24 pp7', '5 Qb');
ylabel('Probability');

%% show burst with const

figure; 
histogram(const_amp_pos,'BinWidth',10); hold on;
histogram(burst_amp_pos,'BinWidth',10,'FaceColor','r');
title('Positive amplitudes'); xlabel('Amp [A.U]');
legend('const', 'burst');

figure; 
histogram(const_amp_neg,'BinWidth',10); hold on;
histogram(burst_amp_neg,'BinWidth',10,'FaceColor','r');
legend('const', 'burst');
title('Negative amplitudes'); xlabel('Amp [A.U]');


figure;
hold on; histogram(const_time_pos);
histogram(burst_time_pos,'FaceColor','r');
legend('const', 'burst');
title('Lengths of positive increments'); xlabel('Time [min]');

figure;
histogram(const_time_neg); hold on;
histogram(burst_time_neg,'FaceColor','r');
title('Lengths of negative increments'); xlabel('Time [min]');
legend('const', 'burst');

figure; 
histogram(const_rate_pos,'BinWidth',10,'Normalization','Probability'); hold on;
histogram(burst_rate_pos,'BinWidth',10,'Normalization','Probability','FaceColor','r');
title('Positive rates'); xlabel('Rate [A.U/min]'); legend('const', 'burst');

figure; 
histogram(const_rate_neg,'BinWidth',10,'Normalization','Probability'); hold on;
histogram(burst_rate_neg,'BinWidth',10,'Normalization','Probability','FaceColor','r');
title('Negative rates'); xlabel('Rate [A.U/min]'); legend('const', 'burst');


%% show 24 pp7 with 4 pp7

figure; 
histogram(const_amp_pos,'BinWidth',10,'Normalization','Probability'); hold on;
histogram(pp7_24_amp_pos,'BinWidth',10,'Normalization','Probability','FaceColor','g');
title('Positive amplitudes'); xlabel('Amp [A.U]');
legend('4 pp7', '24 pp7');

figure; 
histogram(const_amp_neg,'BinWidth',10,'Normalization','Probability'); hold on;
histogram(pp7_24_amp_neg,'BinWidth',10,'Normalization','Probability','FaceColor','g'); 
legend('4 pp7', '24 pp7');
title('Negative amplitudes'); xlabel('Amp [A.U]');


figure;
hold on; histogram(const_time_pos,'Normalization','Probability');
histogram(pp7_24_time_neg,'Normalization','Probability','FaceColor','g');
legend('4 pp7', '24 pp7');
title('Lengths of positive increments'); xlabel('Time [min]');

figure;
histogram(const_time_neg,'Normalization','Probability'); hold on;
histogram(pp7_24_time_neg,'Normalization','Probability','FaceColor','g'); hold on;
title('Lengths of negative increments'); xlabel('Time [min]');
legend('4 pp7', '24 pp7');

figure; 
histogram(const_rate_pos,'BinWidth',10,'Normalization','Probability'); hold on;
histogram(pp7_24_rate_pos,'BinWidth',10,'Normalization','Probability','FaceColor','g'); hold on;
title('Positive rates'); xlabel('Rate [A.U/min]'); legend('4 pp7', '24 pp7');

figure; 
histogram(const_rate_neg,'BinWidth',10,'Normalization','Probability'); hold on;
histogram(pp7_24_rate_neg,'BinWidth',10,'Normalization','Probability','FaceColor','g');
title('Negative rates'); xlabel('Rate [A.U/min]'); legend('4 pp7', '24 pp7');
