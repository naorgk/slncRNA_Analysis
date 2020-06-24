% This script plots comparisons between statistical data gathered from
% different experiments.

clear; close all;
options = statset('MaxIter',1000);

folderPath = './Statistical files/';

Qb5_all_night = load([folderPath,'all_night_induction_5_span_13.mat']);
Qb10_all_night = load([folderPath,'all_night_induction_10_span_13.mat']);
pp724_all_night = load([folderPath,'all_night_induction_24_span_13.mat']);
pp74_all_night = load([folderPath,'all_night_induction_4_span_13.mat']);

pp74_multicopy_intensity = load([folderPath,'all_night_induction_4pp7_multicopy_data.mat']);
pp74_multicopy_intensity = cell2mat(pp74_multicopy_intensity.intensity);
pp74_multicopy_background = load([folderPath,'all_night_induction_4pp7_multicopy_norm.mat']);
pp74_multicopy_background = cell2mat(pp74_multicopy_background.background);
pp74_multicopy_intensity = pp74_multicopy_intensity(pp74_multicopy_background>3500);
pp74_multicopy_background = pp74_multicopy_background(pp74_multicopy_background>3500);
pp74_multicopy_intensity = pp74_multicopy_intensity-pp74_multicopy_background;
pp74_multicopy_background = struct('background_mean',pp74_multicopy_background);
pp74_multicopy_intensity = struct('meanmRNA',pp74_multicopy_intensity);

pp7_control_background = load([folderPath,'all_night_induction_PP7_RBP_only.mat']);
pp7_control_background = double(pp7_control_background.pixels);

pos_amp_lim = 2000;

all_night_amp_pos = struct(...
    'PP74',pp74_all_night.amp_pos_total(pp74_all_night.amp_pos_total<pos_amp_lim),...
    'Qb5',Qb5_all_night.amp_pos_total(Qb5_all_night.amp_pos_total<pos_amp_lim),...
    'Qb10',Qb10_all_night.amp_pos_total(Qb10_all_night.amp_pos_total<pos_amp_lim),...
    'PP724',pp724_all_night.amp_pos_total(pp724_all_night.amp_pos_total<pos_amp_lim));

var_amp_lim = 1000;

all_night_amp_var = struct(...
    'PP74',pp74_all_night.amp_var_total(abs(pp74_all_night.amp_var_total)<var_amp_lim),...
    'Qb5',Qb5_all_night.amp_var_total(abs(Qb5_all_night.amp_var_total)<var_amp_lim),...
    'Qb10',Qb10_all_night.amp_var_total(abs(Qb10_all_night.amp_var_total)<var_amp_lim),...
    'PP724',pp724_all_night.amp_var_total(abs(pp724_all_night.amp_var_total)<var_amp_lim));

neg_amp_lim = 2000;

all_night_amp_neg = struct(...
    'PP74',pp74_all_night.amp_neg_total(abs(pp74_all_night.amp_neg_total)<neg_amp_lim),...
    'Qb5',Qb5_all_night.amp_neg_total(abs(Qb5_all_night.amp_neg_total)<neg_amp_lim),...
    'Qb10',Qb10_all_night.amp_neg_total(abs(Qb10_all_night.amp_neg_total)<neg_amp_lim),...
    'PP724',pp724_all_night.amp_neg_total(abs(pp724_all_night.amp_neg_total)<neg_amp_lim));


all_night_time_pos = struct(...
    'PP74',pp74_all_night.time_stat_pos,...
    'Qb5',Qb5_all_night.time_stat_pos,...
    'Qb10',Qb10_all_night.time_stat_pos,...
    'PP724',pp724_all_night.time_stat_pos);

all_night_time_var = struct(...
    'PP74',pp74_all_night.time_stat_var,...
    'Qb5',Qb5_all_night.time_stat_var,...
    'Qb10',Qb10_all_night.time_stat_var,...
    'PP724',pp724_all_night.time_stat_var);

all_night_time_neg = struct(...
    'PP74',pp74_all_night.time_stat_neg,...
    'Qb5',Qb5_all_night.time_stat_neg,...
    'Qb10',Qb10_all_night.time_stat_neg,...
    'PP724',pp724_all_night.time_stat_neg);

all_night_time_between_pos = struct(...
    'PP74',pp74_all_night.time_between_pos,...
    'Qb5',Qb5_all_night.time_between_pos,...
    'Qb10',Qb10_all_night.time_between_pos,...
    'PP724',pp724_all_night.time_between_pos);


all_night_meanmRNA = struct(...
    'PP74',pp74_all_night.meanmRNA,...
    'Qb5',Qb5_all_night.meanmRNA,...
    'Qb10',Qb10_all_night.meanmRNA,...
    'PP724',pp724_all_night.meanmRNA,...
    'PP74_multi',pp74_multicopy_intensity.meanmRNA);

all_night_background_mean = struct(...
    'PP74',pp74_all_night.background_mean,...
    'Qb5',Qb5_all_night.background_mean,...
    'Qb10',Qb10_all_night.background_mean,...
    'PP724',pp724_all_night.background_mean,...
    'PP74_multi',pp74_multicopy_background.background_mean);

all_night_postVar = struct(...
    'PP74',pp74_all_night.postVar,...
    'Qb5',Qb5_all_night.postVar,...
    'Qb10',Qb10_all_night.postVar,...
    'PP724',pp724_all_night.postVar);


%% Plots


figure;
subplot(2,1,1)
histogram(all_night_amp_pos.PP724,'BinWidth',50,'Normalization','probability'); hold on;
histogram(all_night_amp_pos.Qb10,'BinWidth',50,'FaceColor','m','Normalization','probability');
histogram(all_night_amp_pos.Qb5,'BinWidth',50,'Normalization','probability');
histogram(all_night_amp_pos.PP74,'BinWidth',50,'Normalization','probability');
title('Positive amplitudes'); xlabel('Amp [A.U]'); ylabel('Frequency');
legend('PP7 24', 'Qb 10','Qb 5','PP7 4'); ylim([0,0.15]);
 
subplot(2,1,2);
histogram(all_night_amp_neg.PP724,'BinWidth',50,'Normalization','probability'); hold on;
histogram(all_night_amp_neg.Qb10,'BinWidth',50,'FaceColor','m','Normalization','probability');
histogram(all_night_amp_neg.Qb5,'BinWidth',50,'Normalization','probability');
histogram(all_night_amp_neg.PP74,'BinWidth',50,'Normalization','probability');
legend('PP7 24', 'Qb 10','Qb 5','PP7 4','Location','NorthWest');
ylabel('Frequency'); title('Negative amplitudes'); xlabel('Amp [A.U]');

figure; 
histogram(all_night_amp_var.PP724,'BinWidth',50,'Normalization','probability'); hold on;
histogram(all_night_amp_var.Qb10,'BinWidth',50,'FaceColor','m','Normalization','probability');
histogram(all_night_amp_var.Qb5,'BinWidth',50,'Normalization','probability');
histogram(all_night_amp_var.PP74,'BinWidth',50,'Normalization','probability');
title('Quiescent amplitudes'); xlabel('Amp [A.U]');
legend('PP7 24', 'Qb 10','Qb 5','PP7 4');
xlim([-1000,1000]);

figure;
histogram(all_night_time_neg.PP724,'BinWidth',0.1,'Normalization','probability'); hold on;
histogram(all_night_time_neg.Qb10,'FaceColor','m','BinWidth',0.1,'Normalization','probability');
histogram(all_night_time_neg.Qb5,'BinWidth',0.1,'Normalization','probability');
histogram(all_night_time_neg.PP74,'BinWidth',0.1,'Normalization','probability');
title('Duration of negative increments'); xlabel('Time [min]');
legend('PP7 24', 'Qb 10','Qb 5','PP7 4');
xlim([1,5]);

figure;
histogram(all_night_time_var.PP724,'BinWidth',0.5,'Normalization','probability'); hold on;
histogram(all_night_time_var.Qb10,'FaceColor','m','BinWidth',0.5,'Normalization','probability');
histogram(all_night_time_var.Qb5,'BinWidth',0.5,'Normalization','probability');
histogram(all_night_time_var.PP74,'BinWidth',0.5,'Normalization','probability');
title('Duration of quiescent increments'); xlabel('Time [min]');
legend('PP7 24', 'Qb 10','Qb 5','PP7 4');
xlim([1,40]);

figure;
histogram(all_night_time_between_pos.PP724,'BinWidth',0.5,'Normalization','probability'); hold on;
histogram(all_night_time_between_pos.Qb10,'FaceColor','m','BinWidth',0.5,'Normalization','probability');
histogram(all_night_time_between_pos.Qb5,'BinWidth',0.5,'Normalization','probability');
histogram(all_night_time_between_pos.PP74,'BinWidth',0.5,'Normalization','probability');
title('Duration of quiescent increments'); xlabel('Time [min]');
legend('PP7 24', 'Qb 10','Qb 5','PP7 4');
xlim([1,40]);

figure;
histogram(all_night_time_pos.PP724,'BinWidth',0.1,'Normalization','probability'); hold on;
histogram(all_night_time_pos.Qb10,'FaceColor','m','BinWidth',0.1,'Normalization','probability');
histogram(all_night_time_pos.Qb5,'BinWidth',0.1,'Normalization','probability');
histogram(all_night_time_pos.PP74,'BinWidth',0.1,'Normalization','probability');
title('Duration of positive increments'); xlabel('Time [min]');
legend('PP7 24', 'Qb 10','Qb 5','PP7 4');
xlim([1,5]);

figure;
binwidth = 500;
subplot(4,1,1);
histogram(all_night_background_mean.PP724,'BinWidth',binwidth,'Normalization','probability'); hold on;
histogram(all_night_meanmRNA.PP724,'BinWidth',binwidth,'Normalization','probability'); 
legend('Background mean','Spot mean'); title('PP7x24'); xlim([0,3e4]);

subplot(4,1,2);
histogram(all_night_background_mean.PP74,'BinWidth',binwidth,'Normalization','probability'); hold on;
histogram(all_night_meanmRNA.PP74,'BinWidth',binwidth,'Normalization','probability'); 
legend('Background mean','Spot mean'); title('PP7x4'); xlim([0,3e4]);
subplot(4,1,3);
histogram(all_night_background_mean.Qb5,'BinWidth',binwidth,'Normalization','probability'); hold on;
histogram(all_night_meanmRNA.Qb5,'BinWidth',binwidth,'Normalization','probability'); 
legend('Background mean','Spot mean'); title('Qbx5'); xlim([0,3e4]);

subplot(4,1,4);
histogram(all_night_background_mean.Qb10,'BinWidth',binwidth,'Normalization','probability'); hold on;
histogram(all_night_meanmRNA.Qb10,'BinWidth',binwidth,'Normalization','probability'); 
legend('Background mean','Spot mean'); title('Qbx10'); xlim([0,3e4]);


violin_data = [{pp7_control_background},{all_night_background_mean.PP74},{all_night_background_mean.Qb5},...
    {all_night_background_mean.Qb10},{all_night_background_mean.PP724}];
figure;
violin(violin_data);
xticklabels([{'PP7 RBP only'},{'PP7x4'},{'Qbx5'},{'Qbx10'},{'PP7x24'}]);

binWidth = 2;
xlimit = 50;
figure;
subplot(4,1,1);
histogram(all_night_meanmRNA.PP724/378,'BinWidth',binWidth); 
xlim([0,xlimit]); title('PP7 x 24'); ylabel('Count');
subplot(4,1,2);
histogram(all_night_meanmRNA.Qb10/459,'BinWidth',binWidth);
xlim([0,xlimit]); title('Qb x 10'); ylabel('Count');
subplot(4,1,3);
histogram(all_night_meanmRNA.Qb5/325,'BinWidth',binWidth);
xlim([0,xlimit]); title('Qb x 5'); ylabel('Count');
subplot(4,1,4);
histogram(all_night_meanmRNA.PP74/171.5,'BinWidth',binWidth);
xlim([0,xlimit]); title('PP7 x 4'); ylabel('Count'); xlabel('# Scaffolds');


binwidth = 200;
figure;
histogram(all_night_background_mean.PP74,'BinWidth',binwidth,'Normalization','probability');
hold on;
histogram(all_night_background_mean.PP74_multi,'BinWidth',binwidth,'Normalization','probability');
title('Mean of backgrounds'); xlabel('Intensity [A.U]'); ylabel('probability');
legend('PP7 4','PP7 4 multicopy');

binwidth = 200;
figure;
histogram(all_night_meanmRNA.PP74,'BinWidth',binwidth,'Normalization','probability');
hold on;
histogram(all_night_meanmRNA.PP74_multi,'BinWidth',binwidth,'Normalization','probability');
title('Mean of spot intensity'); xlabel('Intensity [A.U]'); ylabel('probability');
legend('PP7 4','PP7 4 multicopy');

binwidth = 2;
figure;
histogram(all_night_meanmRNA.PP74/171.5,'BinWidth',binwidth,'Normalization','probability');
hold on;
histogram(all_night_meanmRNA.PP74_multi/171.5,'BinWidth',binwidth,'Normalization','probability');
title('estimated scaffold number'); xlabel('# scaffolds'); ylabel('probability');
legend('PP7 4','PP7 4 multicopy');


figure;
categories = [{'Negative'},{'Quiescent'},{'Positive'}];
subplot(2,2,1);
histogram(all_night_postVar.PP724); 
title('PP7 x 24'); ylabel('Count'); xticks([-1,0,1]); xticklabels(categories); ylabel('Count'); 
subplot(2,2,2);
histogram(all_night_postVar.Qb10);
title('Qb x 10'); ylabel('Count'); xticks([-1,0,1]); xticklabels(categories); ylabel('Count'); 
subplot(2,2,3); 
histogram(all_night_postVar.Qb5);
title('Qb x 5'); ylabel('Count'); xticks([-1,0,1]); xticklabels(categories); ylabel('Count'); 
subplot(2,2,4);
histogram(all_night_postVar.PP74); 
title('PP7 x 4'); xticks([-1,0,1]); xticklabels(categories); ylabel('Count'); 



