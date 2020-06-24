% This script plots exponential fits for time between burst events
close all;
clear;

folderPath = ['./Statistical files/'];

Qb5_all_night = load([folderPath,'all_night_induction_5_span_13.mat']);
Qb10_all_night = load([folderPath,'all_night_induction_10_span_13.mat']);
pp724_all_night = load([folderPath,'all_night_induction_24_span_13.mat']);
pp74_all_night = load([folderPath,'all_night_induction_4_span_13.mat']);

quies_duration = struct(...
    'PP74',pp74_all_night.time_stat_var,...
    'Qb5',Qb5_all_night.time_stat_var,...
    'Qb10',Qb10_all_night.time_stat_var,...
    'PP724',pp724_all_night.time_stat_var);
%% 4
[y4_cut,l4_cut] = histcounts(quies_duration.PP74(quies_duration.PP74>2),'BinWidth',1,'Normalization','Probability');
l4_cut = mean([l4_cut(1:end-1);l4_cut(2:end)]);

y4_cut_cut = y4_cut(y4_cut>0);
y4_cut_log = log(y4_cut_cut);
l4_cut_cut = l4_cut(y4_cut>0);

[f4,gof4] = fit(l4_cut_cut',y4_cut_log','poly1');

%% 24
[y24_cut,l24_cut] = histcounts(quies_duration.PP724(quies_duration.PP724>2),'BinWidth',1,'Normalization','Probability');
l24_cut = mean([l24_cut(1:end-1);l24_cut(2:end)]);

y24_cut_cut = y24_cut(y24_cut>0);
y24_cut_log = log(y24_cut_cut);
l24_cut_cut = l24_cut(y24_cut>0);

[f24,gof24] = fit(l24_cut_cut',y24_cut_log','poly1');
%% 5

[y5_cut,l5_cut] = histcounts(quies_duration.Qb5(quies_duration.Qb5>2),'BinWidth',1,'Normalization','Probability');
l5_cut = mean([l5_cut(1:end-1);l5_cut(2:end)]);

y5_cut_cut = y5_cut(y5_cut>0);
y5_cut_log = log(y5_cut_cut);
l5_cut_cut = l5_cut(y5_cut>0);

[f5,gof5] = fit(l5_cut_cut',y5_cut_log','poly1');

%% 10

[y10_cut,l10_cut] = histcounts(quies_duration.Qb10(quies_duration.Qb10>2),'BinWidth',1,'Normalization','Probability');
l10_cut = mean([l10_cut(1:end-1);l10_cut(2:end)]);

y10_cut_cut = y10_cut(y10_cut>0);
y10_cut_log = log(y10_cut_cut);
l10_cut_cut = l10_cut(y10_cut>0);

[f10,gof10] = fit(l10_cut_cut',y10_cut_log','poly1');

%% plot

figure;
subplot(4,1,1);
scatter(l4_cut_cut,y4_cut_log); hold on; plot(f4); 
subplot(4,1,2);
scatter(l5_cut_cut,y5_cut_log); hold on; plot(f5); title('5');
subplot(4,1,3);
scatter(l10_cut_cut,y10_cut_log); hold on; plot(f10); title('10');
subplot(4,1,4);
scatter(l24_cut_cut,y24_cut_log); hold on; plot(f24); title('24');


y4_plot = feval(f4,l10_cut_cut);
y24_plot = feval(f24,l10_cut_cut);
y5_plot = feval(f5,l10_cut_cut);
y10_plot = feval(f10,l10_cut_cut);


figure; 
subplot(2,4,1);
histogram(quies_duration.PP74,'BinWidth',1,'FaceColor','r','FaceAlpha',0.5,'Normalization','Probability');
xlim([1,50]); xlabel('Time [min]'); ylabel('Probability');
ylim([0,0.1]);
subplot(2,4,2);
histogram(quies_duration.Qb5,'BinWidth',1,'FaceColor','c','FaceAlpha',0.5,'Normalization','Probability');
xlim([1,50]); xlabel('Time [min]');
ylim([0,0.1]);
subplot(2,4,3);
histogram(quies_duration.Qb10,'BinWidth',1,'FaceColor','m','FaceAlpha',0.5,'Normalization','Probability');
xlim([1,50]); xlabel('Time [min]'); 
ylim([0,0.1]);
subplot(2,4,4); 
histogram(quies_duration.PP724,'BinWidth',1,'FaceColor','B','FaceAlpha',0.5,'Normalization','Probability'); 
xlabel('Time [min]'); 
xlim([1,50]); 
ylim([0,0.1]);

subplot(2,4,5:8);
plot(l10_cut_cut,y4_plot,'LineWidth',2,'Color','r'); hold on;
plot(l10_cut_cut,y5_plot,'LineWidth',2,'Color','c');
plot(l10_cut_cut,y10_plot,'LineWidth',2,'Color','m');
plot(l10_cut_cut,y24_plot,'LineWidth',2,'Color','B');
legend('PP7x4','Qbx5','Qbx10','PP7x24');
dim = [0.15 0.1 0.3 0.2];
str = {['PP7x4: slope=',num2str(f4.p1,2) ,' rsquare= ',num2str(gof4.rsquare,2)],...
    ['Qb5: slope=',num2str(f5.p1,2) ,' rsquare= ',num2str(gof5.rsquare,2)],...
    ['Qbx10: slope=',num2str(f10.p1,2) ,' rsquare= ',num2str(gof10.rsquare,2)],...
    ['PP7x24: slope=',num2str(f24.p1,2) ,' rsquare= ',num2str(gof24.rsquare,2)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
xlabel('Time [min]'); ylabel('log Probability');


slopes = abs([1/f4.p1,1/f5.p1,1/f10.p1,1/f24.p1]);
e4 = confint(f4,0.95);
e5 = confint(f5,0.95);
e10 = confint(f10,0.95);
e24 = confint(f24,0.95);
% need to properly propegate error to take the reciprocal
error_pos = [abs(abs(f4.p1-e4(1,1))/f4.p1),...
    abs(abs(f5.p1-e5(1,1))/f5.p1),...
    abs(abs(f10.p1-e10(1,1))/f10.p1),...
    abs(abs(f24.p1-e24(1,1))/f24.p1)];
error_neg = [abs(abs(f4.p1-e4(2,1))/f4.p1),...
    abs(abs(f5.p1-e5(2,1))/f5.p1),...
    abs(abs(f10.p1-e10(2,1))/f10.p1),...
    abs(abs(f24.p1-e24(2,1))/f24.p1)];

error_pos = [abs(e4(1,1)/f4.p1),...
    abs(e5(1,1)/f5.p1),...
    abs(e10(1,1)/f10.p1),...
    abs(e24(1,1)/f24.p1)];
error_neg = [abs(e4(2,1)/f4.p1),...
    abs(e5(2,1)/f5.p1),...
    abs(e10(2,1)/f10.p1),...
    abs(e24(2,1)/f24.p1)];

figure;
bar(1:4,slopes); hold on;
errorbar(1:4,slopes,error_neg,error_pos,'.');
xlim([0,5]);
xticklabels({'PP7-4x','Qb-5x','Qb-10x','PP7-24x'});
ylabel('Estimated burst lag time [min]');