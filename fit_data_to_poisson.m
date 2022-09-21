%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads fluorescence intensity measurements (for in vitro) or
% burst amplitude distributions (in vivo) and fits the measuremnts data to
% a modified Poisson distribution.
% The fitting is done on 2 parameters: 
% - lambda (Poisson rate)
% - K0 (normalization factor - can also be considered as the fluorescence
%   of a single protein-bound slncRNA) 

% Principle of operation:
% The fluorescence intensity data is first binned into bins of equal width to
% construct a counting process. 
% For each lambda value, a theoretical Poisson distribution is generated to
% compare against. Then, for each K0 value, the fluorescence data and the bin
% edges are normalized by K0 to get an x-axis on which a Poisson
% distribution is defined (i.e., [0-10]). The mean squared error (MSE) between
% the theoretical Poisson distribution and the experimental, normalized
% counts is computed. 
% The output is the lambda and K0 values which minimize the MSE to the
% theoretical Poisson distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In vivo amplitude data
clear; close all;

PCP_24_data = load('./Data for Poisson fit/PCP_24x.mat');
PCP_4_data = load('./Data for Poisson fit/PCP_4x_QCP_5x.mat');

% Iterate through the following lambda,K0 values
lambda = 0:10;
K0 = 10:1:2000;

% PCP-4x in-vivo amplitude fitting 
[hist_x_axis4pos,nn_full4pos,yq_out4pos,nn4pos] = poisson_sub_func(PCP_4_data.amp_pos_total,lambda,K0,10,'PP7 4 positive');
[hist_x_axis4neg,nn_full4neg,yq_out4neg,nn4neg] = poisson_sub_func(PCP_4_data.amp_neg_total,lambda,K0,20,'PP7 4 negative');

% PCP-24x in-vivo amplitude fitting 
[hist_x_axis24pos,nn_full24pos,yq_out24pos,nn24pos] = poisson_sub_func(PCP_24_data.amp_pos_total,lambda,K0,75,'PP7 24 positive');
% Remove outlier from data
PP7_24_neg_amps_filtered = PCP_24_data.amp_neg_total(PCP_24_data.amp_neg_total>-6000);
[hist_x_axis24neg,nn_full24neg,yq_out24neg,nn24neg] = poisson_sub_func(PP7_24_neg_amps_filtered,lambda,K0,75,'PP7 24 negative');

%% In vitro mean fluorescence data
% Estimated uracil quantity in slncRNA sequences
% Intensity data is noramlized by these numbers to get a standardaized
% measurement and to reduce noise.
% 3x - 56
% 3x/3x - 91
% 4x - 86
% 4x/4x - 117
% 8x - 122
% 14x/15x - 326

% Iterate through the following lambda,K0 values
lambda = 0:10;
K0 = 0.1:0.1:50;

median_measurements_4x = xlsread('./Data for Poisson fit/median_measurements_4.csv');
[hist_x_axis4,nn_full4,yq_out4,nn4] = poisson_sub_func(median_measurements_4x(:,4)/86,lambda,K0,0.5,'PP7-4x-RNA-only granule fluorescence');

% lambda = 0 appears to be a default solution to several measurements
% despite not fitting the data optimally (by manual inspection). Therefore
% the 0 option is removed
lambda = 1:10; 
median_measurements_8x = xlsread('./Data for Poisson fit/median_measurements_8.csv');
[hist_x_axis8,nn_full8,yq_out8,nn8] = poisson_sub_func(median_measurements_8x(:,4)/122,lambda,K0,0.5,'PP7-8x-RNA-only granule fluorescence');

median_measurements_3x_3x = xlsread('./Data for Poisson fit/median_measurements_3_3.csv');
[hist_x_axis33,nn_full33,yq_out33,nn33] = poisson_sub_func(median_measurements_3x_3x(:,4)/91,lambda,K0,0.5,'PP7-3x-MS2-3x-RNA-only granule fluorescence');

median_measurements_4x_4x = xlsread('./Data for Poisson fit/median_measurements_4_4.csv');
[hist_x_axis44,nn_full44,yq_out44,nn44] = poisson_sub_func(median_measurements_4x_4x(:,4)/117,lambda,K0,0.5,'PP7-4x-MS2-4x-RNA-only granule fluorescence');

median_measurements_14x_15x = xlsread('./Data for Poisson fit/median_measurements_14_15.csv');
[hist_x_axis14,nn_full14,yq_out14,nn14] = poisson_sub_func(median_measurements_14x_15x(:,4)/326,lambda,K0,0.5,'PP7-14x-MS2-15x-RNA-only granule fluorescence');


%%%% Generate Figure 1 plots %%%%

% Violin plot of median granule intensity values
figure;
violin([{(median_measurements_4x(:,4))/(86*0.35)},...
    {(median_measurements_3x_3x(:,4))/(91*0.35)},...
    {(median_measurements_4x_4x(:,4))/(117*0.35)},...
    {(median_measurements_8x(:,4))/(122*0.35)},...
    {(median_measurements_14x_15x(:,4))/(326*0.35)}]); 

% Poisson function fits
figure; 
subplot(1,5,1)
scatter(hist_x_axis4,nn_full4);  
hold on;
plot(hist_x_axis4,yq_out4*max(nn4)./max(yq_out4),'r*-');
box on;
title('4'); xlim([1,10]);
xlabel('Intensity [A.U]'); ylabel('Probability');
subplot(1,5,2)
scatter(hist_x_axis33,nn_full33);  
hold on;
plot(hist_x_axis33,yq_out33*max(nn33)./max(yq_out33),'r*-');
box on;
title('3/3'); xlim([1,10]);
xlabel('Intensity [A.U]'); ylabel('Probability');
subplot(1,5,4)
scatter(hist_x_axis8,nn_full8);  
hold on;
plot(hist_x_axis8,yq_out8*max(nn8)./max(yq_out8),'r*-');
box on;
title('8'); xlim([1,10]);
xlabel('Intensity [A.U]'); ylabel('Probability');
subplot(1,5,3)
scatter(hist_x_axis44,nn_full44);  
hold on;
plot(hist_x_axis44,yq_out44*max(nn44)./max(yq_out44),'r*-');
box on;
title('4/4'); xlim([1,10]);
xlabel('Intensity [A.U]'); ylabel('Probability');
subplot(1,5,5)
scatter(hist_x_axis14,nn_full14);  
hold on;
plot(hist_x_axis14,yq_out14*max(nn14)./max(yq_out14),'r*-');
box on;
title('14/15'); xlim([1,10]);
xlabel('Intensity [A.U]'); ylabel('Probability');



% K0 values scatter plot
figure;
scatter(1,5.7,'r');
hold on;
scatter([2,3],[3.1,1.7],'b');
scatter(4,1.4,'k');
scatter(5,0.8,'k');
xlim([0,6]);
ylabel('K0'); xlabel('slncRNA sequence');
