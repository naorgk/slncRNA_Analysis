% This script plot Poisson fits for amplitude data

function [x,y] = poisson_fitting_plots()
clear;
close all;

folderPath = ['./Statistical files/'];

Qb5_all_night = load([folderPath,'all_night_induction_5_span_13.mat']);
Qb10_all_night = load([folderPath,'all_night_induction_10_span_13.mat']);
pp7_24_all_night = load([folderPath,'all_night_induction_24_span_13.mat']);
pp7_4_all_night = load([folderPath,'all_night_induction_4_span_13.mat']);


lim = 4000;
Qb5_all_night_amp_pos = Qb5_all_night.amp_pos_total(Qb5_all_night.amp_pos_total<lim);
Qb10_all_night_amp_pos = Qb10_all_night.amp_pos_total(Qb10_all_night.amp_pos_total<lim);
pp7_24_all_night_amp_pos = pp7_24_all_night.amp_pos_total(pp7_24_all_night.amp_pos_total<lim);
pp7_4_all_night_amp_pos = pp7_4_all_night.amp_pos_total(pp7_4_all_night.amp_pos_total<lim);


Qb5_all_night_amp_neg = Qb5_all_night.amp_neg_total(abs(Qb5_all_night.amp_neg_total)<lim);
Qb10_all_night_amp_neg = Qb10_all_night.amp_neg_total(abs(Qb10_all_night.amp_neg_total)<lim);
pp7_24_all_night_amp_neg = pp7_24_all_night.amp_neg_total(abs(pp7_24_all_night.amp_neg_total)<lim);
pp7_4_all_night_amp_neg = pp7_4_all_night.amp_neg_total(abs(pp7_4_all_night.amp_neg_total)<lim);


lambda = 1:3;
norm = 50:0.5:600;

sub_func(pp7_4_all_night_amp_pos,lambda,norm,'PP7 4 positive');
sub_func(pp7_4_all_night_amp_neg,lambda,norm,'PP7 4 negative');

sub_func(Qb5_all_night_amp_pos,lambda,norm,'Qb 5 positive');
sub_func(Qb5_all_night_amp_neg,lambda,norm,'Qb 5 negative');

sub_func(Qb10_all_night_amp_pos,lambda,norm,'Qb 10 positive');
sub_func(Qb10_all_night_amp_neg,lambda,norm,'Qb 10 negative');

sub_func(pp7_24_all_night_amp_pos,lambda,norm,'PP7 24 positive');
sub_func(pp7_24_all_night_amp_neg,lambda,norm,'PP7 24 negative');


end

 
function [] = sub_func(data,lambda,norm,title_str)
[nn,edges] = histcounts(data,'BinWidth',50,'Normalization','probability');
if edges(1) < 0
    edges = fliplr(edges);
    nn = fliplr(nn);
end
nn_full = nn;
nn = transpose(smooth(nn,10));
hist_x_axis = abs(mean([edges(1:end-1);edges(2:end)]));


mse = zeros(length(lambda),length(norm));
mse_x_axis = zeros(1,length(norm));
for j = 1:length(lambda)
    for i=1:length(norm)
        normalizer = norm(i);
        mse_x_axis(1,i) = normalizer;
        xx = 0:max(hist_x_axis)/normalizer;
        yy = poisspdf(xx,lambda(j));
        yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(nn)));
        mse(j,i) = sum((nn(~isnan(yq))./max(nn)-yq(~isnan(yq))).^2);
    end
end

nn_full = transpose(smooth(nn_full,3));


minMatrix = min(mse(:));
[row,col] = find(mse==minMatrix);
lambda_single = lambda(1,row);
normalizer = mse_x_axis(1,col);
xx = 0:max(hist_x_axis)/normalizer;
yy = poisspdf(xx,lambda_single);
yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(nn)));

% plot fit with lowest MSE
figure;
scatter(hist_x_axis,nn_full);  xlim([0,2000]); ylim([0,0.15]);
hold on;
plot(hist_x_axis,yq*max(nn)./max(yq),'r*-');
box on;
title([title_str,' normalization factor = ',num2str(normalizer),' lambda = ',num2str(lambda_single),' MSE= ',num2str(minMatrix)]);
xlabel('Intensity [A.U]');
ylabel('Probability');

% plot all fits
figure;
scatter(hist_x_axis,nn_full,'DisplayName',[title_str,' experimental data']);  
xlim([0,2000]); ylim([0,0.15]);
box on;
hold on;
minMatrix = min(mse(1:3,:),[],2);
graphic_set = {['--r'];['g-.'];['k:']};
for i = 1:length(minMatrix)
    [row,col] = find(mse==minMatrix(i));
    lambda_single = lambda(1,row);
    normalizer = mse_x_axis(1,col);
    xx = 0:max(hist_x_axis)/normalizer;
    yy = poisspdf(xx,lambda_single);
    yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(nn)));
    legend_txt = ['lambda=',num2str(i),' K0=',num2str(normalizer),' MSE=',num2str(minMatrix(i))];
    plot(hist_x_axis,yq*max(nn)./max(yq),graphic_set{i},'LineWidth',1.5,'DisplayName',legend_txt);
end
xlabel('Intensity [A.U]');
ylabel('Probability');
legend('show');


end
