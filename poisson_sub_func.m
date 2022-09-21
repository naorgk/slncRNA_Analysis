%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principle of operation:
% The fluorescence intensity data is first binned into bins of equal width to
% construct a counting process. 
% For each lambda value, a theoretical Poisson distribution is generated to
% comapre to. Then, for each K0 value, the fluorescence data and the bin
% edges are normalized by K0 to get an x-axis on which a Poisson
% distribution is defined (i.e., [0-10]). The mean squared error (MSE) between
% the theoretical Poisson distribution and the experimental, normalized
% counts is computed. 
% The output is the lambda and K0 values which minimize the MSE to the
% theoretical Poisson distribution.

% Input:
% - Fluorescence data
% - Array of lambda values to iterate on
% - Array of K0 values to iterate on
% - binWidth - width of the bins to divide the data into
% - String to be used as a title on the plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [hist_x_axis,counts_full,yq_out,counts] = poisson_sub_func(data,lambda,K0,binWidth,title_str)
% bin the fluorescence data into bins of equal width
[counts,edges] = histcounts(data,'BinWidth',binWidth,'Normalization','probability');
if edges(1) < 0
    edges = fliplr(edges);
    counts = fliplr(counts);
end
counts_full = counts;
counts = transpose(smooth(counts,10)); % filter to reduce noise
hist_x_axis = abs(mean([edges(1:end-1);edges(2:end)])); % get bin edges


mse = zeros(length(lambda),length(K0));
mse_x_axis = zeros(1,length(K0));
for j = 1:length(lambda)
    for i=1:length(K0)
        normalizer = K0(i);
        mse_x_axis(1,i) = normalizer;
        xx = 0:max(hist_x_axis)/normalizer+1; % +1 to avoid an error on some data sets. should not affect results
        yy = poisspdf(xx,lambda(j)); % generate theoretical distribution
        % interpolate to get same number of points for MSE calculation
        yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(counts)));
        mse(j,i) = sum((counts(~isnan(yq))./max(counts)-yq(~isnan(yq))).^2);
    end
end

% As lambda and K0 arrays can be very large,this section downsamples them 
% for plotting purposes only.
K0_short = K0(1:20:end);
lambda_short = lambda;
mse_short = mse(:,1:20:end);
counts_full = transpose(smooth(counts_full,3));

% Heatmap of MSE values
figure;
heatmap(K0_short,lambda_short,mse_short,'CellLabelColor','none'); colorbar();
colormap(parula); caxis([0,2*round(mean(mse_short(:)))]);
xlabel('Model single molecule intensity');
ylabel('Poisson rate');
title('Poisson fitting parameters');


% Find lambda, K0 values that minimize MSE
minMatrix = min(mse(:));
[row,col] = find(mse==minMatrix);
lambda_single = lambda(1,row);
normalizer = mse_x_axis(1,col);
xx = 0:max(hist_x_axis)/normalizer+1;
yy = poisspdf(xx,lambda_single);
yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(counts)));
% Plot best fit
figure;
scatter(hist_x_axis,counts_full);  
hold on;
plot(hist_x_axis,yq*max(counts)./max(yq),'r*-');
box on;
title([title_str,' normalization factor = ',num2str(normalizer),' lambda = ',num2str(lambda_single),' MSE= ',num2str(minMatrix)]);
xlabel('Intensity [A.U]');
ylabel('Probability');

% Return best fit 
yq_out = yq;


% Plot QQ-plot of best fit for visual inspection
figure;
qqplot(counts_full,yq*max(counts)./max(yq));
title('QQ plot of sample data vs. theoretical distribution');
xlabel('Quantiles of sample data');
ylabel('Quantiles of theoretical distribution');

