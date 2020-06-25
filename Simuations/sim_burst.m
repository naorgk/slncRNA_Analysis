% This function generates simulated signals with intermittent bursts and
% overlaid noise
function [Y_f,norm] = sim_burst(numTracks,exp_const,noise_std,t)
Y_f = zeros(numTracks,360);

ampVec = [80,120];

for i = 1:numTracks
    y = zeros(1,360);
    n = noise_std*randn(1,360);
    signVec = datasample([-1,0,0,0,0,1],40);
    
    for j = 1:9:360
        y(j) = datasample(ampVec,1)*signVec(ceil(j/9));
        y(j+1:j+8) = y(j);
        
        for k = 1:7
           y(j+k) = y(j+k)+k*randi(10);
        end
    end
    
    y = y+randi([2e3,1e4]);
    ye = y.*exp(-exp_const*t);
    yn = ye+n;
    Y_f(i,:) = yn;  
    if i<=5
        figure;
        subplot(2,1,1)
        plot(t,y,'LineWidth',2); hold on; plot(t,ye,'LineWidth',2);
        ylabel('Intensity [A.U]');
        legend('Base signal','With additional exponential component','Location','Southwest');
        subplot(2,1,2);
        plot(t,yn,'g','LineWidth',1.5);
        legend('With additional noise');
        xlabel('Time [min]');
        ylabel('Intensity [A.U]');
    end
end

norm = ones(size(Y_f))*randi([1000,1500]);
norm = norm.*exp(-exp_const*t);
norm = norm+noise_std*randn(numTracks,360);

end

