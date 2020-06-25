% This function generates simulated constant signals with overlaid noise
function [Y_f,norm] = sim_const(numTracks,exp_const,noise_std,t)
Y_f = zeros(numTracks,360);

for i = 1:numTracks 
    coeff1 = 0;
    coeff2 = randi([2e3,1e4]);
    y = coeff1*t+coeff2;
    ye = y.*exp(-exp_const*t);
    n = noise_std*randn(1,360);
    yn = ye+n;
    if i<=5
        figure;
        subplot(2,1,1)
        plot(t,y,'LineWidth',2); hold on; plot(t,ye,'LineWidth',2);
        ylim([ye(end)-10,y(end)+10]);
        ylabel('Intensity [A.U]');
        legend('Base signal','With additional exponential component','Location','Southwest');
        subplot(2,1,2);
        plot(t,yn,'g','LineWidth',1.5);
        legend('With additional noise');
        xlabel('Time [min]');
        ylabel('Intensity [A.U]');
    end
    Y_f(i,:) = yn;
end 

norm = ones(size(Y_f))*randi([1000,1500]);
norm = norm.*exp(-exp_const*t);
norm = norm+noise_std*randn(numTracks,360);

end

