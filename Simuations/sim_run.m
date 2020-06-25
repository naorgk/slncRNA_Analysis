clearvars;
close all;

numTracks = 1000;
t = linspace(0,60,360);
Y_f = zeros(numTracks,360);


for i = 1:numTracks 
    coeff1 = randi(1e2);
%     coeff1 = 0;
    coeff2 = randi(1e3);
    a = coeff1;
    b = coeff2;
    y = a*t+b;
    n = 40*randn(1,360);
    yn = y+n;
    if i<=5
        figure;
        subplot(2,1,1)
        plot(t,y);
        ylabel('Intensity [A.U]');
        title('Base Signal');
        subplot(2,1,2);
        plot(t,yn);
        title('With added noise');
        xlabel('Time [sec]');
        ylabel('Intensity [A.U]');
    end
    Y_f(i,:) = yn;
end 

ampVec = [80,120];

for i = 1:numTracks
    y = zeros(1,360);
    n = 50*randn(1,360); % 570 is the mean of standard deviations from all experiments
    signVec = datasample([-1,0,0,0,0,1],40);
    
    for j = 1:9:360
        y(j) = datasample(ampVec,1)*signVec(ceil(j/9));
        y(j+1:j+8) = y(j);
        
        for k = 1:7
           y(j+k) = y(j+k)+k*randi(10);
        end
       
    end
    y = y+1000;
    yn = y+n;
    Y_f(i,:) = yn;  
    if i<=4
        figure;
        subplot(2,1,1)
        plot(t,y);
        ylabel('Intensity [A.U]');
        title('Base Signal');
        subplot(2,1,2);
        plot(t,yn);
        title('With added noise');
        xlabel('Time [sec]');
        ylabel('Intensity [A.U]');

    end
end

% 
norm = ones(size(Y_f))*randi(100);
norm = norm+40*randn(numTracks,360);
save('sim_run_slope_300.mat','Y_f');
save('sim_run_norm_slope_300.mat','norm')