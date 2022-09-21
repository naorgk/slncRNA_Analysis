%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function plot_tracks_custom_subplot_version recieves an intensity signal,
% and identified segments in that signal and plots them inside figure axes
% provided as input. 
% Used to plot measurements from 2 channels in the
% plot_2_channel_measurements script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_tracks_custom_subplot_version(timeVec,data,plot_index,startPosValid,lengthsValid,p,pVar,ax)


plot(ax,timeVec,data);
title(num2str(plot_index));
xlabel('Time [Min]');
ylabel('Intensity [A.U]');
xlim([0,60]); 
hold on;
% Plot red and green linear lines on significant segments
for tt = 1:length(startPosValid)
    x = linspace(startPosValid(tt),startPosValid(tt)+lengthsValid(tt),...
        lengthsValid(tt)+1);
    y = polyval(p{tt},x);
    if y(end)>=y(1)
        plot(ax,x*(timeVec(2)-timeVec(1)),y,'g','LineWidth',1.2);
    else
        plot(ax,x*(timeVec(2)-timeVec(1)),y,'r','LineWidth',1.2);
    end
end
% Plot blue linear line overlayed on quiescent segments

for tt = 1:length(startPosValid)-1
    xVar = linspace(startPosValid(tt)+lengthsValid(tt)+1,startPosValid(tt+1)-1,...
        startPosValid(tt+1)-startPosValid(tt)-lengthsValid(tt)-1);
    if isempty(xVar)
        continue;
    end
    yVar = polyval(pVar{tt},xVar);
    plot(ax,xVar*(timeVec(2)-timeVec(1)),yVar,'b','LineWidth',1.2);   
end

return

