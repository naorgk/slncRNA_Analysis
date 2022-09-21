function [] = plot_tracks_custom(timeVec,data,plot_index,startPosValid,lengthsValid,p,pVar,f)

if isempty(f)
    figure;
end
plot(timeVec,data);
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
        plot(x*(timeVec(2)-timeVec(1)),y,'g','LineWidth',1.2);
    else
        plot(x*(timeVec(2)-timeVec(1)),y,'r','LineWidth',1.2);
    end
end
% Plot blue linear line overlayed on quiescent segments


for tt = 1:length(startPosValid)-1
    xVar = linspace(startPosValid(tt)+lengthsValid(tt)+1,startPosValid(tt+1)-1,...
        startPosValid(tt+1)-startPosValid(tt)-lengthsValid(tt)-1);
    if isempty(xVar)
        continue;
    end
    % pay attention: if pVar contains edges than this has to be pVar{tt+2}
    yVar = polyval(pVar{tt},xVar);
    plot(xVar*(timeVec(2)-timeVec(1)),yVar,'b','LineWidth',1.2);   
end

return

