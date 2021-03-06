function [time_stat_pos,time_stat_neg,time_stat_var,time_between_pos,...
    rate_pos,rate_neg,amp_pos,amp_neg,amp_var,postPos,postNeg,postVar,...
    numberOfEvents,numberOfSignalWithEvents,eventsPerSignal] = ...
    analysis(spot,cell,timeVecTotal,window_span,match_probability,plot_min,plot_index,plot_max)

time_stat_pos = -Inf;
time_stat_neg = -Inf;
time_stat_var = -Inf;
time_between_pos = -Inf;
rate_pos = -Inf;
rate_neg = -Inf;
amp_pos = -Inf;
amp_neg = -Inf;
amp_var = -Inf;
postPos = [];
postNeg = [];
postVar = [];

numberOfEvents = 0;
numberOfSignalWithEvents = 0;
eventsPerSignal = [];
for i = 1:size(spot,1)
    measure = spot(i,:);
    drop = cell(i,:);
    timeVec = timeVecTotal(1:length(measure));

    [f,gof,output] = fit(timeVec',drop','poly3');
    drop = feval(f,timeVec);
    measure = filloutliers(measure,'nearest','mean');
    measure = smoothdata(measure,'movmean',window_span);
    data = ((measure./drop'-ones(size(measure))))*drop(1);

    [startPosValid,lengthsValid,threshold,threshold_neg] = get_position(data',match_probability,window_span);
    numberOfEvents = numberOfEvents+(length(startPosValid));


    if ~isempty(startPosValid)
        numberOfSignalWithEvents = numberOfSignalWithEvents+1;
        eventsPerSignal = [eventsPerSignal,length(startPosValid)];
    end
    time_between_pos= [time_between_pos; (timeVec(2)-timeVec(1))*diff(startPosValid)'];
    % these three variables are used to determine temporal squence of
    % events. i.e. find out how the segments are organized in time
    posSegmentX = [];
    negSegmentX = [];
    varSegmentX = [];

        % find activity events
        for tt = 1:length(startPosValid)
            x = linspace(startPosValid(tt),startPosValid(tt)+lengthsValid(tt),...
                lengthsValid(tt)+1);
            if x(end)>length(data)
                x = x(1:end-1);
                p{tt} = polyfit(x,data(startPosValid(tt):startPosValid(tt)+lengthsValid(tt)-1),1);
            else
                p{tt} = polyfit(x,data(startPosValid(tt):startPosValid(tt)+lengthsValid(tt)),1); 
            end
            y = polyval(p{tt},x);
            if y(end)-y(1) >= 0
                amp_pos = [amp_pos; y(end)-y(1)];
                time_stat_pos = [time_stat_pos; lengthsValid(tt)*(timeVec(2)-timeVec(1))];
                rate_pos = [rate_pos; (y(end)-y(1))/(lengthsValid(tt)*(timeVec(2)-timeVec(1)))];
                posSegmentX = [posSegmentX; x(1)];
            elseif length(x)>=threshold_neg
                amp_neg = [amp_neg; y(end)-y(1)];
                time_stat_neg = [time_stat_neg; lengthsValid(tt)*(timeVec(2)-timeVec(1))];
                rate_neg = [rate_neg; (y(end)-y(1))/(lengthsValid(tt)*(timeVec(2)-timeVec(1)))];
                negSegmentX = [negSegmentX; x(1)];
            end
        end
        
        % find quiescent segments
        for tt = 1:length(startPosValid)-1
            xVar = linspace(startPosValid(tt)+lengthsValid(tt)+1,startPosValid(tt+1)-1,...
                startPosValid(tt+1)-startPosValid(tt)-lengthsValid(tt)-1);
            if isempty(xVar)
                continue;
            end
            if xVar(end)>length(data)
                xVar = xVar(1:end-1);
                pVar{tt} = polyfit(xVar,data(startPosValid(tt)+lengthsValid(tt)+1:startPosValid(tt+1))-1,1);
            else
                pVar{tt} = polyfit(xVar,data(startPosValid(tt)+lengthsValid(tt)+1:startPosValid(tt+1)-1),1);
            end
            yVar = polyval(pVar{tt},xVar);
            amp_var = [amp_var; (yVar(end)-yVar(1))];
            time_stat_var = [time_stat_var; length(xVar)*(timeVec(2)-timeVec(1))];
            varSegmentX = [varSegmentX; xVar(1)];
        end
        
        positive_class = 1*ones(size(posSegmentX));
        negative_class = -1*ones(size(negSegmentX));
        var_class = 1*zeros(size(varSegmentX));
        
        classes = [positive_class;negative_class;var_class];
        [xxx,iii] = sort([posSegmentX;negSegmentX;varSegmentX]);
        classes = classes(iii);
        
        for k=1:length(classes)-1
            switch classes(k)
                case -1
                    postNeg = [postNeg; classes(k+1)];
                case 0
                    postVar = [postVar; classes(k+1)];
                case 1
                    postPos = [postPos; classes(k+1)];              
            end
        end
        
        
        %%%%%%%%%%%%%%%% segment classification over %%%%%%%%%%%%%%%%
        if ~exist('p','var')
            p = [];
        end
        if ~exist('pVar','var')
            pVar = [];
        end   
        if plot_index >= plot_min && plot_index < plot_max
            plot_tracks_custom(timeVec,data,plot_index,startPosValid,lengthsValid,p,pVar,[])
        end

        plot_index = plot_index+1;
end