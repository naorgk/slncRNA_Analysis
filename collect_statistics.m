% This script extracts traces from  data files, detects burst
% segments and collects relevant statistics

close all;
clear;
window_span = 13;
match_probability = 10; % search for a segment with probability of 1/2^match_probability
timeVecTotal = linspace(0,360/6,360+1);
options = statset('MaxIter',1000);



% spot = struct2cell(load('Trace files./all_night_induction_4pp7_data.mat'));
% cell = struct2cell(load('Trace files./all_night_induction_4pp7_norm.mat'));

spot = struct2cell(load('Trace files./all_night_induction_5Qb_data.mat'));
cell = struct2cell(load('Trace files./all_night_induction_5Qb_norm.mat'));

% spot = struct2cell(load('Trace files./all_night_induction_10Qb_data.mat'));
% cell = struct2cell(load('Trace files./all_night_induction_10Qb_norm.mat'));

% spot = struct2cell(load('Trace files./all_night_induction_24pp7_data.mat'));
% cell = struct2cell(load('Trace files./all_night_induction_24pp7_norm.mat'));


spot = spot{1}; spot = spot(2:end);
cell = cell{1}; cell = cell(2:end);

time_stat_pos = [];
time_stat_neg = [];
time_stat_var = [];
time_between_pos = [];
time_between_neg = [];
time_between_positions = [];
rate_pos = [];
rate_neg = [];
amp_pos_total = [];
amp_neg_total = [];
amp_var_total = [];

mean_raw_measure = [];
meanmRNA = [];
numberOfEvents = 0;
numberOfSignalWithEvents = 0;
segment_length = [];

% collect segments immedietly following segments
postPos = [];
postNeg = [];
postVar = [];

background_std = [];
background_mean = [];

% control single inetensity data plotting
% use these variables to plot batches of traces. 
% usage example: set plot_min = 10; plot_max = 20;
plot_min = 0;
plot_index = 0;
plot_max = 0;

% selectiv_plot variable is used to plot segments of your choice. 
% for example: 
% selective_plot = [ 27,101]; % for the Qb-5x data

% selective_plot = [];

drop_denoise_error = [];
for i = 1:length(spot)
    measure = spot{i};
    drop = cell{i};
    % limits of track lengths. 240 is default. 
    if length(measure)<240 || length(measure)>360
        continue;
    end
    
    raw_measure = measure;
    raw_drop = drop;
    timeVec = timeVecTotal(1:length(measure));
    [f,gof,output] = fit(timeVec',raw_drop,'poly3');
    drop = feval(f,timeVec);
    
    % error statistical data
    drop_denoise_error = [drop_denoise_error; gof.rmse./mean(raw_drop)];
    background_std = [background_std; std(drop)];
    
    measure = filloutliers(measure,'nearest','mean');
    measure = smoothdata(measure,'movmean',window_span);  
    
    % normalize intensity signals
    data = ((measure./drop-ones(size(measure))))*drop(1);

    [startPosValid,...
        lengthsValid,...
        threshold,...
        threshold_neg] = get_position(data,match_probability,window_span);
    
    numberOfEvents = numberOfEvents+(length(startPosValid));
    segment_length = [segment_length,lengthsValid];

    if ~isempty(startPosValid)
        numberOfSignalWithEvents = numberOfSignalWithEvents+1;
        meanmRNA = [meanmRNA,mean(data)];
        mean_raw_measure = [mean_raw_measure,mean(raw_measure)];
        background_mean = [background_mean; mean(raw_drop)];
    end
    
    time_between_positions = [time_between_positions; (timeVec(2)-timeVec(1))*diff(startPosValid)'];
    % these three variables are used to determine temporal squence of
    % events. i.e. find out how the segments are organized in time
    posSegmentX = [];
    negSegmentX = [];
    varSegmentX = [];
    % these following variables are used to convert amplitude signal to a
    % current signal
    posTime = [];
    negTime = [];
    varTime = [];
        % find activity events
        for tt = 1:length(startPosValid)
            x = linspace(startPosValid(tt),startPosValid(tt)+lengthsValid(tt),...
                lengthsValid(tt)+1);
            if x(end)>length(data)
                x = x(1:end-1);
                p{tt} = polyfit(x',data(startPosValid(tt):startPosValid(tt)+lengthsValid(tt)-1),1);
            else
                p{tt} = polyfit(x',data(startPosValid(tt):startPosValid(tt)+lengthsValid(tt)),1); 
            end
            y = polyval(p{tt},x);

            if y(end)-y(1) >= 0
                amp_pos_total = [amp_pos_total; y(end)-y(1)];
                time_stat_pos = [time_stat_pos; lengthsValid(tt)*(timeVec(2)-timeVec(1))];
                posSegmentX = [posSegmentX; x(1)];
                posTime = [posTime; lengthsValid(tt)*(timeVec(2)-timeVec(1))];
            elseif length(x)>=threshold_neg
                amp_neg_total = [amp_neg_total; y(end)-y(1)];
                time_stat_neg = [time_stat_neg; lengthsValid(tt)*(timeVec(2)-timeVec(1))];
                negSegmentX = [negSegmentX; x(1)];
                negTime = [negTime; lengthsValid(tt)*(timeVec(2)-timeVec(1))];
            end
        end
        
        pVar = [];
        for tt = 1:length(startPosValid)-1
            xVar = linspace(startPosValid(tt)+lengthsValid(tt)+1,startPosValid(tt+1)-1,...
                startPosValid(tt+1)-startPosValid(tt)-lengthsValid(tt)-1);
            if isempty(xVar)
                continue;
            end
            if xVar(end)>=length(data)
                xVar = xVar(1:end-1);
                pVar{tt} = polyfit(xVar',data(startPosValid(tt)+lengthsValid(tt)+1:startPosValid(tt+1))-1,1);
            else
                pVar{tt} = polyfit(xVar',data(startPosValid(tt)+lengthsValid(tt)+1:startPosValid(tt+1)-1),1);
            end
            yVar = polyval(pVar{tt},xVar);
            amp_var_total = [amp_var_total; (yVar(end)-yVar(1))];
            time_stat_var = [time_stat_var; length(xVar)*(timeVec(2)-timeVec(1))];
            varSegmentX = [varSegmentX; xVar(1)];
            varTime = [varTime; length(xVar)*(timeVec(2)-timeVec(1))];
        end

        time_between_pos = [time_between_pos; (timeVec(2)-timeVec(1))*diff(posSegmentX)];
        time_between_neg = [time_between_neg; (timeVec(2)-timeVec(1))*diff(negSegmentX)];
        
        positive_class = 1*ones(size(posSegmentX));
        negative_class = -1*ones(size(negSegmentX));
        var_class = 1*zeros(size(varSegmentX));
        
        classes = [positive_class;negative_class;var_class];
        time_for_current = [posTime;negTime;varTime];
        start_positions_full = [posSegmentX;negSegmentX;varSegmentX];
        [~,iii] = sort(start_positions_full);
        classes = classes(iii);
        time_for_current = time_for_current(iii);

        start_positions_significant = [posSegmentX;negSegmentX];
        [start_positions_significant,jj] = sort(start_positions_significant);
        
        for k=1:length(classes)-1
            switch classes(k)
                case -1
                    postNeg = [postNeg; classes(k+1)];
                case 0
                    if time_for_current(k) >= 2.5
                        postVar = [postVar; classes(k+1)];
                    end
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
        if ~exist('selective_plot','var')
            % No multiple track selection, plot tracks independantly
            if plot_index >= plot_min && plot_index < plot_max
                plot_tracks_custom(timeVec,data,plot_index,startPosValid,lengthsValid,p,pVar,[])
            end
        else
            if ~ishandle(1)
               f = figure;
               hold on;
            end
            if ~isempty(find(selective_plot==plot_index))
               
                if length(timeVec)>320
                    plot_tracks_custom(timeVec,data,plot_index,startPosValid,lengthsValid,p,pVar,f)
                end
            end  
        end
        plot_index = plot_index+1;
end

figure;
histogram(time_stat_pos,100);
ylabel('Counts');
xlabel('Time [min]');
title('Duration of positive segments');
xlim([0,10]);

figure;
histogram(time_stat_neg,100);
ylabel('Counts');
xlabel('Time [min]');
title('Duration of negative segments');
xlim([0,10]);

figure;
histogram(time_stat_var,300);
xlabel('Time [min]');
title('Duration of quiescent segments');
xlim([0,25]);


figure;
histogram(time_between_positions(time_between_positions<=45),'BinWidth',1);
xlabel('Time [min]');
ylabel('Counts');
title('Time between significant events');
xlim([0,45]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot histogram of all amplitudes with 3 Gaussians overlaid corresponding
% to Negative, Zero, Positive

figure; 
histogram(amp_var_total(abs(amp_var_total)<3000),'FaceColor','b','BinWidth',100);
hold on;
histogram(amp_pos_total(amp_pos_total<3000),'FaceColor','g','BinWidth',100); 
histogram(amp_neg_total(abs(amp_neg_total)<3000),'FaceColor','r','BinWidth',100);
ylabel('Counts'); xlim([-3000,3000])
title('Segment amplitudes','FontSize',12);
xlabel('Amplitudes [A.U]','FontSize',12);


rate_pos = amp_pos_total./time_stat_pos;
rate_neg = amp_neg_total./time_stat_neg;
rate_var = amp_var_total./time_stat_var;


figure; 
histogram(rate_pos(amp_pos_total<3000),'FaceColor','g','BinWidth',50); hold on;
histogram(rate_var(abs(amp_var_total)<3000),'FaceColor','b','BinWidth',50);
histogram(rate_neg(abs(amp_neg_total)<3000),'FaceColor','r','BinWidth',50);
xlim([-1000,1000])
ylabel('Counts','FontSize',12);
xlabel('Slope [A.U/Min]','FontSize',12);
title('Segment slopes');


categories = [{'Negative'},{'Quiescent'},{'Positive'}];
figure;
histogram(postNeg);
xticks([-1,0,1]);
xticklabels(categories);
title('Following negative segment');

figure;
histogram(postVar);
xticks([-1,0,1]);
xticklabels(categories);
title('Following quiescent segment');

figure;
histogram(postPos);
xticks([-1,0,1]);
xticklabels(categories);
title('Following positive segment');



% save statistical data to file
uisave({'amp_pos_total','amp_neg_total','amp_var_total',...
    'time_stat_pos','time_stat_neg','time_stat_var',...
    'meanmRNA','background_mean','mean_raw_measure','time_between_pos',...
    'time_between_neg','postVar'});