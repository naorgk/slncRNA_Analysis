%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes as input granule intensity files, and background
% intensity files. For each signal (each granule), significant entry and
% exit events are found and their parameters measured.
% The output is a .mat file containing amplitude data and temporal data
% about burst events.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter settings
close all;
clear;

window_span = 13; % Moving average filter window span (Required to reduce signal noise)
match_probability = 12; % Search for a segment with probability of 1/(2^match_probability)
timeVectorTotal = linspace(0,360/6,360+1); % For a full movie there should be 361 time points

[signal_file,path] = uigetfile('','Select signal data file'); % open the signal file 
signal_data = struct2cell(load([path,signal_file]));
signal_data = signal_data{1};
signal_data = signal_data(2:end);

[background_file,path] = uigetfile('','Select background data file'); % open the background file
normalization_data = struct2cell(load([path,background_file]));
normalization_data = normalization_data{1};
normalization_data = normalization_data(2:end);

% control single inetensity data plotting

plot_max = 10;

%% Variable definition
plot_index = 0;
time_stat_pos = []; % Duration of positive events
time_stat_neg = []; % Duration of negative events
time_stat_var = []; % Duration of quiescent events
time_between_pos = []; % Time between succesive positive events
time_between_neg = []; % Time between succesive negative events
amp_pos_total = []; % Positive amplitudes
amp_neg_total = []; % Negative amplitudes
amp_var_total = []; % Quiescent amplitudes

meanIntensity = [];
numberOfEvents = 0;
numberOfSignalWithEvents = 0;


%% Data collection process 

for i = 1:length(signal_data)
    granule_signal = signal_data{i};
    background_signal = normalization_data{i};
    % Limits of track lengths. Discard signals with less than 180 time
    % points (30 minutes of tracking)
    if length(granule_signal)<180 || length(granule_signal)>360
        continue;
    end
    
    % Time vector for current signals. Could be less than 360 time points
    timeVector = timeVectorTotal(1:length(granule_signal));

     % Fit background signal to a 3rd degree polynomial to be used 
     % as normalization. This removes the high frequency noise in the
     % background.
    [f,gof,output] = fit(timeVector',background_signal,'poly3');
    background_signal = feval(f,timeVector);

    % Moving average smoothing of granule signal to remove high frequency
    % noise
    granule_signal = filloutliers(granule_signal,'nearest','mean');
    granule_signal = smoothdata(granule_signal,'movmean',window_span);  
    
    % Generate filtered, normalized granule signal.
    data = ((granule_signal./background_signal-ones(size(granule_signal))))*background_signal(1);

    % Get start positions and lengths of segments classified as burst
    % events
    [startPosValid,lengthsValid,threshold,threshold_neg] = get_position(data,match_probability);
    numberOfEvents = numberOfEvents+(length(startPosValid));

    if ~isempty(startPosValid)
        numberOfSignalWithEvents = numberOfSignalWithEvents+1;
        meanIntensity = [meanIntensity,mean(data)];
    end
    
    % these three variables are used to determine temporal squence of
    % events. i.e. find out how the segments are organized in time
    posSegmentX = [];
    negSegmentX = [];
    varSegmentX = [];

    % find activity events
    for tt = 1:length(startPosValid)
        x = linspace(startPosValid(tt),startPosValid(tt)+lengthsValid(tt),...
            lengthsValid(tt)+1);
        % Due to the linear fit, the segment length might exceed the length
        % of the signal. This code segment checks for and fixes the issue
        if x(end)>length(data)
            x = x(x<length(data)); 
            p{tt} = polyfit(x',data(startPosValid(tt):length(data)-1),1);
        else
            p{tt} = polyfit(x',data(startPosValid(tt):startPosValid(tt)+lengthsValid(tt)),1); 
        end
        % linear fit of the signal, at the appropriate position, to get an
        % increasing, decreasing, or quiescent segment.
        y = polyval(p{tt},x);

        if y(end)-y(1) >= 0
            % Increasing segment - Entry event
            amp_pos_total = [amp_pos_total; y(end)-y(1)];
            time_stat_pos = [time_stat_pos; lengthsValid(tt)*(timeVector(2)-timeVector(1))];
            posSegmentX = [posSegmentX; x(1)];
        elseif length(x)>=threshold_neg % check against negative threshold
            % Decreasign segment - Exit event
            amp_neg_total = [amp_neg_total; y(end)-y(1)];
            time_stat_neg = [time_stat_neg; lengthsValid(tt)*(timeVector(2)-timeVector(1))];
            negSegmentX = [negSegmentX; x(1)];
        end
    end
    % Classify unidentified segments as quiescent
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
        time_stat_var = [time_stat_var; length(xVar)*(timeVector(2)-timeVector(1))];
        varSegmentX = [varSegmentX; xVar(1)];
    end
    time_between_pos = [time_between_pos; (timeVector(2)-timeVector(1))*diff(posSegmentX)];
    time_between_neg = [time_between_neg; (timeVector(2)-timeVector(1))*diff(negSegmentX)];

    %%%%%%%%%%%%%%%% segment classification over %%%%%%%%%%%%%%%%
    if ~exist('p','var')
        p = [];
    end
    if ~exist('pVar','var')
        pVar = [];
    end
    
    % Plot signals overlayed with identified segments
    if plot_index < plot_max
        plot_tracks_custom(timeVector,data,plot_index,startPosValid,lengthsValid,p,pVar,[])
    end
    plot_index = plot_index+1;
end


%% Plots

% Histogram of all amplitudes
figure; 
histogram(amp_var_total(abs(amp_var_total)<3000),'FaceColor','b','BinWidth',100);
hold on;
histogram(amp_pos_total(amp_pos_total<3000),'FaceColor','g','BinWidth',100); 
histogram(amp_neg_total(abs(amp_neg_total)<3000),'FaceColor','r','BinWidth',100);
ylabel('Counts'); xlim([-3000,3000])
title('Segment amplitudes','FontSize',12);
xlabel('Amplitudes [A.U]','FontSize',12);

%Boxplot of duration between events
figure;
plot_var = [time_between_pos;time_between_neg];
group_var = zeros(length(plot_var),1);
group_var(1:length(time_between_pos))=2;
group_var(length(time_between_pos)+1:end)=3;
boxplot(plot_var,group_var);
xticklabels({'Positive','Negative'});
ylabel('Time [min]');
title('Time between burst events');

%% Save data to .mat file
uisave({'amp_neg_total','amp_pos_total','amp_var_total','meanIntensity',...
    'time_between_neg','time_between_pos',...
    'time_stat_neg','time_stat_pos','time_stat_var'});