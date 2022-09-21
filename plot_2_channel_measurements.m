%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads .tiff files (containing microscopy movies), and .mat
% files (containing (x,y,frame) data), from two measured channels (e.g.,
% GFP and mCherry). The detected granules in both channels are matched
% together by position, and intensity signals from each granule are plotted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

%% Parameter definitions 
% File path variables
channel_1_file_path = './8x/green/';
channel_2_file_path = './8x/red/';

% plotting variables
plot_max = 10; % Maximal number of plots to generate

% segment classification variables
match_probability = 12; % Search for a segment with probability of 1/(2^match_probability)
window_span = 13; % Moving average filter window span (Required to reduce signal noise)
%% Iterate over files and convert to array of images

% load all .tiff and all .mat files
channel_1_tiff_files = dir([channel_1_file_path,'*_original.tif']);
channel_1_mat_files = dir([channel_1_file_path,'*.mat']);

channel_2_tiff_files = dir([channel_2_file_path,'*_original.tif']);
channel_2_mat_files = dir([channel_2_file_path,'*.mat']);

sub_frame_length = 7;
% Initialize variables
intensity = cell(2,1);
background = cell(2,1);
% Variable to hold correlation measurements between the two channles
correlation_coefficients = []; 
parameters = struct('start_frame',1,'show_video',0);

for f = 1:size(channel_1_tiff_files,1) % iterate over files
    channel_1_movArray = generate_movArray(channel_1_file_path,channel_1_tiff_files(f).name);
    channel_2_movArray = generate_movArray(channel_2_file_path,channel_2_tiff_files(f).name);
    channel_1_xydata = load([channel_1_file_path,channel_1_mat_files(f).name]);
    channel_1_xdata = channel_1_xydata.x;
    channel_1_ydata = channel_1_xydata.y;
    channel_1_framedata = channel_1_xydata.frame;
    
    channel_2_xydata = load([channel_2_file_path,channel_2_mat_files(f).name]);
    channel_2_xdata = channel_2_xydata.x;
    channel_2_ydata = channel_2_xydata.y;
    channel_2_framedata = channel_2_xydata.frame;
    
    % iterate over green spots and find closest red spots to each
    for k = 1:length(channel_1_xdata)
        spot_x = cell2mat(channel_1_xdata(k)); spot_x = spot_x(2:end);
        channel_1_x_position = spot_x(1,1);
        spot_y = cell2mat(channel_1_ydata(k)); spot_y = spot_y(2:end);
        channel_1_y_position = spot_y(1,1);
        spot_frames = cell2mat(channel_1_framedata(k)); 
        spot_frames = spot_frames(2:end);
        
        % calculate distance between the current green spot and all red
        % spots
        distance_vector = zeros(length(channel_2_xdata),2);
        for kk = 1:length(channel_2_xdata)
            channel_2_spot_x = cell2mat(channel_2_xdata(kk)); channel_2_spot_x = channel_2_spot_x(2:end);
            channel_2_x_position = channel_2_spot_x(1,1);
            channel_2_spot_y = cell2mat(channel_2_ydata(kk)); channel_2_spot_y = channel_2_spot_y(2:end);
            channel_2_y_position = channel_2_spot_y(1,1);
            distance_vector(kk,1) = sqrt((channel_2_x_position-channel_1_x_position)^2+(channel_2_y_position-channel_1_y_position)^2);
            distance_vector(kk,2) = length(channel_2_spot_x);
        end
        % Discard red spots visible for less than 24 frames.
        distance_vector = distance_vector(distance_vector(:,2)>=24);
        [min_distance,index_of_min] = min(distance_vector(:,1));
        channel_2_spot_x = cell2mat(channel_2_xdata(index_of_min)); channel_2_spot_x = channel_2_spot_x(2:end);
        channel_2_spot_y = cell2mat(channel_2_ydata(index_of_min)); channel_2_spot_y = channel_2_spot_y(2:end);
        channel_2_spot_frames = cell2mat(channel_2_framedata(index_of_min)); 
        channel_2_spot_frames = channel_2_spot_frames(2:end);
        
        if (length(spot_x)>=180) & (length(channel_2_spot_x)>=180)
            % iterate over frames and get intensity data
            channel_1_spot_intensity = zeros(length(spot_x),1);
            channel_1_cell_intensity = zeros(length(spot_x),1);
            channel_2_spot_intensity = zeros(length(channel_2_spot_x),1);
            channel_2_cell_intensity = zeros(length(channel_2_spot_x),1);
            % GREEN data
            for n = 1:length(spot_x)
               frame = channel_1_movArray(spot_frames(n)).data; % complete current frame
               % Extract sub frame around the spot coordiantes
               sub_frame_bounds = [max(spot_y(n)-sub_frame_length,1),...
                                   min(spot_y(n)+sub_frame_length,length(frame)),...
                                   max(spot_x(n)-sub_frame_length,1),...
                                   min(spot_x(n)+sub_frame_length,length(frame))];
               sub_frame_bounds = round(sub_frame_bounds);
               sub_frame = frame(sub_frame_bounds(1):sub_frame_bounds(2),sub_frame_bounds(3):sub_frame_bounds(4));
               % subFrame is ready from this point
               
               th = multithresh(sub_frame,4);
               spot_frame = zeros(size(sub_frame));
               cell_frame = zeros(size(sub_frame));
               back_frame = zeros(size(sub_frame));
               test_frame = zeros(size(sub_frame));
               spot_frame(sub_frame>=th(3)) = 1;
               cell_frame(sub_frame>=th(1) & sub_frame<th(3)) = 1;
               back_frame(sub_frame<th(1)) = 1;
               
               % frame partition to spot,background is done
               
               channel_1_spot_intensity(n,1) = mean(sub_frame(spot_frame==1));
               channel_1_cell_intensity(n,1) = mean(sub_frame(cell_frame==1));
               
            end
            
            % RED data
            for n = 1:length(channel_2_spot_x)
               frame = channel_2_movArray(channel_2_spot_frames(n)).data; % complete current frame
               % Extract sub frame around the spot coordiantes
               sub_frame_bounds = [max(channel_2_spot_y(n)-sub_frame_length,1),...
                                   min(channel_2_spot_y(n)+sub_frame_length,length(frame)),...
                                   max(channel_2_spot_x(n)-sub_frame_length,1),...
                                   min(channel_2_spot_x(n)+sub_frame_length,length(frame))];
               sub_frame_bounds = round(sub_frame_bounds);
               sub_frame = frame(sub_frame_bounds(1):sub_frame_bounds(2),sub_frame_bounds(3):sub_frame_bounds(4));
               % subFrame is ready from this point
               
               th = multithresh(sub_frame,4);
               spot_frame = zeros(size(sub_frame));
               cell_frame = zeros(size(sub_frame));
               back_frame = zeros(size(sub_frame));
               test_frame = zeros(size(sub_frame));
               spot_frame(sub_frame>=th(3)) = 1;
               cell_frame(sub_frame>=th(1) & sub_frame<th(3)) = 1;
               back_frame(sub_frame<th(1)) = 1;
               
               % frame partition to spot,background is done
              
               channel_2_spot_intensity(n,1) = mean(sub_frame(spot_frame==1));
               channel_2_cell_intensity(n,1) = mean(sub_frame(cell_frame==1));
               
           end
                      
        intensity = [intensity,[{channel_1_spot_intensity};{channel_2_spot_intensity}]];
        background = [background,[{channel_1_cell_intensity};{channel_2_cell_intensity}]];
        end

    end


end
% Intensity signal gathering complete 
%% Identify segments of exit, entry and quiescence,plot signals.
% Calculate filtered, normalized intensity signal for plotting and segment
% classification

plot_index = 0;
timeVecTotal = linspace(0,360/6,360+1); % For a full movie there should be 361 time points
intensity = intensity(:,2:end);
background = background(:,2:end);
for i = 1:size(intensity,2)
    % Green signal
    granule_signal = intensity{1,i};
    background_signal = background{1,i};
    
    % Time vector for current signals. Could be less than 360 time points
    channel_1_timeVec = timeVecTotal(1:length(granule_signal));
    % Fit background signal to a 3rd degree polynomial to be used 
    % as normalization. This removes the high frequency noise in the
    % background.
    [f,gof,output] = fit(channel_1_timeVec',background_signal,'poly3');
    background_signal = feval(f,channel_1_timeVec);
    % Moving average smoothing of granule signal to remove high frequency
    % noise
    granule_signal = filloutliers(granule_signal,'nearest','mean');
    granule_signal = smoothdata(granule_signal,'movmean',window_span);  
    % Generate filtered, normalized granule signal.
    channel_1_data = ((granule_signal./background_signal-ones(size(granule_signal))))*background_signal(1);
    
    % Red signal
    granule_signal = intensity{2,i};
    background_signal = background{2,i};
    
    channel_2_timeVec = timeVecTotal(1:length(granule_signal));
    % Fit background signal to a 3rd degree polynomial to be used 
    % as normalization. This removes the high frequency noise in the
    % background.
    [f,gof,output] = fit(channel_2_timeVec',background_signal,'poly3');
    background_signal = feval(f,channel_2_timeVec);
    % Moving average smoothing of granule signal to remove high frequency
    % noise
    granule_signal = filloutliers(granule_signal,'nearest','mean');
    granule_signal = smoothdata(granule_signal,'movmean',window_span);  
    % Generate filtered, normalized granule signal.
    channel_2_data = ((granule_signal./background_signal-ones(size(granule_signal))))*background_signal(1);
    
    % Calculate correlation coefficients between the signals
    % Since the signals might not be the same length, the correlation is
    % calculated up to the length of the shorter signal
    min_length = min(length(channel_1_data),length(channel_2_data));
    coeff_mat = corrcoef(smooth(channel_1_data(1:min_length,:),60),...
        smooth(channel_2_data(1:min_length,:),60));
    
    correlation_coefficients = [correlation_coefficients; coeff_mat(1,2)];
    
    
    %%%% START SEGMENT CLASSIFICATION
    
    % Get start positions and lengths of segments classified as burst
    % events
    [channel_1_startPosValid,channel_1_lengthsValid,~,~] = get_position(channel_1_data,match_probability);
    [channel_2_startPosValid,channel_2_lengthsValid,~,~] = get_position(channel_2_data,match_probability);

    % find activity events
    for tt = 1:length(channel_1_startPosValid)
        x = linspace(channel_1_startPosValid(tt),channel_1_startPosValid(tt)+channel_1_lengthsValid(tt),...
            channel_1_lengthsValid(tt)+1);
        if x(end)>length(channel_1_data)
           x = x(x<length(channel_1_data));
           channel_1_p{tt} = polyfit(x',channel_1_data(channel_1_startPosValid(tt):length(channel_1_data)-1),1);
        else
           channel_1_p{tt} = polyfit(x',channel_1_data(channel_1_startPosValid(tt):channel_1_startPosValid(tt)+channel_1_lengthsValid(tt)),1); 
        end
        % linear fit of the signal, at the appropriate position, to get an
        % increasing, decreasing, or quiescent segment.
        y = polyval(channel_1_p{tt},x);
    end
    % Classify unidentified segments as quiescent    
    channel_1_pVar = [];
    for tt = 1:length(channel_1_startPosValid)-1
        xVar = linspace(channel_1_startPosValid(tt)+channel_1_lengthsValid(tt)+1,channel_1_startPosValid(tt+1)-1,...
            channel_1_startPosValid(tt+1)-channel_1_startPosValid(tt)-channel_1_lengthsValid(tt)-1);
        if isempty(xVar)
            continue;
        end
        if xVar(end)>=length(channel_1_data)
            xVar = xVar(1:end-1);
            channel_1_pVar{tt} = polyfit(xVar',channel_1_data(channel_1_startPosValid(tt)+channel_1_lengthsValid(tt)+1:channel_1_startPosValid(tt+1))-1,1);
        else
            channel_1_pVar{tt} = polyfit(xVar',channel_1_data(channel_1_startPosValid(tt)+channel_1_lengthsValid(tt)+1:channel_1_startPosValid(tt+1)-1),1);
        end
        yVar = polyval(channel_1_pVar{tt},xVar);
    end
    
    %%%%%%%%%%%%%%%
    % find activity events -RED
    for tt = 1:length(channel_2_startPosValid)
        x = linspace(channel_2_startPosValid(tt),channel_2_startPosValid(tt)+channel_2_lengthsValid(tt),...
            channel_2_lengthsValid(tt)+1);
        if x(end)>length(channel_2_data)
           x = x(x<length(channel_2_data));
           channel_2_p{tt} = polyfit(x',channel_2_data(channel_2_startPosValid(tt):length(channel_2_data)-1),1);
        else
           channel_2_p{tt} = polyfit(x',channel_2_data(channel_2_startPosValid(tt):channel_2_startPosValid(tt)+channel_2_lengthsValid(tt)),1); 
        end
        % linear fit of the signal, at the appropriate position, to get an
        % increasing, decreasing, or quiescent segment.
        y = polyval(channel_2_p{tt},x);
    end
    % Classify unidentified segments as quiescent    
    channel_2_pVar = [];
    for tt = 1:length(channel_2_startPosValid)-1
        xVar = linspace(channel_2_startPosValid(tt)+channel_2_lengthsValid(tt)+1,channel_2_startPosValid(tt+1)-1,...
            channel_2_startPosValid(tt+1)-channel_2_startPosValid(tt)-channel_2_lengthsValid(tt)-1);
        if isempty(xVar)
            continue;
        end
        if xVar(end)>=length(channel_2_data)
           xVar = xVar(1:end-1);
           channel_2_pVar{tt} = polyfit(xVar',channel_2_data(channel_2_startPosValid(tt)+channel_2_lengthsValid(tt)+1:channel_2_startPosValid(tt+1))-1,1);
        else
           channel_2_pVar{tt} = polyfit(xVar',channel_2_data(channel_2_startPosValid(tt)+channel_2_lengthsValid(tt)+1:channel_2_startPosValid(tt+1)-1),1);
        end
        % linear fit of the signal, at the appropriate position, to get an
        % increasing, decreasing, or quiescent segment.
        yVar = polyval(channel_2_pVar{tt},xVar);
    end


    % Segment classification over 
    if ~exist('p','var')
        p = [];
    end
    if ~exist('pVar','var')
        pVar = [];
    end

    if ~exist('selective_plot','var')
        
        if plot_index < plot_max
           figure;
           subplot(2,1,1);
           plot_tracks_custom_subplot_version(channel_1_timeVec,channel_1_data,plot_index,channel_1_startPosValid,channel_1_lengthsValid,channel_1_p,channel_1_pVar,gca())
           subplot(2,1,2);
           plot_tracks_custom_subplot_version(channel_2_timeVec,channel_2_data,plot_index,channel_2_startPosValid,channel_2_lengthsValid,channel_2_p,channel_2_pVar,gca())
        end
    end
    plot_index = plot_index+1;
    
end

