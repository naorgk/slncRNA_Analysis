%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads all .tiff files containing microscopy movies, and .mat
% files containing (x,y,frame) data about detected fluorescent granules.
% The data from each movie is loaded, converted to an array of images, and
% the fluorescence data (intensity of the granule, and
% intensity of its surrounding environment) from each granule is collected
% and saved for further analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

% Path to folder containing tif and mat files
% file_path = './8x/green/';
file_path = './8x/red/';

% load all .tiff and all .mat files
tiff_files = dir([file_path,'*.tif']);
mat_files = dir([file_path,'*.mat']);

% Length of the subframe to
sub_frame_length = 7;

% Initialize variables
intensity = cell(1,1);
background = cell(1,1);

% Iterate over files, convert to array of images
for f = 1:size(tiff_files,1)
    movArray = generate_movArray(file_path,tiff_files(f).name);
    % load (x,y) coordinates of spots detected in imageJ
    xydata = load([file_path,mat_files(f).name]);
    xdata = xydata.x;
    ydata = xydata.y;
    framedata = xydata.frame;
    
    % Start iterating through spots and gather intensity
    for t = 1:length(xdata)
        spot_frames = cell2mat(framedata(t)); 
        spot_frames = spot_frames(2:end); % remove frame 0 to avoid error
        spot_x = cell2mat(xdata(t)); spot_x = spot_x(2:end);
        spot_y = cell2mat(ydata(t)); spot_y = spot_y(2:end);
        if length(spot_x)>=24 % continue if spot exists for >24 frames
            % iterate over frames and get intensity data
            spot_intensity = zeros(length(spot_x),1);
            background_intensity = zeros(length(spot_x),1);
            
            for n = 1:length(spot_x)
               frame = movArray(spot_frames(n)).data; % complete current frame
               % Extract sub frame around the spot coordiantes
               sub_frame_bounds = [max(spot_y(n)-sub_frame_length,1),...
                                   min(spot_y(n)+sub_frame_length,length(frame)),...
                                   max(spot_x(n)-sub_frame_length,1),...
                                   min(spot_x(n)+sub_frame_length,length(frame))];
               sub_frame_bounds = round(sub_frame_bounds);
               sub_frame = frame(sub_frame_bounds(1):sub_frame_bounds(2),sub_frame_bounds(3):sub_frame_bounds(4));
               % sub frame is ready from this point
               
               % Separate to spot, and 2 levels of background
               th = multithresh(sub_frame,4);
               spot_frame = zeros(size(sub_frame));
               cell_frame = zeros(size(sub_frame));
               back_frame = zeros(size(sub_frame));
               spot_frame(sub_frame>=th(3)) = 1;
               cell_frame(sub_frame>=th(1) & sub_frame<th(3)) = 1;
               back_frame(sub_frame<th(1)) = 1;
               
               % frame partition to spot,background is done
               
               spot_intensity(n,1) = mean(sub_frame(spot_frame==1));
               background_intensity(n,1) = mean(sub_frame(cell_frame==1));
               
            end
                      
           intensity = [intensity,{spot_intensity}];
           background = [background,{background_intensity}];
        end
    end
end
% Intensity signal gathering is done. 

save('PCP_8x_spot.mat','intensity');
save('PCP_8x_background.mat','background');