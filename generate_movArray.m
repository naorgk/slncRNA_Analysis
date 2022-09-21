% Function generate_movArray reads a tiff file containing a sequence of
% images into a matlab array prior to reading.

function [movArray] = generate_movArray(filePath,fileName)


movArray = struct('data',[]);
fname = fullfile(filePath,fileName);
info = imfinfo(fname);
num_images = numel(info);
for k = 1:num_images
    movArray(k).data = imread(fname, k);
end 


end
