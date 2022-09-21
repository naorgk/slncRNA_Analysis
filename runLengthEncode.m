function [startPos,lengths, values] = runLengthEncode(data)
% This function finds segments of uninteruppted increase or decrease
startPos = find(diff([data(1)-1, data]));
lengths = diff([startPos, numel(data)+1]);
values = data(startPos);