function [startPos,lengths, values] = runLengthEncode(data)
startPos = find(diff([data(1)-1, data]));
lengths = diff([startPos, numel(data)+1]);
values = data(startPos);