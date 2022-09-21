function [startPosValid,lengthsValid,threshold,threshold_neg] = get_position(data,probability)
    
    % Function get_position finds increasing or decreasing signal segments 
    % which are statistically unlikly to occur in a random signal.

    dataDiff = diff(data); 
    
    posChance = length(find(dataDiff>0))/length(dataDiff);
    negChance = 1-posChance;
        
    threshold = ceil(-probability/log2(posChance));
    threshold_neg = ceil(-probability/log2(negChance));
    
    % In some cases, the algorithm skips over segments due to single points
    % with a derivative around zero. This section looks for segments of 5
    % points in which the middle point has a derivative in the opposite
    % direction of the segment (i.e. negative derivative in an otherwise
    % increasing segment), and corrects it.
    signVec = sign(dataDiff);
    signVec_corrected = signVec;
    for i = 3:length(dataDiff)-2
        if (signVec(i-2) == signVec(i-1)) && (signVec(i-1) == signVec(i+1)) ...
                && (signVec(i+1) == signVec(i+2)) && (signVec(i-1)~=signVec(i))
            signVec_corrected(i) = signVec(i-1);
        end
    end
    
    signVec_corrected(signVec_corrected==-1) = 0;
    [startPos,lengths, ~] = runLengthEncode(signVec_corrected');
    points = lengths>threshold;
    lengthsValid = lengths(points);
    startPosValid = startPos(points);
end

 