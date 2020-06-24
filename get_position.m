function [startPosValid,lengthsValid,threshold,threshold_neg] = get_position(data,probability,window_span)

    % The gradient function calculates a numerical gradient of the data
    % using the central difference method. Using this, single opposite
    % points are allowed in the sequence.

    dataDiff = diff(data); 
    
    posChance = length(find(dataDiff>0))/length(dataDiff);
    negChance = 1-posChance;
        
    threshold = ceil(-probability/log2(posChance));
    threshold_neg = ceil(-probability/log2(negChance));
    
    if threshold < window_span
        threshold = window_span;
    end
    if threshold_neg < window_span
        threshold_neg = window_span;
    end
       
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

 