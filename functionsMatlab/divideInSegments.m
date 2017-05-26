
function segArray = divideInSegments(dataArray, segNumber)
    segArray = [];
    for tt = 1:segNumber
        segArray(tt,:) = dataArray{:, tt}.';
    end
end