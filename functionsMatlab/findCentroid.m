function trueCenter = findCentroid(estimCenter, bandwidth, dataArray)
    
    count = 0;
    norm = 0;
    for ind = 1:bandwidth/2
       count = count + dataArray(estimCenter-ind)*ind + dataArray(estimCenter+ind)*ind;
       norm = norm + ind;
    end
    
    trueCenter = count/norm;
    
end