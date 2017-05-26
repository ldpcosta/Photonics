function power = getPowerFromFilename(filex)
    [startIndex,endIndex] = regexp(filex,'_P.+\.');
    power = str2num(filex(startIndex+2:endIndex-1));
end

