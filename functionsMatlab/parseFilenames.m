function [filenames, types] = parseFilenames(extension, directory)

%FALTA ACABAR OS TYPES
    direc = dir(strcat(directory,'\*',extension));
    filenames = {};
    types = {};
    
    for i = 1:size(direc,1)
        filename = direc(i).name;
        if ~(strcmp(filename,'.')) && ~(strcmp(filename,'..'))
            filenames = [filenames {filename}];
        end
        
        data = strsplit(filename,'_');
                
    end
   
    
end