function [y,ind] = isclose2(Y, x)
% function [y ind] = isclose2(Y, x)
% 
% Function that determines the value and index of Y closest to x

szY = size(Y);
szx = size(x);

if (szY(1) == 1 && szY(2) == 1)
    disp('The first input should be a vector or matrix');
    y = Y;
    ind = 1;
else
    if ((szY(1) == 1 && szY(1) == 1 && szY(2) == szx(2)) || ...
        (szY(2) == 1 && szY(2) == 1 && szY(1) == szx(1)))
        % Y and x are two vectors with the same dimension
%         [val,ind] = min(abs(x-n));
%         n1 = x(ind(1));
%         y = Y(x == n1);
        error('The second input must be a number');
    elseif ((szx(1) == 1 && szx(2) == szY(2)) || ...
            (szx(2) == 1 && szx(1) == szY(2)))
        [val,ind] = min(abs(Y-x));
        y = Y(ind(1));
%         y = Y(:,x == n1);
    elseif ((szx(2) == 1 && szx(1) == szY(1)) || ...
            (szx(1) == 1 && szx(2) == szY(1)))
        [val,ind] = min(abs(Y-x));
        y = Y(ind(1));
%         y = Y(x == n1,:);
    else
        error('Please, review the size of inputs');   
    end
end


    
