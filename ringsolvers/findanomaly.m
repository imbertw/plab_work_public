function [anomaly] = findanomaly(x,y)
Y = y; %initial data to manipulate and delete bad points
Y1 = y; %final data to manipulate, delete bad points, and fit to good points
X = x;
newflips = [];
flips = [];
flag = 1;

%Below is the code that gets rid of flips where there are multiple bad
%points in a row. After the first pass, the points in the "valley" or
%"plateau" are not identified as flips.
%Therefore you have to delete the initial flips, and then search again for
%flips, until there are none left.
%---\        /-----    ----\      /----     -----\    /------   --------------
%    \      /       =>      \    /       =>       \  /        =>
%     \____/                 \__/                  \/

while flag == 1
    Y(newflips) = 0;
    idx_to_keep = find(Y~=0);
    Y = Y(idx_to_keep);
    X = X(idx_to_keep);
    dy = diff(Y);
    dx = diff(X);
    dy_dx = dy./dx;
    d2y_dx2 = diff(dy_dx)./dx(2:end);
    % Find anomalies at the end of the vector
    endupflip = [];
    enddownflip = [];
    if dy_dx(end) > 0 & dy_dx(end-1) < 0
        endupflip = length(dy_dx)-1;
    elseif dy_dx(end) < 0 & dy_dx(end-1) > 0
        enddownflip = length(dy_dx)-1;
    end
    newupflips = find(dy_dx(2:end-1) > 0 & dy_dx(3:end) < 0 & d2y_dx2(1:end-1) > 0 & d2y_dx2(2:end) < 0);
    newdownflips = find(dy_dx(2:end-1)<0 & dy_dx(3:end) > 0 & d2y_dx2(1:end-1) < 0 & d2y_dx2(2:end) > 0); 
    newupflips = [newupflips endupflip];
    newdownflips = [newdownflips enddownflip];
    newflips = [newupflips newdownflips];
    newflips = newflips + 2;
    newflips = unique(newflips);
    d = ismember(y,Y([newflips])); %Find where new flips are indexed in the original data
    g = find(d);
    if(~isempty(newflips))    %
        flips = [flips g]; %Add new flips to global ledger
    else
        % Find outliers which have continuous derivatives and second
        % derivatives
        hampel_Y = hampel(Y);
        std_Y = sqrt((hampel_Y-Y).^2);
        newupflips = find(dy_dx(1:end-1) > 0 & dy_dx(2:end) < 0 & std_Y(2:end-1) > 0.1*mean(std_Y) );
        newdownflips = find(dy_dx(1:end-1) < 0 & dy_dx(2:end) > 0 & std_Y(2:end-1) > 0.1*mean(std_Y) );
        newflips = [newupflips newdownflips];
        newflips = unique(newflips);
        if(~isempty(newflips))
            %% Adding more flips to global ledger
            newflips = newflips + 1;
            d = ismember(y,Y([newflips])); % This bit of code is why we have both y and Y. y houses the original data, and we want to see where the indices of the flips are in the original data, not the modified data.
            g = find(d); % find(array) gives you the indices where the array is not zero
            flips = [flips g];
        else
            %% One more check--let's take a look at the third derivative to see any discontinuities
            d3y_dx3 = diff(d2y_dx2)./dx(3:end);
            newflips = find(dy_dx(2:end-2).*dy_dx(3:end-1) < 0 & d2y_dx2(2:end-1).*d2y_dx2(3:end) < 0 & d3y_dx3(1:end-1).*d3y_dx3(2:end) < 0 );
            newflips = newflips + 3;
            if(~isempty(newflips))
                d = ismember(y,Y([newflips])); % This bit of code is why we have both y and Y. y houses the original data, and we want to see where the indices of the flips are in the original data, not the modified data.
                g = find(d); % find(array) gives you the indices where the array is not zero
                flips = [flips g];
            else
                flag = 0;
            end
        end
    end
end
%% Post-selection of anomalies
% Y1(flips) = 0;
% idx_to_keep = find(Y1~=0);
anomaly = flips;
%yfit = interp1(x(idx_to_keep),y(idx_to_keep),x, 'linear','extrap');
%plot(x,yfit);
%anomaly = find((abs(yfit-y) > 0.1*std(yfit))==1); % Find where the points deviate majorly from the fit
end