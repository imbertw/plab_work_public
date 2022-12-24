function [yfit] = fixanomaly(x,y,anomaly)
nans = find(isnan(y));
anomaly = [anomaly nans];
yfit = interp1(x(setdiff(1:end,anomaly)),y(setdiff(1:end,anomaly)), x, 'linear','extrap');
%ytestfit = y(setdiff(1:end,anomaly));
%ytestfit = hampel(ytestfit,5);
%yfit = hampel(ytestfit,5);
end