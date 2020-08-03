function [ydata_binned_matrix, xdata_binned, xbin_center,idx_binned]=sort_ydata_by_xdata(ydata, xdata,xdata_binwidth, xdata_binrange)
% this function can acutually be applied to not only windspeed..think xdata
% as another quantity used to connect to the to be binned data.
% output: [ydata_binned_matrix, xdata_binned, xbin_center,idx_binned]
if nargin == 3
    % automatically find xdata range:
    xdatamin=min(xdata);
    xdatamax=max(xdata);
else
    xdatamin=xdata_binrange(1);
    xdatamax=xdata_binrange(2);
end
xbin=[xdatamin:xdata_binwidth:xdatamax];
xbin_center=0.5*(xbin(1:end-1)+xbin(2:end));
disp(['number of wspd bins:' num2str(length(xbin))]);

clear ydata_binned xdata_binned idx_binned
for ibin=1:length(xbin)-1
    mask= NaN(size(ydata));
    mask(xdata>=xbin(ibin))=1;
    mask(xdata>=xbin(ibin+1))=NaN;
    id = find(isnan(mask)==0);
    ydata_binned{ibin}=ydata(id);
    xdata_binned{ibin}=xdata(id);
    idx_binned{ibin}=id;
end

% combine cell into a matrix: (Updated on Jul 4, 2018)
%%% __ combine bins into 1 matrix: __ %%%
num_of_bins=length(ydata_binned);
%%%%%% searching for the maximum dimension:
maxnl=length(ydata_binned{1});
for i=1:num_of_bins-1;
    maxnl=max(maxnl, length(ydata_binned{i+1}) );
end
ydata_binned_matrix=NaN(maxnl,num_of_bins);
for i=1:num_of_bins
    nl=length(ydata_binned{i});
    ydata_binned_matrix(1:nl,i)=ydata_binned{i};
end

end