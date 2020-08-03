function [IQR,stdv,bias]=variability_analysis(data, given_mean)
% input: data (a 1 dimensional vector) data does not include NaN;
%        given_mean (a mean to be compared the data against)
N=length(data);
IQR=iqr(data);
stdv=sqrt(sum((data-given_mean).^2)./N);   % standard deviation from the given mean value.
bias=sum(data - given_mean)/N;

end