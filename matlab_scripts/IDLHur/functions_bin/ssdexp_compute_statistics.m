function [data_stat]=ssdexp_compute_statistics(indata);
data_median=median(indata,'omitnan');           % median
data_mean=mean(indata,'omitnan');
data_wisEdge=prctile(indata,[2.5 100-2.5],1);   % get the data value corresponding to 2.5% and 97.5% in the distribution.
data_boxEdge=prctile(indata,[25 75],1);   % get the data value corresponding to 2.5% and 97.5% in the distribution.

data_IQR = iqr(indata,1);                       % interquartile range (75%-25%)
rel_95Var=(data_wisEdge(2,:)-data_wisEdge(1,:))./data_mean.*100;
rel_IQRVar=data_IQR./data_mean.*100;

rel_ave95var= mean(abs((data_wisEdge - repmat(data_mean,2,1))) ./repmat(data_mean,2,1),1)*100;
rel_ave50var= mean(abs((data_boxEdge - repmat(data_mean,2,1))) ./repmat(data_mean,2,1),1)*100;


data_stat.median=data_median;
data_stat.mean=data_mean;
data_stat.IQR = data_IQR;
data_stat.wisEdge=data_wisEdge;
data_stat.boxEdge=data_boxEdge;
data_stat.rel95var=rel_95Var;
data_stat.rel50var=rel_IQRVar;

data_stat.averel95var=rel_ave95var;
data_stat.averel50var=rel_ave50var;

return
