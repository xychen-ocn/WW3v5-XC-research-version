function [Cd_stdv2,Cd_wmean,wspd_cen, Cd_binned, U10_binned]=compute_stdv_from_Cd_scatters_at_each_windspd_bin(Cd,U10,U10_binwidth, U10_binrange)
% purpose: this function is used to compute standard deviation of the Drag
% Coefficient scatters.
% 1. sort Cd scatter by wind speed bins.
[Cd_binned, U10_binned,wspd_cen, idx_binned]=sort_data_by_windspeed(Cd, U10,U10_binwidth, U10_binrange);
disp(['number of bins:' num2str(length(Cd_binned))]);
Cd_stdv=[]; Cd_wmean=[];
for i=1:length(Cd_binned)
    % 2. weighted average Cd in each bin:
    [counts, centers]=hist(Cd_binned{i},5);   %sort Cd_binned based on its value into 5 bins.
%      figure(2)
%      bar(centers,counts)
     totnum=sum(counts);
    % compute weighted average:
    Cd_wmean(i)=sum(counts./totnum.*centers);
    %Cd_mean(i)=mean(Cd_binned{i},'omitnan');

    % 3. compute standard deviation:
   % Cd_stdv(i)=std(Cd_binned{i},1);         % 1 -> averaged over N instead of N-1;
    % used the weighted averaged mean to compute standard deviation:
    sqdif=(Cd_binned{i}-Cd_wmean(i)).^2;
    Cd_stdv2(i)=sqrt(sum(sqdif)/totnum);
end


% compare 2 different definitions of Cd: same here.
%     figure(3); clf;
% %     plot(wspd_cen,Cd_wmean,'-r');
% %     hold on
% %     plot(wspd_cen,Cd_mean,'-k');
%     plot(Cd_stdv,'ob');
%     hold on
%     plot(Cd_stdv2,'*r');


end


%%%  nested function  %%%
% function [Cd_binned, U10_binned, center_wspd]=sort_data_by_windspeed(Cd, U10,U10_binwidth)
% U10min=min(U10);
% U10max=max(U10);
% wspd=[U10min:U10_binwidth:U10max];
% center_wspd=0.5*(wspd(1:end-1)+wspd(2:end));
% disp(['number of wspd bins:' num2str(length(wspd))]);
% 
% clear Cd_binned U10_binned
% for ibin=1:length(wspd)-1
%     mask= NaN(size(Cd));
%     mask(U10>=wspd(ibin))=1;
%     mask(U10>=wspd(ibin+1))=NaN;
%     id = find(isnan(mask)==0);
%     Cd_binned{ibin}=Cd(id);
%     U10_binned{ibin}=U10(id);
% end
% end
%%%  nested function  %%%

    