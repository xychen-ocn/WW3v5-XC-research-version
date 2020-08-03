function [Cd_binned, U10_binned, center_wspd,idx_binned]=sort_data_by_windspeed(Cd, U10,U10_binwidth, U10_binrange)
if nargin == 3
    % automatically find U10 range:
    U10min=min(U10);
    U10max=max(U10);
else
    U10min=U10_binrange(1);
    U10max=U10_binrange(2);
end
wspd=[U10min:U10_binwidth:U10max];
center_wspd=0.5*(wspd(1:end-1)+wspd(2:end));
disp(['number of wspd bins:' num2str(length(wspd))]);

clear Cd_binned U10_binned idx_binned
for ibin=1:length(wspd)-1
    mask= NaN(size(Cd));
    mask(U10>=wspd(ibin))=1;
    mask(U10>=wspd(ibin+1))=NaN;
    id = find(isnan(mask)==0);
    %pause
    
    if ~isempty(id)
        Cd_binned{ibin}=Cd(id);
        U10_binned{ibin}=U10(id);
        idx_binned{ibin}=id;
    else
        Cd_binned{ibin}=NaN;
        U10_binned{ibin}=NaN;
        idx_binned{ibin}=NaN;
    end
end
end