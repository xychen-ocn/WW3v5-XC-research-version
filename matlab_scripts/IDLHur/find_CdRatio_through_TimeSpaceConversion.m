% This is an updated version of the previous function.
% Now, the function will interpolate results at the exact depth instead of
% finding results at the nearest depth. 
%
function [data_t, data_Ratio, TT,YT,coordN]=find_CdRatio_through_TimeSpaceConversion(data,DW, XX,YY,dpt_q, dpt_vec,timenum_data,stmx,stmy,timenum_stm)
% specifically build for horizontalling moving stationary hurricane.
% input: data, XX, YY (data and its corresponding space coordinate info.)
%        stmx: storm center in x-dir (in units of meter!)
%        timenum_data, timenum_stm (julian datenum of the data and storm track
%        info.) 
%        BE AWARE: stmx and XX, YY should be in the same units! METER!
% output: data_t: time evolution of data: (dimension: Ydist, time/sudo-space)
%         TT, YT : space-time coordinate
%         XX_sudo: space-sudospace coordinate;
%
% Xuanyu Chen; Jul 3, 2018 (generalized from previous scripts)
%
%  [data_t, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate(indata.(varn), ...
%            coord.XX,coord.YY,depID, ww3_timenum,stmx*1000,stmy*1000, stmtime);
%

% find when the storm "hit" the selected xlocation: torig_id
% that time is set to the "origin" for the selected xlocation.
%
% New explaination: This code basically extract the data time series at a given
% location (deps, (a.k.a. x-dist)), then find the corresponding
% psudo-xdistance and use that to replace the time dimension.
% June 15, 2019: added deep water input to find the ratio more accurately.
%  It is required that the DW dataset has at least the same number of time
%  entry of the data.
% 


    xvec=XX(1,:);
    yvec=YY(:,1);
    mps2kph=3.6;
    data_t=zeros(size(data,1), size(data,3));
% 1. determine the depID according to the depth at request:

[tmp,depID]=min(abs(dpt_q-dpt_vec));
%tmp
%depID
if tmp~=0
    % interpolation needed:
    % find the corresponding cross shore distance at which depth=dpt_q
   % depID
    nonanIDs=find(isnan(dpt_vec)==0);
    x_at_dptq=interp1(dpt_vec(nonanIDs),xvec(nonanIDs),dpt_q);
    %% data extraction:
    xvec_q=x_at_dptq.*ones(size(yvec));
    
    for it=1:size(data,3)
        data_t(:,it)=interp2(XX,YY,data(:,:,it),xvec_q,yvec);
        % make sure the deep water dataset has the same dimension as the
        % shallow water counterparts:
        DWdata_intp=interp2(DW.XX, DW.YY, DW.data(:,:,it),XX,YY);
        data_tDW(:,it)=interp2(XX,YY,DWdata_intp, xvec_q, yvec);
    end
    
    % it is likely that one of these two are problematic.
    data_Ratio=data_t./data_tDW;      
    
else
    idx=depID;
    x_at_dptq=xvec(idx);
    %% data extraction:
    data_t= squeeze(data(:,idx,:));
    
    for it=1:size(data,3);
        DWdata_intp(:,:,it)=interp2(DW.XX, DW.YY, DW.data(:,:,it),XX,YY);
    end
    data_tDW=squeeze(DWdata_intp(:,idx,:));
    
    disp('exact depth is found')
    % it is likely that one of these two are problematic.
    data_Ratio=data_t./data_tDW;
       
end

%% build coordinate:
    tspd=(stmx(2:end)-stmx(1:end-1))./((timenum_stm(2:end)-timenum_stm(1:end-1))*86400);
    tspd_ave=mean(tspd,'omitnan');   % in m/s;
    disp(['mean_tspd=', num2str(tspd_ave,'%6.2f') 'm/s']);
    
    %% find out when the storm center reach the selected depth.
    [xdif, torg_id]=min(abs(x_at_dptq-stmx));     % xvec and stmx need to be in the same unit, in meter!
    xdif=sign(x_at_dptq-stmx(torg_id))*xdif;
    
    if xdif <0 
        tspd_inst=tspd(torg_id); 
    else
        tspd_inst=tspd(torg_id-1);
    end 
    
    % add minor correction: if xdif ~= 0
    % find the correct tspd at the time of closesest to the depth contour.
    dt = xdif./abs(tspd_inst);    % in second
    dt_in_day=dt./86400;
    datestr(timenum_stm(torg_id))
    timnum_stm_at_xloc=timenum_stm(torg_id);%+dt_in_day;   %this assumes that timenum_stm and stmx has the same output freq. and aligned.
    
    % compute the relative time to the time when storm hit xloc:
    time2=(timenum_data-timnum_stm_at_xloc)*24;    %in hour
    
    % build time-space mesh:
    [TT,YT]=meshgrid(time2, yvec);
    XX_sudo=TT.*abs(tspd_ave)*mps2kph;      % in km.
    %(stmx(torg_id)+dt_in_dy*86400*tspd_ave)/1000
    coordN.XX=XX_sudo ; %in km
    coordN.YY=(YT - stmy(torg_id))/1000;     %in km
    
%     figure
%     subplot(2,2,1);
%     pcolor(coordN.XX, coordN.YY, data_t);shading interp;
%     hold on;
%     plot([0 0],[-600 600]);
%     colorbar;
%     
%     subplot(2,2,2);
%     pcolor(coordN.XX, coordN.YY, data_Ratio);shading interp;
%     colorbar;
%         hold on;
%     plot([0 0],[-600 600]);
% 
%     pause; close gcf;

    
return