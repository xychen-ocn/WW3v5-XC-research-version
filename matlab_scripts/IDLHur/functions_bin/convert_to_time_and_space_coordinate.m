function [data_t, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate(data,XX,YY,idx,timenum_data,stmx,stmy,timenum_stm)
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

% find when the storm "hit" the selected xlocation: torig_id
% that time is set to the "origin" for the selected xlocation.
xvec=XX(1,:);
yvec=YY(:,1);
mps2kph=3.6;

tspd=(stmx(2:end)-stmx(1:end-1))./((timenum_stm(2:end)-timenum_stm(1:end-1))*86400);
tspd_ave=mean(tspd,'omitnan');   % in m/s;
%disp(['mean_tspd=', num2str(tspd_ave,'%6.2f') 'm/s']);

[xdif, torg_id]=min(abs(xvec(idx)-stmx));     % xvec and stmx need to be in the same unit, in meter!

% add minor correction: if xdif ~= 0 
dt = xdif/abs(tspd_ave);    % in second
dt_in_dy=dt./86400;
timnum_stm_at_xloc=timenum_stm(torg_id)+dt_in_dy;   %this assumes that timenum_stm and stmx has the same output freq. and aligned.

% compute the relative time to the time when storm hit xloc:
time2=(timenum_data-timnum_stm_at_xloc)*24;    %in hour

% build time-space mesh:
[TT,YT]=meshgrid(time2, yvec);
XX_sudo=TT.*abs(tspd_ave)*mps2kph;      % in km. 
%(stmx(torg_id)+dt_in_dy*86400*tspd_ave)/1000
coordN.XX=XX_sudo ; %in km
coordN.YY=(YT - stmy(torg_id))/1000;     %in km

data_t= squeeze(data(:,idx,:));


return