% This function is used to convert Hurricane track information (location and time)
% from ASCII to matlab file, and from spherical to cartesian.
% X.C. Jan 2, 2019
%
function [stmx,stmy,datenum_stm]=read_ASCII_trackInfo(trackdir,trackfile)
% input: trackdir, trackfile (string)
% output: 
% stmx, stmy: the center location of the storm in cartesian coordinate
% datenum_stm: the time number of the storm. the frequency is the same as
% the wind forcing.
%
%%% get TC center information, and convert spherical coordinate to cartesian
%%% coordinate.
center=load([trackdir filesep trackfile]);
stmlon=center(:,2);
stmlat=center(:,1);

x0=4000-(140-stmlon(1))*100;
stmx=x0+(stmlon-stmlon(1))*100;

y0=1800-(23-stmlat(1))*100;
stmy=y0+(stmlat-stmlat(1))*100;

center_time=load([trackdir filesep 'date.out']);
stmhrly=center_time(1:4:end,:);
for i=1:length(stmhrly)
    stmdatestr_hrly(i,:)=[num2str(stmhrly(i,1)) num2str(stmhrly(i,2),'%6.6i')];
end
for i=1:length(center_time);
    stmdatestr(i,:)=[num2str(center_time(i,1)) num2str(center_time(i,2),'%6.6i')];
end
datenum_stmhrly=datenum(stmdatestr_hrly,'yyyymmddHHMMSS');
datenum_stm=datenum(stmdatestr,'yyyymmddHHMMSS');

save([trackdir filesep 'stmcenter_cartesian_and_time.mat'],'stmx','stmy','datenum_stm', 'datenum_stmhrly');
return