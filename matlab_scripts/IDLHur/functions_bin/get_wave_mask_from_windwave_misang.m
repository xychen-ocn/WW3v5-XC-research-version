% This is a function used to get the wave_mask for the IDLHur runs:
function wave_mask=get_wave_mask_from_windwave_misang(dth)
% Note: dth is the misalignment angle between wind and dominant wave. it is
% in RADIAN!!!

% 1. use dth_QCed_spatial to construct a wask mask for plotting
% wave_mask==1 -> wind sea (WS)
% wave_mask==2 -> Crossing Swell Positive (CSP)
% wave_mask==3 -> Crossing Swell Negative (CSN)
% wave_mask==4 -> Opposing Swell (OS)



% build wave mask:
wave_mask=NaN(size(dth));
wave_mask( cos(dth) >= sqrt(2)/2 ) = 1;    % Wind Sea
wave_mask( cos(dth) <  sqrt(2)/2 ) = 2;    % CSP
wave_mask( cos(dth) <  0         ) = 3;    % CSN
wave_mask( cos(dth) <= -sqrt(2)/2 ) =4;    % OS


end