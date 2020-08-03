% ================================================================================================ %
% Purpose: 
%   This script define the experimental parameters for IDLHur experiments and the associated pathes;
%   This shell will call function/scripts to perform data extraction, processing and visualization ;
%   In this latest version, global variables are abandoned.
%   Instead, functions are built for similar purpose.
% Creator:  Xuanyu Chen;
% Date:  2019/1/18

clear all; clc; close all;
global mwstr

%============ Secion I =====================%
%@@---------- Define pathes: -------------@@%
TRS_url_root='http://tds.renci.org:8080/thredds/dodsC/dhs-crc-unc/xychen/work_ST4/IDLHur';
GDrive_root='/Volumes/Cumulonimbus 1/Research_backup/matlab_projects/ssd_shallow_water/IDLHur';
%GDrive_root='/Users/xychen/Desktop/research_work/matlab_projects/shallow_water_ssd_exp/analysis/results/GDrive';
local_root='/Users/xychen/Desktop/research_work/matlab_projects';
scptdir=[local_root filesep 'shallow_water_ssd_exp/analysis/script/matlab/IDLHur_Exps'];
addpath([scptdir filesep 'toward_final/update/functions_bin']);
%@@---------- Define experimental parameters: -------------@@%
tspd_arry=[10];                     % [5 10]
FLD_arry=[1];                       % [1,2]
B_arry=1; %[2 6 12];                         %[2 6 12];
slpfac_arry=200; %[200,2000];                    %[200, 2000];
TR_arry=1;                          %[1 0]
mwstr='65';                         % maximum wind speed.

  %@@@@ some relevant parameters for calling different functions
    wind_fac=1.0/1.14 ;                 % or 1.0;  gust-factor.

%@@ experiment groups:
ExpGrps.main={'DW','SW'};
ExpGrps.sub={'coarse','fine'};

% start time and end time of simulation.
time_start.DW='20171125 180000';
%time_end.DW ='20171126 180000';

time_start.SW='20171125 180000';

% for different angle of approach:
% time_end.SW.tspd10.slp2000='20171201 180000'; %'20171128 130000';


time_end.SW.tspd5.slp2000 ='20171201 073000';
time_end.SW.tspd10.slp2000='20171128 130000'; %'20171128 130000';

time_end.SW.tspd5.slp200 ='20171129 060000 ';
time_end.SW.tspd10.slp200='20171128 000000';

time_end.DW.tspd5='20171201 073000';    %'1130 000000'
time_end.DW.tspd10='20171128 000000';   %'1128 000000'

%time_end.DW.tspd5='20171127 000000';

%============= Section II ====================%
%@ add pathes to functions
addpath([scptdir filesep 'toward_final/function_bin']);

 close all;       
for tspd=tspd_arry
    tspdstr=num2str(tspd);
    
    %TRS_branch=[TRS_url_root filesep 'mw' mwstr '_rmw70_tspd' tspdstr];
    LOC_branch=[GDrive_root filesep 'mw' mwstr '_rmw70_tspd' tspdstr];
    TRS_branch=LOC_branch;
    
    exp_parms={FLD_arry; B_arry; slpfac_arry; TR_arry; wind_fac; tspd; mwstr};                        %[2 6 12];

    % call wrapper function:
    %@@ ------- Procedure 1: data extraction and processing ------- @@%
    % loop over the experiment cases to perform the follow procedures  %
    
    % note: make sure to change the folder directory name.
  % process_IDLHurExp_field_ouput_wrapper(TRS_branch, LOC_branch, exp_parms, ExpGrps, time_start,time_end);
    
    
    %@@ ------- Procedure 2: data visualization ----------- @@%
    % looping will be determined by the need of my plots.
    % option: 1. scatter (1D)
    % optino: 2. spatial (2D)
    % note: make sure to change the folder directory name.
    DataDir= Results_Visualization_Wrapper(LOC_branch, exp_parms, ExpGrps);
    
end




