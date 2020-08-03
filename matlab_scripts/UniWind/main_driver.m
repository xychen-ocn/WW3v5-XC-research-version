% Purpose: this script is used to drive data processing and figure making
%          for the Uniform Wind Experiments.
% Time   : Aug 23rd, 2018
% Coder  : Xuanyu Chen; 
%

clear all; close all; clc;

%% define global variables:
%%% pathes:
global TRS_abspath LOC_abspath local_root scptdir trackdir GDrive_root

%%% looping variables:
global wspd FLDID fetch expID wdir deps_of_interest wind_fac
global wspd_arry FLD_arry tailev_arry fetch_arry dep_selected 
%%% others
global fieldFN specFN slope fetch_arry_DW WDL_arry_SW  TRID slp


%% assign values to global variables:
casename='UniWnd';  %UniWnd; IDLHur; swell_gen
%casename='UniWnd';
TRS_url_root=['http://tds.renci.org:8080/thredds/dodsC/dhs-crc-unc/xychen/work_ST4' filesep casename];
GDrive_root=['/Volumes/Cumulonimbus 1/Research_backup/matlab_projects/ssd_shallow_water' filesep casename];
%GDrive_root='/Users/xychen/Desktop/research_work/matlab_projects/shallow_water_ssd_exp/analysis/results/GDrive/UniWnd';
local_root='/Users/xychen/Desktop/research_work/matlab_projects';
%datadir='';
scptdir=[local_root filesep 'shallow_water_ssd_exp/analysis/script/matlab/UniWind_Exps'];

wspd_arry=[10:5:65]; %[15 35 65]; %[20 25 35 40 45]; %[15 30 50]; %; %5:5:50; [10 30 50];[20 25 35 40 45];
slp_arry=[100 200 400 1000 2000]; %[200 400 1000 2000]; %[400 1000]; %[100 200 400 1000 2000];
wdir=180;
fetch_arry=[200 400 600]; %[200 400 600]; %[25 50 100 200 400 750]; %[100 750]; %[25 50 100 200 400 750];
tailev_arry=1; %[2 6 12];
FLD_arry=[1 2];
srcterm='ST4';
wind_fac=1/1.14;
%srcterm={'ST6','ST4','ST2'};%,'ST2'

expID='SW';   %DW or SW;
%expID='cntrl';   % or exp

fieldFN='ww3.201711.nc';
specFN='ww3.201711_spec.nc';


addpath('/Users/xychen/Desktop/research_work/matlab_projects/shallow_water_ssd_exp/analysis/script/matlab/UniWind_Exps');

%%% 1. Process data: (either DW or SW data) 
process_scrptn='read_ncfiles_single_case.m';

%%% 2. make figures of different kinds for each wind speed group.
% ------------------------------------------------------------------- %
%   Types: 1. 2D spatial map at steady state;
%        : 2. 1D cross-section (Hs, LM, mss, wvslp, cd, wind stress);
%        : 3. wind profiles and stress profiles;
%        : 4. form drag composition.
%        : 5. saturation spectrum.
% ------------------------------------------------------------------- %
plotting_scptn{1}='spatial_map_of_standard_output_variables';
plotting_scptn{2}='alongshore_cross_section_of_variables';
plotting_scptn{3}='profiles_of_stress_and_wind';
plotting_scptn{4}='form_drag_wavenumber_components';
plotting_scptn{5}=['compare_saturation_spectrum_between_srceterms_' expID];
plotting_scptn{6}='plot_2DWaveSpec_and_GrowthRate.m';

%% Process data:
for TRID=[1]
    % ids=[1 3 4];
    for slp=slp_arry
        if strcmp(expID, 'SW')
            SW_flag=1;
            TRS_domain_name=[expID '_slp' num2str(slp)];
            % TRS_domain_name=['slp' num2str(slp) ];
        else
            SW_flag=0;
            TRS_domain_name=[expID '_4km_DefST4'];
            % TRS_domain_name=['slp' num2str(slp) filesep 'cntrl_4km'];
        end
        %  TRS_domain_name=['slp' num2str(slp) filesep 'cntrl'];
        ext='';
        %ext='';
        for wspd = wspd_arry
            wspdstr_2digit=num2str(wspd,'%2.2i');
            wdirstr=num2str(wdir);
            upper_TRS = ['wspd' num2str(wspd) '_wdir' wdirstr filesep TRS_domain_name ext];
            upper_loc=['wspd' wspdstr_2digit '_wdir' wdirstr filesep TRS_domain_name ext];
            if strcmp(expID,'DW')
                TRS_abspath=[TRS_url_root filesep upper_TRS];
                LOC_abspath=[GDrive_root filesep upper_loc];
                TRS_abspath=LOC_abspath;
                if exist(LOC_abspath,'dir') == 0
                    eval(['mkdir ' LOC_abspath]);
                end
                % TRS_abspath=LOC_abspath;
                %%% run script:
                run(process_scrptn);
            else   % SW:
                for fetch = fetch_arry
                   %                     end
                    fetchstr_3digit=num2str(fetch,'%2.2i');
                    fetchstr_2digit=num2str(fetch);
                  
                    folder_loc=['DL' fetchstr_3digit 'km'];  % DefST4 %_ww3v6
                    
                    %expstr={'A','B','C'};
                    %        'indx1.0_misalgn','indx1.0_misalgn270'};
                    %  swell_cases={'3_misalgn240','3_misalgn270'};
                    %  swell_cases={'windsea','windsea_misalgn','windsea_misalgn270'};
                    %            for ic=1:2 %:length(swell_cases);
                    %  %                casename=swell_cases{ic};
                    %   for iexp=1:length(expstr)
                    %folder_loc=['F750km_wspd30_HurST4_TR1_swell' casename filesep 'probe_exp/Exp' expstr{iexp}];%indx1.0_misalgn270';
                    % folder_loc='swell2_TR1';
                    % folder_loc=['F750km_TR1_' casename];
                    
                    %folder_loc=['DL' fetchstr_3digit 'km' filesep 'probe_exp/stpness_v2'];
                    %folder_loc=[filesep 'F' fetchstr_3digit 'km_TR0_10km'];
                    
                    LOC_abspath=[GDrive_root filesep upper_loc filesep folder_loc];
                    if exist(LOC_abspath,'dir') == 0
                        eval(['mkdir ' LOC_abspath]);
                    end
                    TRS_abspath =LOC_abspath;
                    %%% run script:
                    disp(['Working on: wspd=' wspdstr_2digit 'm/s, Incoming Fetch=' fetchstr_2digit 'km...']);
                    % run(process_scrptn);
                    %% Quality control and find peak:
                    % the following code require visual examination.
                     run('find_peakCdRatio_locations_FPIresults.m');
                    % end
                    %            end
                end       % end fetch loop
            end
            %%% calling scripts:
            % run(scrpt1);
            % run(scrpt2);
            % run(scrpt3);
            % run(scrpt4);
        end
    end
end




if 1==0
%% collapse all data:
for FLDID=FLD_arry
    run('collapse_all_data.m');
end

for FLDID=FLD_arry
    run('aggregate_Cd_in_deepwater_for_allwspd.m');
end

%% extract SSDCd for each wind speed run:
deps_of_interest=[10:5:50]; % depths at which Cd will be extracted from the data record.

fetch_arry_DW=fetch_arry;
WDL_arry_SW=fetch_arry;


for slope=slp_arry(1); %;
    
    for wspd=wspd_arry(12);
        wspdstr_2digit=num2str(wspd,'%2.2i');
        wdirstr=num2str(wdir);
        
        ext=['_slp' num2str(slope)];
        upper_loc=['wspd' wspdstr_2digit '_wdir' wdirstr filesep expID ext ];
        LOC_abspath=[GDrive_root filesep upper_loc];
        
        for FLDID = FLD_arry
            run('extract_SSDCd_at_DOI_for_all_cases_in_newExpSetup2019.m');
        end
    end
end

% select a depth to make plot:

for slope=slp_arry(1); %(1:end)
    dcnt=0;    
    run('aggregate_CdRatio_at_allwspd_and_make_plots_newExpSetup2019.m');
    %pause(1)

end


%% plot 2-dim wave spectrum:
for wspd=wspd_arry(5);
    wspdstr_2digit=num2str(wspd,'%2.2i');
    wdirstr=num2str(wdir);
    
    upper_loc=['wspd' wspdstr_2digit '_wdir' wdirstr];
    LOC_abspath=[GDrive_root filesep upper_loc];

    for fetch=fetch_arry(5);
        
        run(plotting_scptn{6});
    end
end

end

