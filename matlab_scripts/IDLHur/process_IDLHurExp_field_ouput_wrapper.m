function process_IDLHurExp_field_ouput_wrapper(remote_branch, local_branch, exp_parms, exp_grps, ...
                                               time_start,time_end)
% Description:
%
%
%


% ============== Start of function ============ %
%% ----------- Section I -------------- %%
%@@ define I/O variables:
InpFN_core='ww3.201711';
InpFN_ext='.nc';
DW_folder='DW';

%@@ experiment groups:
ExpGrps=exp_grps.main;
ExpSubGrps=exp_grps.sub;

%@@ following serve as input to the actual data processing function.
process_freq.xinc=1;
process_freq.yinc=1;
process_freq.tinc=1;


DepParms.deps_q=[4000 40 30 15];
DepParms.dw_thres=500;
DepParms.sw_thres=50;
DepParms.dep_binwidth=0.1;                            % used to define a depth range within this binwidth.



% OutVars={'hs','dpt','dth','pwdir','dirspr','wvage','wnd','ust','mss_tot',...
%     'mss_measured','misang_wust','cpU10','Cd','dirspr'};  

MatFname{1}='spatially_masked_quantites.mat';
%MatFname{1}='spatially_masked_quantites_trackonly.mat';
MatFname{2}='scattered_data_QCed_and_sorted_by_wave_category.mat';

% release exp. parmeters:
FLD_arry=exp_parms{1};
B_arry=exp_parms{2};
slpfac_arry=exp_parms{3};
TR_arry=exp_parms{4};
wind_fac=exp_parms{5};
tspd=exp_parms{6};
tspdn=['tspd' num2str(tspd)];
mwstr=exp_parms{7};

%TrackInfo.trackFN=['IDL_mw65_rmw70_' tspdn '_center.out'];  % mw and rmw value doesn't matter.
TrackInfo.trackFN='stmcenter_cartesian_and_time_angle135.mat';
TrackInfo.dir=[local_branch filesep 'trackinfo'];


%% ------------- Section II ----------- %%
for TRID = TR_arry;
    TRIDstr=num2str(TRID);
       
    for FLDID=FLD_arry;
        FLDIDstr=num2str(FLDID);
        
        for tailev=B_arry;
            if length(B_arry)==3
                tailevstr=num2str(tailev,'%2.2i');
            elseif length(B_arry)==1 && tailev==1
                tailevstr='GFDL';
                %tailevstr='ADCIRC';
            end
           % remote_work_folder=['FLD' FLDIDstr '_B' tailevstr '_TR' TRIDstr '_smooth']; 
         %  remote_work_folder=['FLD' FLDIDstr '_' tailevstr '_TR' TRIDstr '_smooth'];
           
            local_work_folder=['B' tailevstr '_FLD' FLDIDstr '_TR' TRIDstr '_angle135_new'];
            remote_work_folder=local_work_folder;
            disp(local_work_folder);
            for ig=1:length(ExpGrps)
                expID=ExpGrps{ig};
                disp(expID);
                switch expID
                    
                    case 'DW'                     
                        %local_work_folder=['B' tailevstr '_FLD' FLDIDstr '_TR' TRIDstr '_angle45'];

                        % build up the path to input data in the remote computer:
                        InpFN=[InpFN_core InpFN_ext];
                        remote_dir=[remote_branch filesep DW_folder filesep remote_work_folder];
                        absFN=[remote_dir filesep InpFN];
                        
                        % build up local path for data storage:
                        local_dir=[local_branch filesep local_work_folder filesep expID];
                        if strcmp(remote_branch, local_branch)
                            absFN=[local_dir filesep InpFN];
                        end

                        process_time_range={time_start.DW; time_end.DW.(tspdn)};
                        process_freq.x=1;
                        process_freq.y=1;
                        process_freq.t=1;
                        
                      % if ~exist([local_dir filesep MatFname{1}],'file')
                         %   call the true function for data processing
                            process_IDLHur_field_output(absFN, process_time_range, process_freq, ...
                                wind_fac, DepParms, MatFname, local_dir, TrackInfo, expID,tailev);
                      % end
                    case 'SW'
                        for slpfac=slpfac_arry
                            slpfacstr=num2str(slpfac);
                            slpn=['slp' slpfacstr];
                            for isub=1:length(ExpSubGrps)
                                subgrd=ExpSubGrps{isub};
                                % build up the path to input data in the remote
                                % computer:
                                InpFN=[InpFN_core '_' subgrd InpFN_ext];
                                SW_folder=[expID '_slp' slpfacstr '_multigrd'];
                                subgrd_folder=['grd_' subgrd];
                                remote_dir=[remote_branch filesep SW_folder filesep ...
                                    remote_work_folder filesep  subgrd_folder];
                                absFN=[remote_dir filesep InpFN];
                                
                                % build up local path for data storage:
                               local_dir=[local_branch filesep local_work_folder filesep ...
                                   expID '_slp' slpfacstr filesep subgrd_folder];

%                                local_dir=[local_branch filesep local_work_folder filesep ...
%                                          'DW_slp' slpfacstr filesep subgrd_folder];
%                                 
                                if strcmp(remote_branch, local_branch);
                                    absFN=[local_dir filesep InpFN];
                                end
                                
                                process_time_range={time_start.SW; ...
                                                time_end.SW.(tspdn).(slpn)};
                                process_freq.x=1;
                                process_freq.y=2;
                                process_freq.t=1;
                                
                                disp(subgrd);
                               %if ~exist([local_dir filesep MatFname{1}],'file')
                                    % call the data processing function
                                   process_IDLHur_field_output(absFN, process_time_range, process_freq, ...
                                       wind_fac,  DepParms,  ...
                                       MatFname, local_dir, TrackInfo, expID,tailev);
                                %end
                            end    % loop for subgrid
                            
                        end        % loop for slope
                        
                end          %  End switch
                
            end              % end for Exp. Groups
            
        end     % end for tail level
        
    end           % end for FLD
    
end           % end for Triad interaction


% End of Wrapper

return
          
