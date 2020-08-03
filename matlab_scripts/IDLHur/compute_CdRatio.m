% This script is used to compute the ratio between Cd in the shallow water
% and that in the deep water. (to bulk Cd and also SSD Cd in the deep
% water.)
% Xuanyu Chen;
% Jun 15, 2019;
% 
function Cdts=compute_CdRatio(DataDir, MatFname,DepParms)

dw_thres=DepParms.dw_thres;
sw_thres=DepParms.sw_thres;
deps_q=DepParms.deps_q;

DragDir='/Users/xychen/Desktop/research_work/matlab_projects/Drag_Formula';
% 1. load DW dataset:
DW=load([DataDir.DW filesep MatFname]);

% 2. load SW dataset:
SWC=load([DataDir.SW.coarse filesep MatFname]);
SWF=load([DataDir.SW.fine filesep MatFname]);

meanCdFname='GFDL_Cd.mat';
bulkchoice='mean'; 
svdataFN='CdRatio_atDifferentDeps_to_DW.mat';

Ref_Cdbulk=load([DragDir filesep meanCdFname]);
Ref_Cdbulk.Cd=Ref_Cdbulk.Cd*1000;

varn='Cd';

% 3. compute the ratio at request depth.;
for id=1:length(deps_q);
    
    d = deps_q(id);
    
    % load in different data set based on the requested depth:
        if (d>=dw_thres) 
           indata=DW.QCed_fieldVar;
           coord=DW.coord;
           stmx=DW.stmx./1E3;
           stmy=DW.stmy./1E3;
           
        elseif (d>=sw_thres)   % load coarse resolution data.
           indata=SWC.QCed_fieldVar;
           coord=SWC.coord;
           stmx=SWC.stmx./1E3;
           stmy=SWC.stmy./1E3;
           stmtime=SWC.stmtime;
           ww3_timenum=SWC.ww3_timenum;
          [tmp,dpt_maxtrix,tmp,tmp]=check_qtty_in_netCDFfile([DataDir.SW.coarse filesep 'ww3.201711_coarse.nc'],0,'dpt');

           
        elseif (d<sw_thres)    % load fine resolution data.
           indata=SWF.QCed_fieldVar;
           coord=SWF.coord;
           stmx=SWF.stmx./1E3;
           stmy=SWF.stmy./1E3;
           stmtime=SWF.stmtime;
           ww3_timenum=SWF.ww3_timenum;
           [tmp,dpt_maxtrix,tmp,tmp]=check_qtty_in_netCDFfile([DataDir.SW.fine filesep 'ww3.201711_fine.nc'],0,'dpt');

           
        end      
       
       ww3_time_sw=SWC.ww3_timenum;     % the time for the SW run is the same regardless of the resolution.
       ww3_time_dw=DW.ww3_timenum;      %

      %% ------ generate the sudo-space--space 2D matrix: ------ %%
       DWdata.data=DW.QCed_fieldVar.Cd;
       DWdata.XX=DW.coord.XX;
       DWdata.YY=DW.coord.YY;
       % display information:
       disp(['DW NT=' num2str(size(DWdata.data,3))]);
       disp(['DW dt=' num2str((DW.ww3_timenum(2)-DW.ww3_timenum(1))*60)]);
       
       disp(['SWC NT=' num2str(length(SWC.ww3_timenum))]);
       disp(['SWC dt=' num2str((SWC.ww3_timenum(2)-SWC.ww3_timenum(1))*60)]);
       
       disp(['SWF NT=' num2str(length(SWF.ww3_timenum))]);
       disp(['SWF dt=' num2str((SWF.ww3_timenum(2)-SWF.ww3_timenum(1))*60)]);

      %  size(indata.(varn),3)
       
         if d<dw_thres
             [ny,nx,nt]=size(indata.dpt);
             dpt_vec=dpt_maxtrix(round(ny/2),:,1);
                
%              size(indata.(varn),3)
%              size(DWdata.data,3)
%              pause
             [data_t, Cd2deep_ratio, TT,YT,coordN]=find_CdRatio_through_TimeSpaceConversion(indata.(varn), DWdata, ...
                 coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);

            % compute ratio:
            % interpolation from the DW grid to the grid of the shoaling domain:
            % this assumes that the resolution did not impact the Cd. 
            Cd_bulk_interp=interp2(double(coordN_DW.XX),double(coordN_DW.YY),Cd_bulk_t, ...
                double(coordN.XX),double(coordN.YY));
            Cd2bulk_ratio= data_t./Cd_bulk_interp;

           Cdts.(['D' num2str(d)]).value=data_t;
           Cdts.(['D' num2str(d)]).ratio2deepSSD=Cd2deep_ratio;
           Cdts.(['D' num2str(d)]).ratio2bulk=Cd2bulk_ratio;
           Cdts.(['D' num2str(d)]).coord=coordN;


         else  
            % It would be better to sync the deep water results with the
            % shallow water results. 
            % I need to know where the TC is in the shallow water...
            [tmp,sid]=min(abs(ww3_time_dw-ww3_time_sw(24)));
            data_t=indata.(varn)(:,:,sid);
            wnddata_t=indata.wnd_mag(:,:,sid);
            
            coordN.XX = (coord.XX - stmx(sid)*1000)./1000;
            coordN.YY= (coord.YY - stmy(sid)*1000)./1000;
            
                        
            %@ 3. additional computation
            Cd_bulk_t=compute_Cd_bulk_DW(Ref_Cdbulk,wnddata_t);

            coordN_DW=coordN;
            % compute the ratio also in the deep water domain.
            Cd2bulk_ratio= data_t./Cd_bulk_t;
           
           Cdts.(['D' num2str(d)]).value=data_t;
          % Cdts.(['D' num2str(d)]).ratio2deepSSD=Cd2deep_ratio;
           Cdts.(['D' num2str(d)]).ratio2bulk=Cd2bulk_ratio;
           Cdts.(['D' num2str(d)]).coord=coordN;

           Cdts.(['D' num2str(d)]).wnd=wnddata_t; 
           
        end
     
   


end


%if exist([DataDir.SW.main filesep svdataFN],'file')==0
    % results will be saved into a matlab dataset:
    save([DataDir.SW.main filesep svdataFN],'Cdts');
%end

return