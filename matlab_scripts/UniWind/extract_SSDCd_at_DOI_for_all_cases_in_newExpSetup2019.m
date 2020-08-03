%  This script is used to construct Cd/wind stress/frictional velocity
%  results in the same form to compare with Jame Edson's observation.
%  Requirement: data in NetCDF file has been retrieved already.
%               Aug 31: added deep water dataset to get reference Cd.
%  Xuanyu Chen; (:o_o:)
%  Aug 28, 2018;
%  updated for the new experiment setup: May 9, 2019;
%

%%% :o_o: global variables
global TRS_abspath LOC_abspath local_root scptdir trackdir
global  wspd  FLDID expID tailev_arry fetch_arry wspd_arry
global fieldFN specFN deps_of_interest slope
global fetch_arry_DW WDL_arry_SW

%%% :o_o: variables used locally:
FLDstr=num2str(FLDID);
infile1=[expID '_WW3spec_at_outpnts.mat'];
if length(tailev_arry)==1
    infile2=[expID '_FLD' FLDstr '_BGFDLcor_diagout.mat'];
else
    infile2=[expID '_FLD' FLDstr '_diagout.mat'];    
end
%%% add deep water result:
infile3='DW_WW3spec_at_outpnts.mat';
if length(tailev_arry)==1
    infile4=['DW_FLD' FLDstr '_BGFDLcor_diagout.mat'];
else
    infile4=['DW_FLD' FLDstr '_diagout.mat'];
end
%DW_path=[LOC_abspath(1:end-2) 'DW'];
if slope>=1000
    DW_path=[LOC_abspath(1:end-10) filesep 'DW_4km'];
else
    DW_path=[LOC_abspath(1:end-10+1) filesep 'DW_4km'];
    
end
%%% load deep water dataset:
DWpnt=load([DW_path filesep infile3]);
DWfld=load([DW_path filesep infile4]);
x0_dw=min(DWpnt.pntloc.x(:));
xN_dw=max(DWpnt.pntloc.x(:));

%deps_of_interest=[5 10 15 20 30 40]; % depths at which Cd will be extracted from the data record.
%fetch_arry=[50 100 200 400];

svdir=[LOC_abspath filesep 'matdata'];
if exist(svdir,'dir') ==0
    eval(['mkdir ' svdir]);
end

if length(tailev_arry)==1
    svfname=['FLD' FLDstr 'BGFDLcor_SSD_Cd_atSelectedDeps_slp' num2str(slope) '.mat'];
else
    svfname=['FLD' FLDstr 'SSD_Cd_atSelectedDeps_slp' num2str(slope) '.mat'];
end
%%% separate the deep water data retrieving and shoaling domain data retrieving.
%%% in case that the fetch in the deep water is not the same as those from the shoaling domain.

iWDL=0; 
for WDL=WDL_arry_SW
    WDLstr=num2str(WDL,'%3.3i');
    iWDL=iWDL+1;
    %%% load necessary input file:
    % note LOC_abspath here only to domainID.
   % work_path=[LOC_abspath filesep 'F' fetchstr 'km'];
    work_path=[LOC_abspath filesep 'DL' WDLstr 'km'];
    load([work_path filesep infile1]);    % check the path here, I need to specify fetch directory.
    load([work_path filesep infile2]); 
    icheck2=1;
    icheck=0;
    for id = 1:length(deps_of_interest);
        
        doi=deps_of_interest(id);
        % correlate dpt & slope to the location of this depth on the x
        % axis.
        xloc_dw=x0_dw+(doi/slope)/1E3;
        
                   
            kd_dw=DWpnt.spec.kp.*DWpnt.wvpar.dpt;
            kd_dw_DOI(id)=interp1(DWpnt.pntloc.x, kd_dw, xloc_dw);
      
    
        for itlev=1:length(tailev_arry)
            
            %% quantities in shallow water domain:
            %%% now find Cd_deep for reference: (a bit tricky part)
            Cd_sw=FLD_drv{itlev}.Cd(:,end);
            % use kp*d as the criterion.
            kd_sw=spec.kp(:,end).* wvpar.dpt(:,end);    % shallowness 
            
            %% quantities in deep water domain:
            % these are found by interpolation. 
            Cd_dw=DWfld.FLD_drv{itlev}.Cd;
            Cpi_dw=DWfld.FLD_drv{itlev}.Cp;
            
            Cd_dw_DOI(id,itlev)=interp1(DWpnt.pntloc.x, Cd_dw, xloc_dw);
            Cpi_dw_DOI(id,itlev)=interp1(DWpnt.pntloc.x, Cpi_dw, xloc_dw);
            
           
            % make plots to visualize for examination. 
            icheck=icheck+1;
            if icheck<=3
                figure(1)
                if slope==2000 || (slope==1000 && WDL==200)
                    plot(pntloc.x(:,end)+xN_dw-WDL, Cd_sw,'or');
                else
                    plot(pntloc.x(:,end),Cd_sw,'or');
                end
                hold on
                plot(DWpnt.pntloc.x, Cd_dw,'*b');
                hold off
                ylim([0 5])
                pause(0.5)
           end
%             pause
            
            %% these are all shallow water variables:
            % the ratio needs to be computed by interpolation also.
            % double check is needed.
            if slope==2000 || (slope==1000 & WDL==200)
                offset=xN_dw-WDL;
            
            else
                offset=0;
            end
            Cd_dw_interped=interp1(DWpnt.pntloc.x, Cd_dw, pntloc.x(:,end)+offset);
            
            
            if icheck2<=3
                figure(2);clf
                plot(DWpnt.pntloc.x,ones(size(DWpnt.pntloc.x)),'*b');
                hold on;
                plot(pntloc.x(:,end)+offset,ones(size(pntloc.x(:,end))),'or');
                pause(0.5)
                icheck2=icheck2+1;
            end
            
            Cd_ratio = Cd_sw./Cd_dw_interped;
            
            % extract results at depth of interest. Ues interpolation to
            % find instead. 
            Cd_ratio_DOI{id}(iWDL,itlev)=interp1(wvpar.dpt(:,end), Cd_ratio, doi);
            Cd_DOI{id}(iWDL,itlev)=interp1(wvpar.dpt(:,end),Cd_sw, doi);
            
            CpU10_DOI{id}(iWDL,itlev)=interp1(wvpar.dpt(:,end), FLD_drv{itlev}.CpU10(:,end),doi);
            wvage_DOI{id}(iWDL,itlev)=interp1(wvpar.dpt(:,end),FLD_drv{itlev}.wvage(:,end),doi);
  %          CpiU10(iWDL)= Cpi/wspd;        % not a function of tailev.
            
            
            
            
        end
    end
end


% %% DW data retrieving:
% for iWDL=1:length(WDL_arry_DW)
%     for itlev=1:length(tailev_arry)
%         % use fetch /distance to find the correct Cd. 
%         % The Cd outputed is every grid point. 
%        % Cd_dw_all(iWDL,itlev)=DWfld.FLD_drv{itlev}.Cd(iWDL);
%     end
% end


% save the results into a matlab file:
save([svdir filesep svfname],'Cd_ratio_DOI','CpU10_DOI','wvage_DOI','Cd_DOI');
