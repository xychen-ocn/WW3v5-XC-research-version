function [DataDir]=Results_Visualization_Wrapper(local_branch, exp_parms, exp_grps)
global mwstr
%% -------------- Section I --------------------- %%
%@@ 0. extract 
FLD_arry=exp_parms{1};
B_arry=exp_parms{2};
slpfac_arry=exp_parms{3};
TR_arry=exp_parms{4};
wind_fac=exp_parms{5};

%@@ experiment groups:
ExpGrps=exp_grps.main;
ExpSubGrps=exp_grps.sub;

if str2num(mwstr)==35
    DepParms.deps_q=[4000 25 20 15];
    
else
    DepParms.deps_q=[4000 40 30 15 10 4 2]; %[4000 15 10 5]; %
end
DepParms.dw_thres=500;
DepParms.sw_thres=50;
DepParms.dep_binwidth=0.1;                            % used to define a depth range within this binwidth.

%@@ 1. Define local parameters
InMatFile{1}='scattered_data_QCed_and_sorted_by_wave_category.mat';
InMatFile{2}='spatially_masked_quantites.mat';

SSDVarNames={'Cd'}; %'misang_wust'}; %'hs','lm','fp'};    %, 'fp'         %{'Cd','mss','misang'};
SSDVarUnits={'x1000'}; %{'m','m','s'};   %,'s'

% SSDVarNames={'misang_wust'};
% SSDVarUnits={'degree'}; 
% 

misang_type='misang_wust';
%SSDVarUnits={'x1000'};

StdWvParms={'hs','lm'};  %,'dirspr'};'hs','lm',
StdWvUnits={'m','m'}; %'m', 'm',
vecflags=[0, 1];
% define contour levels for different variables:
conlev.hs=[0:2:24];
conlev.lm=[0:50:500];
conlev.fp=[0:2:20];
conlev.wnd_mag=[30:10:60];
conlev.dirspr=[0:10:90];

conlev.Cd=[1:0.2:3 3.5 4 4.5 5];
conlev.mss=[0.001,0.005,0.01,0.05,0.1,0.15];

conlev.Cd_ratio=[0.5:0.1:2.0];
conlev.misang_wndwv=[-180:20:180];
%@@ 3. control variables in figure plotting:
isScatVis = true;        % flag controlling the visibility of scatters.
isRefData = false;        % flag for using RefData or not.


%@@ 2. build reference data structure:

local_root='/Users/xychen/Desktop/research_work/matlab_projects';

DragFormula_dir=[local_root filesep 'Drag_Formula'];
MSSFormula_dir=[local_root filesep 'MeanSquareSlope'];

RefData.Cd=load([DragFormula_dir filesep 'RefCd_dataset.mat']);
RefData.mss=load([MSSFormula_dir filesep 'Gleason_MSS.mat']);
RefData.misang=NaN;

% note that there is a bunch of observed Cd;
% need to pick no more than 3 to display as reference.
idx=3; %[1, 3];      % 1: Large and Pond, 3: COARE3.5



%% -------------- Section II --------------------- %%
% 1. set up information on different variables (used during figure plotting):
 %%@@: Cd
    Ref_xdata.Cd=RefData.Cd.Ref_U10data(:,idx); %Ref_U10data(:,idx);                % Observed_Cd.U10N;
    Ref_ydata.Cd=RefData.Cd.Ref_Cddata(:,idx).*1000;           % Observed_Cd.Cd * 1000
    jcnt=0;
    for j=idx
        jcnt=jcnt+1;
        RefData_ppt.tag{jcnt}.Cd=RefData.Cd.RefData_tag{j};
        RefData_ppt.color{jcnt}.Cd=RefData.Cd.RefData_color{j};
        RefData_ppt.lnsty{jcnt}.Cd=RefData.Cd.RefData_lnsty{j};
        RefData_ppt.lw(jcnt).Cd=RefData.Cd.RefData_lw(j);
        RefData_ppt.marker{jcnt}.Cd=RefData.Cd.RefData_marker{j};
    end
    Ylabelstr.Cd='Cd (x1000)';
    YMin.Cd=0;                                                 % or 0
    YMax.Cd=5;                                                 % or 4
    
    
    Ylabelstr.(misang_type)='Misalignment Angle';
    YMin.(misang_type)=-20;
    YMax.(misang_type)=20;
 %%@@: mean square slope  
 if 1==0
    Ref_xdata.mss=RefData.mss.U10;     % Observed_MSS.U10;
    Ref_ydata.mss=RefData.mss.MSS_Glson;    % Observed_MSS.MSS
    
    RefData_ppt.tag.mss='Gleason et al. 2018';
    RefData_ppt.color.mss=[0.45 0.45 0.45];
    RefData_ppt.lnsty.mss='-';
    RefData_ppt.lw.mss=1.2;
    RefData_ppt.marker.mss='';
    
    Ylabelstr.mss='Mean Square Slope';
    YMin.mss=0.001;
    YMax.mss=0.15;
    

 %%@@: misalignment angle:   
    Ref_xdata.misang=NaN;
    Ref_ydata.misang=NaN;
    
    RefData_ppt.color.misang=NaN;
    RefData_ppt.lnsty.misang=NaN;
    RefData_ppt.lw.misang=NaN;
    RefData_ppt.marker.misang=NaN;   
    RefData_ppt.tag.misang='';
 % end
 end
 






%% -------------- Section III --------------------- %%
% loop through experimental cases to make figures with functions:
% figure will have several subplots: 
%   subplot dimension (row = # of tail levels, column = # of deps)
%
if 1==1
for slpfac = slpfac_arry;
    slpfacstr=num2str(slpfac);
    
for TRID = TR_arry;
    TRIDstr=num2str(TRID);
       
    for FLDID=FLD_arry;
        FLDIDstr=num2str(FLDID);  
        
        if FLDID ==1
            YMin.(misang_type)=-5;
            YMax.(misang_type)=5;
            conlev.(misang_type)=[-2:0.5:2];
        else
            YMin.(misang_type)=-50;
            YMax.(misang_type)=50;
            conlev.(misang_type)=[-20:5:20];
        end
      
        for tailev=B_arry;
            %tailevstr=num2str(tailev,'%2.2i');
            tailevstr='GFDL';
            local_work_folder=['B' tailevstr '_FLD' FLDIDstr '_TR' TRIDstr '_angle135_new' ];  %'_FPI'
            
            DataDir.DW=[local_branch '/' local_work_folder '/DW'];  %/FPI_1.25x_3x
            DataDir.SW.main=[local_branch '/' local_work_folder '/SW_slp' slpfacstr ];
            
            for isub=1:length(ExpSubGrps);
                DataDir.SW.(ExpSubGrps{isub})=[DataDir.SW.main '/grd_' ExpSubGrps{isub}];
            end
            
            
            %%% 1. first, call a function to compute Cd Ratio:
          % if exist([DataDir.SW.main filesep 'CdRatio_atDifferentDeps_to_DW.mat'],'file')==0
                Cdts=compute_CdRatio(DataDir,InMatFile{2},DepParms);
                disp('finished finding CdRatio.');
         %end
        end
    end
end
end
end
         

%for angle=[45, 135]
 %   anglestr=num2str(angle);
    
 for slpfac = slpfac_arry;
     slpfacstr=num2str(slpfac);
     
     for TRID = TR_arry;
         TRIDstr=num2str(TRID);
         
         for FLDID=FLD_arry;
             FLDIDstr=num2str(FLDID);
             
             %         if FLDID ==1
             %             YMin.(misang_type)=-50;
             %             YMax.(misang_type)=50;
             %             conlev.(misang_type)=[-10:5:10];
             %         else
             YMin.(misang_type)=-30;
             YMax.(misang_type)=30;
             if FLDID==2
             conlev.(misang_type)=[-30:5:30];
             else
                 conlev.(misang_type)=[-30:1:30];
             end
             %         end
             
             for tailev=B_arry
                 %tailevstr=num2str(tailev,'%2.2i');
                 tailevstr='GFDL';
                 local_work_folder=['B' tailevstr '_FLD' FLDIDstr '_TR' TRIDstr '_angle135_new'];  %'_FPI'_angle' anglestr
                 
                 DataDir.DW=[local_branch '/' local_work_folder '/DW'];  %/FPI_1.25x_3x
                 DataDir.SW.main=[local_branch '/' local_work_folder '/SW_slp' slpfacstr ];
                 
                 for isub=1:length(ExpSubGrps);
                     DataDir.SW.(ExpSubGrps{isub})=[DataDir.SW.main '/grd_' ExpSubGrps{isub}];
                 end
                 
                 %%% 2. make scatter plots and spatial plots of sea-state dependent variables.
                 % close all;
                 if 1==0
                 % call function to make figures
                 % variation of slope:
                 make_IDLHur_scatter_plots(DataDir, InMatFile{1},SSDVarNames, slpfac, DepParms, isScatVis, isRefData, ...
                     Ref_xdata, Ref_ydata, RefData_ppt,Ylabelstr, YMin, YMax);
                 %
                 % variation of angle of attack:
                 %             make_IDLHur_scatter_plots(DataDir, InMatFile{1},SSDVarNames, angle, DepParms, isScatVis, isRefData, ...
                 %                 Ref_xdata, Ref_ydata, RefData_ppt,Ylabelstr, YMin, YMax);
                 % % % % %
                 end
                 
                 
                 if 1==1
                     % make plot for tail level dependent parameters (Cd, misang, mss):
                     if length(SSDVarNames)==1 && strcmp(SSDVarNames{1},'Cd')
                         make_IDLHur_SpaceTime_plots_for_Cd(DataDir, InMatFile{2}, SSDVarNames, SSDVarUnits, slpfac, DepParms, conlev);
                         %     make_IDLHur_SpaceTime_plots_for_Cd(DataDir, InMatFile{2}, SSDVarNames, SSDVarUnits, angle, DepParms, conlev);
                         
                     elseif strcmp(SSDVarNames, 'misang_wust');
                         make_IDLHur_SpaceTime_plots_for_misang(DataDir, InMatFile{2}, SSDVarNames,StdWvUnits,DepParms,conlev,slpfac);
                     else
                         make_IDLHur_SpaceTime_plots_for_SSD_variables(DataDir, InMatFile{2}, SSDVarNames, SSDVarUnits, tailev, DepParms, conlev);
                     end
                     
                 end
                 
             end     % tail level
             
         end          % methods switch
         
     end             % triad switch
     
 end             % slope

%end             % angle.

if 1==0
%close all;
%@@ make spatial plots (in space and time) for standard wave parameters:;
FLDIDstr='1'; 
tailevstr='GFDL';
for slpfac = slpfac_arry;
    slpfacstr=num2str(slpfac);
    for TRID = TR_arry;
        TRIDstr=num2str(TRID);
        
        local_work_folder=['B' tailevstr '_FLD' FLDIDstr '_TR' TRIDstr '_smooth'];
        
        DataDir.DW=[local_branch '/' local_work_folder '/DW'];
        DataDir.SW.main=[local_branch '/' local_work_folder '/SW_slp' slpfacstr];
        
        for isub=1:length(ExpSubGrps);
            DataDir.SW.(ExpSubGrps{isub})=[local_branch '/' local_work_folder ...
                '/SW_slp' slpfacstr '/grd_' ExpSubGrps{isub}];
        end
        
        % make plot
        %close all;
      make_IDLHur_SpaceTime_plots_for_standard_wave_parameters(DataDir, InMatFile{2}, StdWvParms,StdWvUnits,DepParms,conlev,slpfac,vecflags);
      
    end
end
end

%% ----  End of function ----- %%
return
