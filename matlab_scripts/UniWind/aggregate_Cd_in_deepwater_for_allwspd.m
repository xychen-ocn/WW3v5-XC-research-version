% This script is used to extract Cd in deep water experiment:
% Xuanyu Chen
% Nov 15, 2018

%%% :o_o: global variables
global TRS_abspath LOC_abspath local_root scptdir trackdir GDrive_root
global  wspd  FLDID expID tailev_arry fetch_arry wspd_arry wind_fac
global fieldFN specFN 

%%% :o_o: variables used locally:
FLDstr=num2str(FLDID);
%%% add deep water result:
infile3='DW_WW3spec_at_outpnts.mat';
infile4=['DW_FLD' FLDstr '_diagout.mat'];


for iwspd=1:length(wspd_arry)
    wspdstr_2digit=num2str(wspd_arry(iwspd),'%2.2i');
    wdirstr=num2str(wdir);
    
    upper_loc=['wspd' wspdstr_2digit '_wdir' wdirstr];
    LOC_abspath=[GDrive_root filesep upper_loc];
    
    DW_path=[LOC_abspath filesep 'DW'];
    %%% load deep water dataset:
    DWpnt=load([DW_path filesep infile3]);
    DWfld=load([DW_path filesep infile4]);
    
    U10_allwspd(iwspd,:)=DWpnt.wvpar.wnd(:,end)./wind_fac;
    for itlev=1:length(tailev_arry)
        Cd_dw{itlev}(iwspd,:)=DWfld.FLD_drv{itlev}.Cd(:,end);
        wvage{itlev}(iwspd,:)=DWfld.FLD_drv{itlev}.wvage(:,end);
        mean_Cd_dw{itlev}(iwspd)=mean(Cd_dw{itlev}(iwspd,:),'omitnan');
        mean_Cd_dw2{itlev}(iwspd)=mean(Cd_dw{itlev}(iwspd,[3,6]),'omitnan');
    end
end


%% save data:
svfname=['FLD' FLDstr '_DeepWaterCd_allwspd.mat'];
save([GDrive_root filesep 'mat' filesep svfname],'Cd_dw','U10_allwspd', ...
        'wvage','mean_Cd_dw','mean_Cd_dw2');




