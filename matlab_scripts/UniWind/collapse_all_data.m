% This script is used to plot all data into a dimensionless space:
% Cd/Cd_deep vs. kp*d or (kpi*d)
% Xuanyu Chen (:>o_o<:)
% Nov 10, 2018

%%% :o_o: global variables
global TRS_abspath LOC_abspath local_root scptdir trackdir GDrive_root
global wspd_arry fetch_arry expID tailev_arry wdir dep_selected FLD_arry
global fieldFN specFN FLDID

%%% :o_o: variables used locally:
FLDstr=num2str(FLDID);
% shallow water matlab files:
infile1=[expID '_WW3spec_at_outpnts.mat'];
infile2=[expID '_FLD' FLDstr '_diagout.mat'];
%%% add deep water result:
infile3='DW_WW3spec_at_outpnts.mat';
infile4=['DW_FLD' FLDstr '_diagout.mat'];
%DW_path=[LOC_abspath(1:end-2) 'DW'];
svdir='/Users/xychen/Desktop/research_work/matlab_projects/shallow_water_ssd_exp/analysis/results/test_collapsed_data';
jetcolor=jet(50);
color=jetcolor(1:5:50,:);
color=color(2:end,:);

grav=9.8;
% attempt to use wave age Cp/U10 or Cp/U* to separate the data:
% Cp/U* failed to group data. 

%%% structure:
cnt=0; cnt2=0;
Cd_ratio_all=NaN(200,162);
kpd_all=[];
kpid_all=[];
wvage_deep=[];
CpU10_deep=[];
X=[];
pos=customize_subplot_size(9,6,0.01,0.01);
for iwspd=1:length(wspd_arry)    % different wind speed
    
    wspdstr_2digit=num2str(wspd_arry(iwspd),'%2.2i');
    wdirstr=num2str(wdir);
    upper_loc=['wspd' wspdstr_2digit '_wdir' wdirstr];
    
    %% load DW data:
    DW_path=[GDrive_root filesep upper_loc filesep 'DW'];
    %%% load deep water dataset:
    DWpnt=load([DW_path filesep infile3]);
    DWfld=load([DW_path filesep infile4]);
    
    for ifetch=1:length(fetch_arry)     % different fetch
        fetchstr=num2str(fetch_arry(ifetch),'%3.3i');
        SW_path=[GDrive_root filesep upper_loc filesep 'SW' filesep 'F' fetchstr 'km'];
        %% load SW data
        load([SW_path filesep infile1]);
        SWfld=load([SW_path filesep infile2]);
        
        cnt2=cnt2+1;
        for itlev=1:3      % different tail levels
            cnt=cnt+1;
            Cd_shallow=SWfld.FLD_drv{itlev}.Cd(:,end);
            Cd_shallow(Cd_shallow==0)=NaN;
            dep_shallow=wvpar.dpt(:,end);
            kp_shallow=spec.kp(:,end);           
            
            Cd_deep = DWfld.FLD_drv{itlev}.Cd(ifetch);
            wvage_deep(cnt)= DWfld.FLD_drv{itlev}.wvage(ifetch);
            CpU10_deep(cnt)=DWfld.FLD_drv{itlev}.CpU10(ifetch);
            kp_in= DWpnt.spec.kp(ifetch);
            
            % Cd ratio:
            Cd_ratio = Cd_shallow./Cd_deep;
            kpd=kp_shallow.*dep_shallow;
            kpid=kp_in.*dep_shallow;
            
            %% apply mask to select only data at which the kpd is decreasing.
            mask_criterion=diff(kpd)./diff(dep_shallow);
            mskid=find(mask_criterion<0);
            if length(mskid)>1
                %for im=1:length(mskid)
                valid= find( dep_shallow(mskid)<50);
                mskid=mskid(valid);
            end
            kpd_masked=kpd;
            kpd_masked(mskid:end)=NaN; 
           % Cd_ratio(mskid:end)=NaN;
            kpid(mskid:end)=NaN;
            
            
            kpd_all(1:length(kpd_masked),cnt)=kpd_masked;
            kpid_all(1:length(kpid),cnt)=kpid;
            Cd_ratio_all(1:length(Cd_ratio),cnt)=Cd_ratio;
            X(1,cnt)=grav.*fetch_arry(ifetch).*1000/wspd_arry(iwspd).^2;

            
            
            %% make plot
            if 1==0
            figure(1);  
            subplot(9,6,cnt2)% local kd.
            plot(kpd_masked,Cd_ratio,'-','color',color(iwspd,:)); hold on;
            ylim([0 2])
            xlim([0 2])
            axis('square');
            set(gca,'pos',pos{cnt2});
            end
%             
%             subplot(1,2,2)
%             plot(kpid,Cd_ratio,'-','color',color(iwspd,:)); hold on;
%             ylim([0 2])
%             xlim([0 2])
%             axis('square');
            % plot(Cd_ratio,kpi);
            
            
            

            
            
        end
        %pause
        
    end
    
end
%
%figure(2)
%hist(CpU10_deep,5);
xbin=[0.2:0.2:1.2];
%% regroup:
% use waveage to regroup data:
group_crit=CpU10_deep;    %log10(X)
%group_crit=log(X);    %log10(X)

% figure(4)
% h=histogram(group_crit);
%xbin=[5:5:35];
xbin=[0.2:0.2:1.2];
%xbin=[1.5:0.5:5];
%xbin=[4:1:12];

% Group 1: wave age <10
clear Gy Gx Gx2 
for i=1:length(xbin)-1
    range_mask=zeros(size(group_crit));
    range_mask(group_crit<=xbin(i+1))=1;
    range_mask(group_crit<=xbin(i))=0;
    IDs=find(range_mask==1);
    Gy{i}=Cd_ratio_all(:,IDs);
    Gx{i}=kpd_all(:,IDs);
    Gx2{i}=kpid_all(:,IDs);
    
end


figname=['FLD' num2str(FLDID) '_collapsed_data_sorted_by_wave_age_CpU10.jpg'];
figure(6)
for ig=1:length(Gx);
    subplot(3,3,ig);
    % eval(['y=G' num2str(ig) ';']);
   
    plot(Gx2{ig},Gy{ig});
    xlabel('k_pd');
    ylabel('Cd\_ratio');
    axis([0 2 0 2]);
    title(['Cp/U10=' num2str(xbin(ig))]);
end

 % save figure
 set(gcf,'paperunits','inches','paperposition',[0 0 12 8]);
 print(gcf,'-r300','-djpeg',[svdir filesep figname]);




