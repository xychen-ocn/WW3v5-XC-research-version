function make_IDLHur_scatter_plots(DataDir, MatFname,VarNames, slp, DepParms, isScatVis, isRefData, ...
                           Ref_xdata, Ref_ydata, RefData_ppt,Ylabelstr, YMin, YMax)
% Descriptions:
%     RefData_ppt: Reference data properties. (This needs to be specified again in the wrapper function. 
%                  (not yet accomplished.)
%     Function needs to be built to. (get_IDLHur_QCed_data_at_depths).
%
%   This script will output (generate) 
%         a data set with statistics computed on the variability of requested variables in shallow water.
%         figures showing the variable in scatters at different depth and with different spectral tail levels.
global mwstr
%-----------------------Data I/O -------------------------------%
% Need to load in data to begin with..
DW=load([DataDir.DW filesep MatFname]);
SWC=load([DataDir.SW.coarse filesep MatFname]);
SWF=load([DataDir.SW.fine filesep MatFname]);
dw_thres=DepParms.dw_thres;
sw_thres=DepParms.sw_thres;

OutMatFname='Cd_variability_of_variables_at_several_deps.mat';


%% ---------------- Section I -------------------- %%
%@@ 0. define parameters used for figure plotting:
xvarn='wnd_mag';
binwidth_xdata=5;
binrange_xdata=[10,str2num(mwstr)];

% color code data: can I compute it?
cc_flag=1;
%cc_varn='inpwvage';                            % color code quantity, I can also consider Cp/U10 as representation of wave age.
cc_varn_labelstr= 'Wave Age (c_{p}cos\theta/u_{*})';
ccdata_range=[-30 30];  %[-3 0.5]; %
cinc=10;
colorTicks=[ccdata_range(1):cinc:ccdata_range(2)];
num_color=length(colorTicks)-1;
sjet=shade(jet(64),0.9);
cmap=sjet(10:2:10+2*(num_color*4-1),:);

%cbar_font=14;
Xlabelstr='U_{10} (m/s)';
ScaSize=10;

fontsz.cbar=12;
fontsz.axis=12;
fontsz.title=18;

Xrange=[0 70];
ytickinc=1;
YScaleStr='linear';


if isScatVis
    extstr='scatter_on';
else
    extstr='scatter_off';
end

%@@ 1. specify the dimension of subplots:
deps_q=DepParms.deps_q;
nrow=2;              % number of tail levels.
ncol=4;
nfig=ceil(length(deps_q)/4);
% 
pos=customize_subplot_size(nrow,ncol,0.08,0.08);
% modify subplot 1 position:
subpos1=pos{1};
plothgt=subpos1(4);
pos{1}(2)=0.5-0.5*plothgt;


% switch B
%    case 2
%        rowid=1;
%    case 6
%        rowid=2;
%    case 12
%        rowid=3;
% end

switch slp
    case 2000
    %case 45
        rowid=1;
    case 200
   % case 135
        rowid=2;
end


%@@ 2. create folder to save figures:
figsvdir=[DataDir.SW.main filesep 'fig'];
if ~exist(figsvdir,'dir')
    feval('mkdir', figsvdir);
end

GFDL=load('/Users/xychen/Desktop/research_work/matlab_projects/Drag_Formula/GFDL_Cd.mat');

%% ---------------- Section II -------------------- %%
%@@ 1. get necessary data at request depth first:
for id=1:length(deps_q);
    
    d =deps_q(id);
    
    if (d>=dw_thres) 
       indata=DW.DataSubVec;
    elseif (d>=sw_thres)
       indata=SWC.DataSubVec;
    else
       indata=SWF.DataSubVec;
    end
       
    % getting the variables for y-axis:
    for iv=1:length(VarNames);
        varn=VarNames{iv};
        [ydata_tmp]=sort_data_by_depth(indata.(varn),indata.dpt, ...
                                      d, DepParms.dep_binwidth);
        ydata(id).(varn)=ydata_tmp{1};
    end
    
    % variables for x-axis (wind speed)
    xdata{id}=sort_data_by_depth(indata.(xvarn), indata.dpt, ...
                                 d, DepParms.dep_binwidth);
    
    % variables for color-code:
    wvage=indata.cp.* cosd(indata.misang_ustwv) ./indata.ust_mag;
    ccdata{id}=sort_data_by_depth(wvage, indata.dpt, ...
                                 d, DepParms.dep_binwidth);     %indata.(cc_varn)
   
    % 
    relx=sort_data_by_depth(indata.XX, indata.dpt, ...
                                 d, DepParms.dep_binwidth);
                             
    rely=sort_data_by_depth(indata.YY, indata.dpt, ...
                                 d, DepParms.dep_binwidth);
                             
    %% further sorting: 
    % sort Cd data by quadrants:
    RIDs=find(rely{1}>=0);
    LIDs=find(rely{1}<0);
    
   RHS.ydata(id).(varn)=ydata(id).(varn)(RIDs);
    LHS.ydata(id).(varn)=ydata(id).(varn)(LIDs);
    RHS.xdata{id}=xdata{id}{1}(RIDs);
    LHS.xdata{id}=xdata{id}{1}(LIDs);
    RHS.ccdata{id}=ccdata{id}{1}(RIDs);
    LHS.ccdata{id}=ccdata{id}{1}(LIDs);
    
   
    
end
     % save the data for each typhoon cases, and plot it separately. 
     save([DataDir.SW.main filesep 'Cd_split_LeftRight.mat'],'RHS','LHS'); 
     
%@@ 2. Data binning and mean/median computing:
for iv=1:length(VarNames)
    varn=VarNames{iv};
for id=1:length(deps_q)
    %%% --> a. bining data in different quadrants:
    % do right hand side for now:
    xdata_tmp=RHS.xdata{id};
    y=RHS.ydata(id).(varn);
    [ybinned, xbinned, xbin_center, xidx_binned]=sort_ydata_by_xdata(y,xdata_tmp, ...
                                                           binwidth_xdata,binrange_xdata);
    ydata_binned(id).(varn)=ybinned;
    ydata_mean(id).(varn)=mean(ybinned,'omitnan');
    ydata_median(id).(varn)=median(ybinned,'omitnan');
    
    xdata_binned{id}=xbinned;
    xbin_cen{id}=xbin_center;
  
    %%% --> b. compute statistics to represent variability from the binned ydata:
    ydata_stat(id).(varn)=ssdexp_compute_statistics(ybinned);  
    


        
end
end
clear xdata_tmp xidx_binned
%%% --> save data:
% filename=[DataDir.SW.main filesep OutMatFname];
% save(filename, 'ydata_stat','xbin_cen', 'DepParms');
%save([DataDir.SW.main filesep 'JPDF_for_misang_scattered_data.mat'],'PDF_data','PDF_val','JPDFparm');
% save([DataDir.SW.main filesep 'JPDF_for_Cd_scattered_data.mat'],'PDF_data','PDF_val','JPDFparm');

% h =  findobj('type','figure');
% fignum = length(h);
% fignum=fignum+1;

if 1==1
%% ---------------- Section III -------------------- %%
%@@  now, ready to make plots:
for iv=1:length(VarNames);
    varn=VarNames{iv};
    for id=1: length(deps_q);
        baseid=mod(id,ncol);
        if baseid ==0
            baseid=ncol;
        end
       % ifig=ceil(id/ncol);
        ifig=1;
        % index for sub-panel:
        ip=baseid+(rowid-1)*ncol;
        
        hfig{ifig}=figure(iv+ifig);
        
        set(hfig{ifig},'name',[varn 'scatters']);
        
        
        %%% make subplots:
        if ip~=1+ncol
        hsb=subplot(nrow,ncol,ip);
        xdata_tmp=RHS.xdata{id};
        ccdata_tmp=RHS.ccdata{id};
%         ccdata_tmp=log10(PDF_data(id).(varn));
        
        %% add polynominal fit:
        %         [p,S]=polyfit(xdata_tmp{1}, ydata(id).(varn), 5);
        %         yfit=polyval(p,xbin_cen{id});
        if cc_flag==0
            hsc=scatter(hsb,xdata_tmp,RHS.ydata(id).(varn), ScaSize,[0.68 0.68 0.68]);
        else
            hsc=scatter(hsb,xdata_tmp,RHS.ydata(id).(varn), ScaSize, ccdata_tmp);
            colormap(cmap);
            hb=colorbar;
            hdl.colorbar=hb;
        end
        hold on
   
        hdl.scatter=hsc;
%        
        % reference median Cd in deep water.
        if strcmp(varn,'Cd')
            hl(1)=plot(GFDL.U10,GFDL.Cd*10^(3), '--b','linewidth',1.2);
            lgdstr{1}='GFDL';
        end

        % plot statistic box:
       bplot(ydata_binned(id).(varn),'nooutliers','whisker',2.5,'width',2.8,'boxloc',xbin_cen{id});

%         hdl=scatter_boxplot(hsb,xdata_tmp{1},ydata(id).(varn),ccdata_tmp{1}, ScaSize, ...
%                             ydata_binned(id).(varn),xbin_cen{id}, isScatVis);
%         hold on;
%        
        clear xdata_tmp ccdata_tmp
       

        
%            hl(1)=plot(xbin_cen{1},ydata_mean(1).(varn),'-','color',[0.55 0.55 0.55],'linewidth',1.2);
%            lgdstr{1}='DW mean';
%            
%         if id~=1
%             hl(2)=plot(xbin_cen{id},ydata_mean(id).(varn),'--','color','r','linewidth',1.2);
%             lgdstr{2}='SW mean';
%         end
        
        
        % add additional reference data:
        if isRefData   
            n=length(hl);
            
            href=plot(Ref_xdata.(varn),Ref_ydata.(varn));

            %%% set line properties:
         for i=1:length(href)
             %href(i).Color=RefData_color{idx(i)};
             %href(i).LineStyle=RefData_lnsty{idx(i)};
             %href(i).LineWidth=RefData_lw(idx(i));
            % href(i).Marker=RefData_marker{idx(i)};
            href(i).Color=RefData_ppt.color{i}.(varn);
            href(i).LineStyle=RefData_ppt.lnsty{i}.(varn);
            href(i).LineWidth=RefData_ppt.lw(i).(varn);
            href(i).Marker=RefData_ppt.marker{i}.(varn);
            lgdstr{n+i}=RefData_ppt.tag{i}.(varn);
         end
         
         hll=[hl href'];
         
        end
        
        
        % add legend:
        
        if cc_flag==1
            % add add colorbar and labels:
            caxis(hsb,ccdata_range)
            hdl.colorbar.Ticks=colorTicks;
            %set(get(hdl.colorbar,'title'),'String','10^{\^}','units','data','pos',[0.5 0.51 0]);
            
            if mod(ip,ncol) ==0
                hdl.colorbar.XLabel.String=cc_varn_labelstr;
                hdl.colorbar.XLabel.FontSize=fontsz.cbar;
            else
                hdl.colorbar.Visible='off';
                
            end
        end

        % touch on ylabel:
        %if mod(ip,ncol)==1
            ylabel(Ylabelstr.(varn));
       % end

        % touch on xlabel:
        %if rowid==2   %3
            xlabel(Xlabelstr);           
        %end
        set(gca,'fontsize',fontsz.axis);
        
        % touch on title
        %if rowid==1 %1
            title([num2str(deps_q(id)) 'm'],'fontsize',fontsz.title);
        %end

        ylim([YMin.(varn) YMax.(varn)]);             
        xlim(Xrange);
        set(gca,'xtick',[Xrange(1):10:Xrange(2)],'yscale',YScaleStr);
        set(gca,'ytick',[YMin.(varn):ytickinc:YMax.(varn)]);
        axis('square')
        grid on
        
        %         if ip==1
        %             hlgd=legend(hll,lgdstr);
        %             set(hlgd,'loc','northoutside','fontsize',11);
        %             %pause
        %         end
        if ip< 1+ncol
            icnt=ip;
        else
            icnt=ip-1;
        end
       label_panels(hsb,'northwest',14,icnt);

        % reposition subplot afterwards.
        set(hsb,'pos',pos{ip},'box','off');
        end
       % pause
        
    end      % end of depth loop;    
    if rowid==2
        for i=1:length(hfig)
            figure(hfig{i});
            if cc_flag==0
                figname=[varn '_scatters_at_selected_depth_' extstr '_RHS' date];
            else
                figname=[varn '_ColorCode_scatters_at_selected_depth_' extstr '_RHS' date];
            end
            set(gcf,'paperunits','inches','paperposition',[0 0 12 6]);
            print(gcf,'-djpeg','-r350',[figsvdir filesep figname '.jpg']);
            savefig([figsvdir filesep figname '.fig']);
            close gcf
        end
    end
end
end

                           
%%%% Obsolete:        
% plot JPDF again.
% contour(JPDFparm(id).Xcen, JPDFparm(id).Ycen, log10(PDF_val(id).(varn)'),[-10:0.01:1]);
                           
                           

return