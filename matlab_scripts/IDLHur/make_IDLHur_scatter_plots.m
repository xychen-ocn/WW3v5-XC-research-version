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
    
    %% save the data obtain above into a matrix:
    scadata.U10{id}=xdata{id}{1};
    scadata.wvage{id}=ccdata{id}{1};
    scadata.Cd{id}=ydata(id).(varn);

end


    
%@@ 2. Data binning and mean/median computing:
for iv=1:length(VarNames);
    varn=VarNames{iv};
    for id=1:length(deps_q);
        %%% --> a. bin ydata by xdata:
        xdata_tmp=xdata{id};
        [ybinned, xbinned, xbin_center, xidx_binned]=sort_ydata_by_xdata(ydata(id).(varn),xdata_tmp{1}, ...
            binwidth_xdata,binrange_xdata);
        ydata_binned(id).(varn)=ybinned;
        ydata_mean(id).(varn)=mean(ybinned,'omitnan');
        ydata_median(id).(varn)=median(ybinned,'omitnan');
        
        xdata_binned{id}=xbinned;
        xbin_cen{id}=xbin_center;
        
        %%% --> b. compute statistics to represent variability from the binned ydata:
        ydata_stat(id).(varn)=ssdexp_compute_statistics(ybinned);
        
        %%% --> compute JPDF for data at different depths:
        % parm. for JPDF:
        %         JPDFparm(id).Xedges=[10:0.1:ceil(max(xdata_tmp{1}))];
        % %         inc=(ceil(max(ydata(id).(varn)))-floor(min(ydata(id).(varn))))*0.05;
        % %         JPDFparm(id).Yedges=[floor(min(ydata(id).(varn))):inc:ceil(max(ydata(id).(varn)))];
        %         inc=(ceil(max(ydata(id).(varn)))-floor(min(ydata(id).(varn))))*0.01;
        %         JPDFparm(id).Yedges=[floor(min(ydata(id).(varn))):inc:ceil(max(ydata(id).(varn)))];
        %  %        JPDFparm(id).Yedges=[0.5:0.01:4.5];
        %
        %
        %         JPDFparm(id).Xcen=0.5*(JPDFparm(id).Xedges(1:end-1)+JPDFparm(id).Xedges(2:end));
        %         JPDFparm(id).Ycen=0.5*(JPDFparm(id).Yedges(1:end-1)+JPDFparm(id).Yedges(2:end));
        
        
        figure(10);
        %     hjpdf=histogram2(gca,xdata_tmp{1},ydata(id).(varn),  ...
        %                JPDFparm(id).Xedges, JPDFparm(id).Yedges,'Normalization','pdf');
        %     pdf_val=hjpdf.Values;
        
        %%% call a function to colorcode the Cd:
        %     valid_x=xdata_tmp{1}(isnan(xdata_tmp{1})==0);
        %     valid_y=ydata(id).(varn)(isnan(xdata_tmp{1})==0);
        %     PDF_data(id).(varn)=find_pdf_for_each_point(valid_x, valid_y, JPDFparm(id).Xedges, JPDFparm(id).Yedges, pdf_val);
        %     PDF_val(id).(varn)=pdf_val;
        
        
        scadata.Cd_binned{id}=ybinned;
        scadata.U10_cen{id}=xbin_center;
        
    end
    

end



clear xdata_tmp xidx_binned
%%% --> save data:
filename=[DataDir.SW.main filesep OutMatFname];
save(filename, 'ydata_stat','xbin_cen', 'DepParms');


%% save data for plots in the manuscript:
svdatadir='/Users/xychen/Desktop/research_work/paper_writing/ssd_shallowater/figuresInPaper/Figure_PartII/Fig10/data';

if ~exist(svdatadir,'dir')
    feval('mkdir', svdatadir);
end
save([svdatadir filesep 'Cd_scatters_at_selected_depths_slp' num2str(slp) '.mat'],'scadata','deps_q');




%save([DataDir.SW.main filesep 'JPDF_for_misang_scattered_data.mat'],'PDF_data','PDF_val','JPDFparm');
% save([DataDir.SW.main filesep 'JPDF_for_Cd_scattered_data.mat'],'PDF_data','PDF_val','JPDFparm');





if 1==0
%% ---------------- Section III -------------------- %%
%@@  now, ready to make plots:
for iv=1:length(VarNames)
    varn=VarNames{iv};
    
    for id=1: length(deps_q)
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
        xdata_tmp=xdata{id};
        ccdata_tmp=ccdata{id};
%         ccdata_tmp=log10(PDF_data(id).(varn));
        
        %% add polynominal fit:
        %         [p,S]=polyfit(xdata_tmp{1}, ydata(id).(varn), 5);
        %         yfit=polyval(p,xbin_cen{id});
        if cc_flag==0
            hsc=scatter(hsb,xdata_tmp{1},ydata(id).(varn), ScaSize,[0.68 0.68 0.68]);
        else
            hsc=scatter(hsb,xdata_tmp{1},ydata(id).(varn), ScaSize, ccdata_tmp{1});
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
            lgdstr{1}='GFDL Cd';
        end

        % plot statistic box:
       bplot(ydata_binned(id).(varn),'nooutliers','whisker',2.5,'width',2.8,'boxloc',xbin_cen{id});
       
       %% add 1.5x IQR: not necessary --> But it is a good exercise for myself. 
       % add 1.5x IQR above 3rd and below 1st quartile.
       IQR=ydata_stat(id).(varn).IQR;
       boxedge=ydata_stat(id).(varn).boxEdge;
       
       ythres_up=1.5.*IQR + boxedge(2,:);
       ythres_low=-1.5.*IQR + boxedge(1,:);
       % upper outliers:
       hl(2)=plot(xbin_cen{id}, ythres_up,'linestyle','-.','linewidth',1.2,'color',[0.5 0.5 0.5]);
       % lower outliers:
       plot(xbin_cen{id}, ythres_low,'linestyle','-.','linewidth',1.2,'color',[0.5 0.5 0.5]);
       
%         hdl=scatter_boxplot(hsb,xdata_tmp{1},ydata(id).(varn),ccdata_tmp{1}, ScaSize, ...
%                             ydata_binned(id).(varn),xbin_cen{id}, isScatVis);
%         hold on;
%        
        clear xdata_tmp ccdata_tmp
       %%

        
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
        else
            hll=hl;
            lgdstr{2}='outlier threshold';
         
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
        
        if ip==1
            hlgd=legend(hll,lgdstr);
            set(hlgd,'loc','southoutside','fontsize',11);
             lgdpos=get(hlgd,'pos');
             set(hlgd,'loc','northoutside','fontsize',11);
            %pause
        end
        if ip< 1+ncol
            icnt=ip;
        else
            icnt=ip-1;
        end
       label_panels(hsb,'northwest',14,icnt);

        % reposition subplot afterwards.
        set(hsb,'pos',pos{ip},'box','off');
        
        
        if 1==0
        %% build an extra legend for the box plot:
        %%% make legend plot in another axes
        if ip==1
            a1=gca;
           
            a2pos=[lgdpos(1)+0.28*lgdpos(3), lgdpos(2)*0.15, lgdpos(3), lgdpos(4)*2];
            a2=axes('pos',a2pos);
            
            % make distribution:
           rpd= makedist('Rayleigh',0.5);
           randata=random(rpd,200,1);
           mean(randata)
           median(randata)
           

            [T, boxEdge, wisEdge]=bplot(randata,'nooutliers','whisker',2.5,'width',2.);
            
            text(2.5,mean(randata)*1.1,'mean','fontsize',10,'color','r');
            text(2.5,median(randata)*0.95,'median','fontsize',10,'color','k');
            text([2.5 2.5], wisEdge, {'2.5%';'97.5%'});
            text([0.5-2 0.5-2], boxEdge, {'25%';'75%'},'horizontalalignment','center');
            
            a2.XAxis.Visible='on';
            a2.YAxis.Visible='on';
            a2.XTick='';
            a2.YTick='';
            a2.XTickLabel='';
            a2.YTickLabel='';
            a2.Box='on';
            axis('square');
            
            pause
            
            
        end
        end
        end
       % pause
        
    end      % end of depth loop;    
    if rowid==2
        for i=1:length(hfig)
            figure(hfig{i});
            if cc_flag==0
                figname=[varn '_scatters_at_selected_depth_' extstr '_updated' date];
            else
                figname=[varn '_ColorCode_scatters_at_selected_depth_' extstr '_updated' date];
            end
            set(gcf,'paperunits','inches','paperposition',[0 0 12 6]);
            print(gcf,'-djpeg','-r350',[figsvdir filesep figname '_fig' num2str(i) '_' date '.jpg']);
            % savefig([figsvdir filesep figname '.fig']);
            %close gcf
        end
    end
end
end

                           
%%%% Obsolete:        
% plot JPDF again.
% contour(JPDFparm(id).Xcen, JPDFparm(id).Ycen, log10(PDF_val(id).(varn)'),[-10:0.01:1]);
                           
                           

return