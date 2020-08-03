% This script serves to aggregate SSD Cd (for all incoming waves) from all 
% wind speed and prepare dataset in scatter format and save for final
% plotting.
% Xuanyu Chen (:o_o:)
% Aug 28, 2018
% 

%%% :o_o: global variables
global TRS_abspath LOC_abspath local_root scptdir trackdir GDrive_root
global  wspd_arry  expID tailev_arry wdir dep_selected FLD_arry
global fieldFN specFN deps_of_interest slope WDL_arry_SW


%%% :o_o: variables used locally:
%dep_selected=30;
grav=9.8;
d=dep_selected;
CpiU10_mature=1.2;   %This is not getting used now.
svfigdir=[GDrive_root filesep 'fig'];
svdatdir=[GDrive_root filesep 'mat'];
%figure plotting:
nrow=2; ncol=3;
colorn={'b';'k';'r'};

if exist(svfigdir,'dir')==0
    eval(['mkdir ' svfigdir]);
end

if exist(svdatdir,'dir')==0
    eval(['mkdir ' svdatdir]);
end

addpath([local_root filesep 'shallow_water_ssd_exp/analysis/script/matlab/IDLHur_Exps/developing_stage/scripts_for_spectral_and_FLD_data']);
%marker={'o','d','<','>','^','p'};

subjet=jet(64);
subjet=subjet([1:30 41:64],:);
Ncolor=length(subjet);
cinc=floor(Ncolor/length(deps_of_interest));
depcolor=subjet(end:-cinc:1,:);
FLDName={'URI';'Miami'};



for FLDID=FLD_arry
    FLDstr=num2str(FLDID);
    infile1=['FLD' FLDstr 'BGFDLcor_SSD_Cd_atSelectedDeps_slp' num2str(slope) '.mat'];
    svdataFN=['FLD' FLDstr 'BGFDLcor_Cd_related_data_for_all_windspeed_onSlope' num2str(slope) '.mat'];
    
    dcnt=0;
    % find dep_selected id:
    for dep_selected=deps_of_interest(1:end);
        id=find(deps_of_interest==dep_selected);
        
        if isempty(id)
            disp('data has not been extracted for this depth!');
            return
        else
            dcnt=dcnt+1;
            lgdstr{dcnt}=['d' num2str(dep_selected)];
            
            
            for itlev=1:length(tailev_arry);
                for iwspd = 1:length(wspd_arry)
                    wspdstr=num2str(wspd_arry(iwspd),'%2.2i');
                    wdirstr=num2str(wdir);
                    %%% construct path to load in data:
                    upper_loc=['wspd' wspdstr '_wdir' wdirstr filesep expID '_slp' num2str(slope)];
                    workpath=[GDrive_root filesep upper_loc ];
                    load([workpath filesep 'matdata' filesep infile1]);
                    NWDL=size(Cd_ratio_DOI{1},1);
                    
                    %%% lump data into wind speed bins:
                    Cd_ratio_box(iwspd,:)=Cd_ratio_DOI{id}(:,itlev);
                    Cd_box(iwspd,:)=Cd_DOI{id}(:,itlev);
                    Cd_dw_box(iwspd,:)=Cd_DOI{id}(:,itlev)./Cd_ratio_DOI{id}(:,itlev);
                    % mean_Cd_sw_box(iwspd)=mean(Cd_DOI{id}(:,itlev));
                    wvage_box(iwspd,:) = wvage_DOI{id}(:,itlev);
                    CpU10_box(iwspd,:) = CpU10_DOI{id}(:,itlev);
                    % CpiU10_box(iwspd,:)= CpiU10;
                    %     Cd_dw_box(iwspd,:) = Cd_dw_all(:,itlev);
                    wspd_box(iwspd,:) = repmat(wspd_arry(iwspd),1,NWDL);
                    dep_box(iwspd,:)=repmat(dep_selected, 1, NWDL);
                    % kpD_box(iwspd,:)=repmat(
                    %U10=wspd_arry(iwspd);
                    % CpiU10 .* U10
                    %                 if itlev==1
                    %                     kd_exp(iwspd,:)=grav*d ./(CpiU10 .* U10).^2;    %the second dimenson stores the variation due to differences in the incoming waves.
                    %                 end
                    
                end
                %%% now reshape data into a column vector for scatter plots.
                scatters.(lgdstr{dcnt}).Cd{itlev}(1,:)=reshape(Cd_box,1,[]);
                scatters.(lgdstr{dcnt}).Cd_ratio{itlev}(1,:)=reshape(Cd_ratio_box,1,[]);
                scatters.(lgdstr{dcnt}).wvage{itlev}(1,:)=reshape(wvage_box,1,[]);
                scatters.(lgdstr{dcnt}).CpU10{itlev}(1,:)=reshape(CpU10_box,1,[]);
                scatters.(lgdstr{dcnt}).Cd_dw{itlev}(1,:)=reshape(Cd_dw_box,1,[]);
                scatters.(lgdstr{dcnt}).wspd(1,:)=reshape(wspd_box,1,[]);
                % scatters.dep(1,:)=reshape(dep_box,1,[]);
                
                
            end
            %
            save([svdatdir filesep svdataFN],'scatters');
            
            if 1==0
                %% Now, we are able to make plots based on the above preped data:
                itlev=0;
                figure(2)
                %lgdstr={'Cd_{deep}';'1.2Cd_{deep}';'Cd_{shallow}'};
                %         jet2=jet(50);
                %         jetnew=jet2(10:50,:);
                %         jetnew=shade(jetnew,0.9);
                
                for tailev=tailev_arry
                    itlev=itlev+1;
                    
                    % specify xdata and ydata and other alias
                    xdata=wspd_scatter{itlev};
                    ydata=Cd_ratio_allwspd{itlev};
                    cdata=CpU10_allwspd{itlev};
                    
                    %  mean_Cd_dw12=1.2*mean_Cd_dw{itlev};   % 20% increase of Cd in dw.
                    
                    ip=itlev+(FLDID-1)*ncol;
                    subplot(nrow,ncol,ip);  %hsub(ip)=
                    for iftch=1:length(WDL_arry_SW)
                        indx=[1:length(wspd_arry)]+length(wspd_arry)*(iftch-1);
                        %  scatter(xdata(indx), ydata(indx), 30, cdata(indx),'filled');
                        if iftch==1
                            h(dcnt)=plot(xdata(indx),ydata(indx),'linewidth',1+0.1*iftch, ...
                                'color',shade(depcolor(id,:),0.90+0.05*(iftch-2)));
                        else
                            plot(xdata(indx),ydata(indx),'linewidth',1+0.1*iftch, ...
                                'color',shade(depcolor(id,:),0.90+0.05*(iftch-2)));
                        end
                        
                        hold on
                    end
                    
                    % colormap(jetnew);
                    plot([5 70],[1 1],'--','color',[0.5 0.5 0.5]);
                    
                    ylim([0.5 1.5]);
                    set(gca,'fontsize',12,'xtick',[5:5:70]);
                    xlim([5 70]);
                    xlabel('U_{10N} (m/s)');
                    if mod(ip,ncol)==1
                        ylabel('Cd_{sh}/Cd_{deep}');
                    end
                    
                    if ip==ncol*nrow    && id==length(deps_of_interest)
                        hlgd=legend(h,lgdstr);
                        set(hlgd,'location','southeast','fontsize',11);
                    end
                    
                    % if ip <=ncol
                    title([FLDName{FLDID} ':B=0.0' num2str(tailev,'%2.2i')],'fontsize',18);
                    %end
                    box on
                    grid on
                    % set(hsub(ip),'pos',pos{ip});
                    
                end
                
                
                
                
                
                
                pause(1.0)
                %%% save figure:
                if id==length(deps_of_interest)
                    set(gcf,'paperunit','inches','paperposition',[0 0 12 10]);
                    figname=['Slp' num2str(slope) '_Cd_allwspd_at_severalDeps_newsetup' date '.jpg'];
                    print(gcf,'-djpeg','-r350',[svfigdir filesep figname]);
                    print(gcf,'-djpeg','-r350',[figname]);
                    close all
                end
            end
        end
    end
    
end