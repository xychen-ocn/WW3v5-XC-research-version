% This script will create a plotting script 

function make_IDLHur_SpaceTime_plots_for_Cd(DataDir, MatFname, VarNames, unitstr, slp, DepParms,conlev)
% Description:
%
%    some figure plotting parameters need to be setup in the wrapper before enter into
%    this function.
%    some parameters are dependent on tail levels and some are not...  (the plotting script 
%    needs to be unified.)
isNom=true;       % false
dmsz=3;           % 300
 

%%------------ Section I ----------------%%
%@ 1. load data
DW=load([DataDir.DW filesep MatFname]);
SWC=load([DataDir.SW.coarse filesep MatFname]);
SWF=load([DataDir.SW.fine filesep MatFname]);
load([DataDir.SW.main filesep 'CdRatio_atDifferentDeps_to_DW.mat']);

dw_thres=DepParms.dw_thres;
sw_thres=DepParms.sw_thres;
rho_a=1.225;
meanCdFname='DW_Cd.mat';
bulkchoice='mean';              %or median.

%@ 2. define parameters for figure plotting.

deps_q=DepParms.deps_q;
nrow=2;
ncol=length(deps_q);
nfig=ceil(length(deps_q)/4);
pos=customize_subplot_size(nrow,ncol,0.1,0.05);
% modify subplot 1 position:
subpos1=pos{1};
plothgt=subpos1(4);
pos{1}(2)=0.5-0.5*plothgt;

Rmax=70;
if ~isNom
    labelstr.x='km';
    labelstr.y='km';
else
    labelstr.x='X/Rmax'; %'Normalized Dist.';  
    labelstr.y='Y/Rmax'; %'Normalized Dist.';  

end

fontsz.cbar=12;
fontsz.axis=12;
fontsz.title=18;

% absolute coord.
if ~isNom
  %  dmsz=300;
    tickinc=100;
    
    if dmsz > 400
        domainstr='LargeDM';
    else
        domainstr='ZoomIn';
    end
else
    % normalized coord.:
   % dmsz=3;   % 6x Rmax
    tickinc=1;
    
    if dmsz > 3
        domainstr='LargeDM_normalized';
    else
        domainstr='ZoomIn_normalized';
    end
end



switch slp
   case 2000
   % case 45
       rowid=1;
   case 200
  %  case 135
       rowid=2;
end


figsvdir=[DataDir.SW.main filesep 'fig'];

if ~exist(figsvdir,'dir') 
   feval('mkdir', figsvdir);
end


% load and/or set colormap here:
colormap_repo='/Users/xychen/Desktop/research_work/matlab_projects/color_schemes/homemade_colormap';
load([colormap_repo filesep 'cool_warm.mat']);
idx=[6:10,12,15:19];
coolwarm_sub=coolwarm(idx,:);

%load([colormap_repo filesep 'newjet.mat']);

%jet_sub=jet(64);
% tuning the colormap; (really tedious..)
jet_sub=jet(64);
jet_sub=jet_sub(11:2:58,:);


%jet_sub=jet(length(conlev.Cd)*2);



%@ 4. generate figures
%% plot quantities dependent of tail levels:
% h =  findobj('type','figure');
% fignum = length(h);
% 
%  fignum=fignum+1;
for iv = 1:length(VarNames);
    varn=VarNames{iv};

    
    % contour levels:
       data_conlev=conlev.(varn);
       
    % colorbar properties:  
        cbar_ppt.tag=[varn ' (' unitstr{iv} ')'];
        cbar_ppt.cmap=jet_sub;
        cbar_ppt.cmax=3.5;   %max(conlev.(varn));
        cbar_ppt.cmin=min(conlev.(varn));
        cbar_ppt.ytick=[1:0.5:3.5];
       
    % data deriving
    if strcmp(varn,'misang_wust')
        vec_flag=1;                % use vectors in the figures I made.
        
        udata(1)= indata.uwnd;     %indata
        vdata(1)= indata.vwnd;
        
        udata(2)=rho_a .* indata.ust_mag .* indata.uust;
        vdata(2)=rho_a .* indata.ust_mag.* indata.vust;
        
        xinc=4;
        yinc=8;
    else
        vec_flag=0;
    end

    for id = 1:length(deps_q)
        baseid=mod(id,ncol);
        if baseid ==0
            baseid=ncol;
        end
        ifig=ceil(id/ncol);
        
        d = deps_q(id);
        disp(['depth=' num2str(d)]);
        
        if (d>=dw_thres) 
           indata=DW.QCed_fieldVar;
           coord=DW.coord;
           stmx=DW.stmx;
           stmy=DW.stmy;
           
        elseif (d>=sw_thres)
           indata=SWC.QCed_fieldVar;
           coord=SWC.coord;
           stmx=SWC.stmx;
           stmy=SWC.stmy;
           stmtime=SWC.stmtime;
           ww3_timenum=SWC.ww3_timenum;
          [tmp,dpt_maxtrix,tmp,tmp]=check_qtty_in_netCDFfile([DataDir.SW.coarse filesep 'ww3.201711_coarse.nc'],0,'dpt');

           
        elseif (d<sw_thres)
           indata=SWF.QCed_fieldVar;
           coord=SWF.coord;
           stmx=SWF.stmx;
           stmy=SWF.stmy;
           stmtime=SWF.stmtime;
           ww3_timenum=SWF.ww3_timenum;
           [tmp,dpt_maxtrix,tmp,tmp]=check_qtty_in_netCDFfile([DataDir.SW.fine filesep 'ww3.201711_fine.nc'],0,'dpt');

           
        end      
       
       ww3_time_sw=SWC.ww3_timenum;
       ww3_time_dw=DW.ww3_timenum;
        
        %% ------ generate the sudo-space--space 2D matrix: ------ %%
        %%% where to add the Cd ratio dataset?
        depn=['D' num2str(d)];
        data_t=Cdts.(depn).value;
        
        Cd2bulk_ratio=Cdts.(depn).ratio2bulk;
        if d>=dw_thres
            Cd2deep_ratio=data_t;
        else
            Cd2deep_ratio=Cdts.(depn).ratio2deepSSD;
        end
        
        coordN=Cdts.(depn).coord;
        
        % Cd itself:
        svdata.Cd{id}=data_t;
        svdata.Cdratio2dw{id}=Cd2deep_ratio;
        svdata.Cdratio2blk{id}=Cd2bulk_ratio;
        
               
        
        if d< dw_thres
            [ny,nx,nt]=size(indata.dpt);
            dpt_vec=dpt_maxtrix(round(ny/2),:,1);
            
            %  find the wind data:
            [wnddata_t, TT,YT,XX_sudo,coordN0]=convert_to_time_and_space_coordinate_updated(indata.wnd_mag, ...
                coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);
            svdata.wndmag{id}=wnddata_t;
        else
           %[tmp,sid]=min(abs(ww3_time_dw-ww3_time_sw(end)));
            wnddata_t=Cdts.(depn).wnd; %indata.wnd_mag(:,:,sid);
            svdata.wndmag{id}=wnddata_t;
        end

      if vec_flag==1
        % wind vector:
        [udata_t{1}, TT,YT,XX_sudo,coordN0]=convert_to_time_and_space_coordinate_updated(udata(1), ...
            coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime);            
        [vdata_t{1}, TT,YT,XX_sudo,coordN0]=convert_to_time_and_space_coordinate_updated(vdata(1), ...
            coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime);
            
        % wind stress vector:
        [udata_t{2}, TT,YT,XX_sudo,coordN0]=convert_to_time_and_space_coordinate_updated(udata(2), ...
            coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime);           
        [vdata_t{2}, TT,YT,XX_sudo,coordN0]=convert_to_time_and_space_coordinate_updated(vdata(2), ...
            coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime); 
      end
        
%%% modify coordinate:
if isNom
    coordN.XX=coordN.XX./Rmax;
    coordN.YY=coordN.YY./Rmax;
end
svdata.coordN{id}=coordN;
%%% create time-space map of data:

% figure(1)
% pcolor(coordN.XX, coordN.YY, data_t); shading interp;
% pause

labelstr.title=['Depth=' num2str(d) 'm'];
  ip=baseid+ncol*(rowid-1);
%%% make plot:
if 1==1
        hfig{1}=figure(iv+length(VarNames)+1);
        set(hfig{1},'name',varn);
        
      

        if ip~=1+ncol*(rowid-1) || ip==1
            if vec_flag == 1
                generate_SpaceTime_map_at_deps(hfig{ifig},ncol,ip, coordN,data_t,wnddata_t,data_conlev,conlev.wnd_mag,...
                    cbar_ppt,Rmax,labelstr,fontsz,dmsz,tickinc,pos,udata_t,vdata_t)
            else
                generate_SpaceTime_map_at_deps(hfig{ifig},ncol,ip, coordN,data_t,wnddata_t,data_conlev,conlev.wnd_mag,...
                    cbar_ppt,Rmax,labelstr,fontsz,dmsz,tickinc,pos)
            end
            
        end

        
        %% plot Cd ratio:
        if 1==1
           cbar_ppt2.cmap=coolwarm_sub;      %'redblue';
           cbar_ppt2.cmax=1.5+1/22; %max(conlev.Cd_ratio);
           cbar_ppt2.cmin=0.5-1/22; %min(conlev.Cd_ratio);
           cbar_ppt2.tag='Cd^{bulk}/Cd^{bulk}_{o}';
           cbar_ppt2.ytick=[0.5:0.1:1.5];
           
           
           hexf{1,1}=figure(iv*10+1);
           set(hexf{1,1},'name','Cd2bulk_ratio');

           if ip~=1+ncol*(rowid-1) ||ip==1
               generate_SpaceTime_map_at_deps(hexf{1,1}, ncol,ip, coordN,Cd2bulk_ratio,wnddata_t,conlev.Cd_ratio, conlev.wnd_mag,...
                   cbar_ppt2,Rmax,labelstr,fontsz,dmsz,tickinc,pos)
%            elseif ip==1
%                
%                generate_SpaceTime_map_at_deps(hfig{ifig},ncol,ip, coordN,data_t,wnddata_t,data_conlev,conlev.wnd_mag,...
%                    cbar_ppt,Rmax,labelstr,fontsz,dmsz,tickinc,pos)
%                
               
           end
           % end
           
          cbar_ppt2.tag='Cd_{sh}/Cd_{deep}';
           hexf{2,1}=figure(iv*20+1);   
           set(hexf{2,1},'name','Cd2deep_ratio');
           if ip~=1+ncol*(rowid-1) 
               generate_SpaceTime_map_at_deps(hexf{2,1}, ncol,ip, coordN,Cd2deep_ratio,wnddata_t,conlev.Cd_ratio,conlev.wnd_mag,...
                   cbar_ppt2,Rmax,labelstr,fontsz,dmsz,tickinc,pos)
           elseif ip==1
               
               generate_SpaceTime_map_at_deps(hexf{2,1},ncol,ip, coordN,data_t,wnddata_t,data_conlev,conlev.wnd_mag,...
                   cbar_ppt,Rmax,labelstr,fontsz,dmsz,tickinc,pos)
               
               
           end
        end
end
    end% depth loop
    
    if 1==1
    if rowid==2
        % save figure:
        for i=1:length(hfig);
            svfigname=[varn '_at_different_deps_with_desgnBGFDL'];
            figure(hfig{i});
            set(gcf,'paperunits','inches','paperposition',[0 0 12 6]);
            print(gcf,'-djpeg','-r350',[figsvdir filesep svfigname '_' domainstr '_fig' num2str(i) '_otherDeps' date '_v2.jpg']);
           % close gcf
        end
        
        if 1==1
        for i=1:size(hexf,2)
            svfigname='Cd2bulk_ratio_at_different_deps_with_desgnBGFDL';
            figure(hexf{1,i});
            set(gcf,'paperunits','inches','paperposition',[0 0 12 6]);
            print(gcf,'-djpeg','-r350',[figsvdir filesep svfigname '_' domainstr '_fig' num2str(i) '_otherDeps' date '.jpg']);
%            % close gcf
   
            svfigname='Cd2deep_ratio_at_different_deps_with_desgnBGFDL';
            figure(hexf{2,i});
            set(gcf,'paperunits','inches','paperposition',[0 0 12 6]);
            print(gcf,'-djpeg','-r350',[figsvdir filesep svfigname '_' domainstr '_fig' num2str(i) '_otherDeps' date '.jpg']);
         %  close gcf
       end
        end
        
    end
    end
    
  
  
end        % variable loop
    %% save the dataset to a location;
    svdatadir='/Users/xychen/Desktop/research_work/paper_writing/ssd_shallowater/figuresInPaper/Figure_PartII/Fig06_07/data';
   
    if ~exist(svdatadir,'dir')
        feval('mkdir', svdatadir);
    end
    save([svdatadir filesep 'Cd_and_ratios_maps_at_selected_depths_slp' num2str(slp) '_extra.mat'],'svdata','deps_q');


return