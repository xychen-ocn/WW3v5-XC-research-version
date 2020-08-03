% This script is used to plot Misalgnment angle at different depths.

function make_IDLHur_SpaceTime_plots_for_misang(DataDir, MatFname, VarNames,unitstr, DepParms,conlev, slp)
% Description:
%
%    some figure plotting parameters need to be setup in the wrapper before enter into
%    this function.
%    some parameters are dependent on tail levels and some are not...  (the plotting script 
%    needs to be unified.)
isNom=true;
dmsz=3;   % 6x Rmax
tickinc=1;
%%------------ Section I ----------------%%
%@ 1. load data
DW=load([DataDir.DW filesep MatFname]);
SWC=load([DataDir.SW.coarse filesep MatFname]);
SWF=load([DataDir.SW.fine filesep MatFname]);
dw_thres=DepParms.dw_thres;
sw_thres=DepParms.sw_thres;
rho_a=1.225;

%@ 2. define parameters for figure plotting.
deps_q=DepParms.deps_q;

nrow=2;
ncol=length(deps_q);


nfig=ceil(length(deps_q)/4);
%pos=customize_subplot_size(nrow,ncol,0.05,0.05);

Rmax=70;
if ~isNom
    labelstr.x='km';
    labelstr.y='km';
else
    labelstr.x='X/Rmax';
    labelstr.y='Y/Rmax';

end

fontsz.cbar=12;
fontsz.axis=12;
fontsz.title=18;

if ~isNom
   % dmsz=600;
   % tickinc=200;
    
    if dmsz > 400
        domainstr='LargeDM';
    else
        domainstr='ZoomIn';
    end
else
    % normalized coord.:

   % dmsz=6;   % 6x Rmax
   % tickinc=2;
    
    if dmsz > 3
        domainstr='LargeDM';
    else
        domainstr='ZoomIn';
    end
end

pos=customize_subplot_size(nrow,ncol,0.1,0.05);
% modify subplot 1 position:
subpos1=pos{1};
plothgt=subpos1(4);
pos{1}(2)=0.5-0.5*plothgt;

switch slp
    case 2000
        rowid=1;
    case 200
        rowid=2;
end


%svfigname='Standard_WaveParms_at_different_deps_novec_2slopes';

figsvdir=[DataDir.SW.main filesep 'fig'];
if ~exist(figsvdir,'dir') 
   feval('mkdir', figsvdir);
end

% load colormap here:
load('/Users/xychen/Desktop/research_work/matlab_projects/color_schemes/homemade_colormap/cool_warm.mat');
idx=[1:10,12,14:23];
coolwarm_sub=coolwarm(idx,:);

%@ 3. generate figures
%% plot quantities independent of tail levels:

fignum=30;

for iv = 1:length(VarNames);

    varn=VarNames{iv};
    
    % contour levels of the data:
    data_conlev=conlev.(varn);
    
    % colorbar properties:
    if strcmp(varn,'wnd_mag');
        cbar_ppt.tag=['wind speed (m/s)'];
    else
        cbar_ppt.tag='\phi' ; %'\Delta \theta_{wu_*}'; %[varn ' (' unitstr{iv} ')'];
    end
    cbar_ppt.cmap=coolwarm_sub;
    cbar_ppt.cmax=max(conlev.(varn));
    cbar_ppt.cmin=min(conlev.(varn));
    %inc=round((cbar_ppt.cmax-cbar_ppt.cmin)/8);
    inc=10;
    cbar_ppt.ytick=[cbar_ppt.cmin: inc :cbar_ppt.cmax];
    
    %%% create time-space map of data:    
    for id=1: length(deps_q)
        baseid=mod(id, ncol);
        if baseid==0
            baseid=ncol;
        end
        
        ifig=ceil(id/ncol);
        
        d =deps_q(id);
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

        else
           indata=SWF.QCed_fieldVar;
           coord=SWF.coord;
           stmx=SWF.stmx;
           stmy=SWF.stmy;
           stmtime=SWF.stmtime;
           ww3_timenum=SWF.ww3_timenum;
           [tmp,dpt_maxtrix,tmp,tmp]=check_qtty_in_netCDFfile([DataDir.SW.fine filesep 'ww3.201711_fine.nc'],0,'dpt');
           
        end      
        % data deriving
        if strcmp(varn,'misang_wust')
            vec_flag=1;                % use vectors in the figures I made.
            
            udata{1}= indata.uwnd;     %indata
            vdata{1}= indata.vwnd;
            
            udata{2}=rho_a .* indata.ust_mag .* indata.uust;
            vdata{2}=rho_a .* indata.ust_mag.* indata.vust;
            
            xinc=4;
            yinc=8;
        else
            vec_flag=0;
        end
        
       ww3_time_sw=SWC.ww3_timenum;
       ww3_time_dw=DW.ww3_timenum;

      %% ------ generate the sudo-space--space 2D matrix: ------ %%
       
         if d<dw_thres
              [ny,nx,nt]=size(indata.dpt);
       %   [~,depID]=min(abs(d-indata.dpt(round(ny/2),:,round(nt/2))));
           dpt_vec=dpt_maxtrix(round(ny/2),:,1);

           
           [data_t, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(indata.(varn), ...
               coord.XX,coord.YY,d,dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);
           
           [wnddata_t, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(indata.wnd_mag, ...
               coord.XX,coord.YY,d,dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);
           
%             [udata_t{1}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(indata.ucp, ...
%                 coord.XX,coord.YY,d,dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);
%             
%             [vdata_t{1}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(indata.vcp, ...
%                 coord.XX,coord.YY,d,dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);

        if vec_flag==1
             % wind vector:
             [udata_t{1}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(udata{1}, ...
                 coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime);
             [vdata_t{1}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(vdata{1}, ...
                 coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime);
             
             % wind stress vector:
             [udata_t{2}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(udata{2}, ...
                 coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime);
             [vdata_t{2}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(vdata{2}, ...
                 coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime);
         end
         
         
         
         else
             [tmp,sid]=min(abs(ww3_time_dw-ww3_time_sw(24)));
             data_t=indata.(varn)(:,:,sid);
             wnddata_t=indata.wnd_mag(:,:,sid);
             %             udata_t{1}=indata.ucp(:,:,sid);
             %             vdata_t{1}=indata.vcp(:,:,sid);
             if vec_flag==1
                 for ii=1:2
                     udata_t{ii}=udata{ii}(:,:,sid);
                     vdata_t{ii}=vdata{ii}(:,:,sid);
                 end
             end
             coordN.XX = (coord.XX - stmx(sid)*1000)./1000;
             coordN.YY= (coord.YY - stmy(sid)*1000)./1000;
         end
         
         %% save data for plots in manuscript;
         svdata.ma{id}=data_t;
         svdata.wndmag{id}=wnddata_t;
         
         %%% modify coordinate:
      if isNom
          coordN.XX=coordN.XX./Rmax;
          coordN.YY=coordN.YY./Rmax;
      end
      svdata.coordN{id}=coordN;

        
      if 1==0
        % make plots at different depths:
        ip=baseid+ncol*(rowid-1);
        if ip~=1+ncol
            labelstr.title=['depth=' num2str(d) 'm'];
            hfig{iv}=figure(fignum+iv);
            if vec_flag==1
                generate_SpaceTime_map_at_deps(hfig{iv},ncol,ip, coordN,data_t,wnddata_t,data_conlev,conlev.wnd_mag,...
                    cbar_ppt,Rmax,labelstr,fontsz,dmsz,tickinc,pos); %, udata_t,vdata_t);
                % if add vector, also add legend of the vector.
            else
                generate_SpaceTime_map_at_deps(hfig{iv},ncol,ip, coordN,data_t,wnddata_t,data_conlev,conlev.wnd_mag,...
                    cbar_ppt,Rmax,labelstr,fontsz,dmsz,tickinc,pos)
            end
        end
      end
        
    end   % loop for different depths.
  
    if 1==0
        svfigname=[upper(varn) '_at_different_deps_withvec_2slopes'];
        % save figure:
       % for i=1:length(hfig);
       if rowid==2
            figure(hfig{iv})
            set(gcf,'paperunits','inches','paperposition',[0 0 12 6]);
            if vec_flag==0
                print(gcf,'-djpeg','-r350',[figsvdir filesep svfigname '_' domainstr '_fig' date '.jpg']);
            else
                print(gcf,'-djpeg','-r350',[figsvdir filesep svfigname '_' domainstr '_fig_vecon_' date '.jpg']);
                
            end
      %  end
       end
 
    end
end     % loop for different variables

    %% save the dataset to a location;
    svdatadir='/Users/xychen/Desktop/research_work/paper_writing/ssd_shallowater/figuresInPaper/Figure_PartII/Fig14/data';
   
    if ~exist(svdatadir,'dir')
        feval('mkdir', svdatadir);
    end
    save([svdatadir filesep 'Misang_maps_at_selected_depths_slp' num2str(slp) '.mat'],'svdata','deps_q');



return