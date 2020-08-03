function make_IDLHur_SpaceTime_plots_for_standard_wave_parameters(DataDir, MatFname, VarNames,unitstr, DepParms,conlev, slp, vecflags)
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
%rho_a=1.225;

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

pos=customize_subplot_size(nrow,ncol,0.1,0.1);
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

%@ 3. generate figures
%% plot quantities independent of tail levels:

fignum=30;

for iv = 1:length(VarNames);

    varn=VarNames{iv};
    vec_flag=vecflags(iv);
    
    % contour levels of the data:
    data_conlev=conlev.(varn);
    
    % colorbar properties:
    if strcmp(varn,'wnd_mag')
        cbar_ppt.tag=['wind speed (m/s)'];
    else
        cbar_ppt.tag=[varn ' (' unitstr{iv} ')'];
    end
    cbar_ppt.cmap='jet';
    cbar_ppt.cmax=max(conlev.(varn));
    cbar_ppt.cmin=min(conlev.(varn));
   % inc=round((cbar_ppt.cmax-cbar_ppt.cmin)/6);
  %  cbar_ppt.ytick=[cbar_ppt.cmin: inc :cbar_ppt.cmax];
    if strcmp(varn,'hs')
        inc=4;
    elseif strcmp(varn,'lm')
        inc=100;
    end
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
           indata_unmask=DW.ODATA;
           coord=DW.coord;
           stmx=DW.stmx;
           stmy=DW.stmy;
           
        elseif (d>=sw_thres)
           indata=SWC.QCed_fieldVar;
           indata_unmask=SWC.ODATA;
           coord=SWC.coord;
           stmx=SWC.stmx;
           stmy=SWC.stmy;
           stmtime=SWC.stmtime;
           ww3_timenum=SWC.ww3_timenum;
           [tmp,dpt_maxtrix,tmp,tmp]=check_qtty_in_netCDFfile([DataDir.SW.coarse filesep 'ww3.201711_coarse.nc'],0,'dpt');

        else
           indata=SWF.QCed_fieldVar;
           indata_unmask=SWF.ODATA;
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
       
         if d<dw_thres
              [ny,nx,nt]=size(indata_unmask.dpt);
       %   [~,depID]=min(abs(d-indata.dpt(round(ny/2),:,round(nt/2))));
           dpt_vec=dpt_maxtrix(round(ny/2),:,1);

           
           [data_t, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(indata_unmask.(varn), ...
               coord.XX,coord.YY,d,dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);
           
           [wnddata_t, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(indata_unmask.wnd.mag, ...
               coord.XX,coord.YY,d,dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);
           
            [udata_t{1}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(indata.ucp, ...
                coord.XX,coord.YY,d,dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);
            
            [vdata_t{1}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(indata.vcp, ...
                coord.XX,coord.YY,d,dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);
            
            
            svdata.(varn){id}=data_t;
            if iv==2
                svdata.ucp{id}=udata_t{1};
                svdata.vcp{id}=vdata_t{1};
                svdata.wndmag{id}=wnddata_t;
            end
            

            
        else  
            [tmp,sid]=min(abs(ww3_time_dw-ww3_time_sw(24)));
            data_t=indata_unmask.(varn)(:,:,sid);
            wnddata_t=indata_unmask.wnd.mag(:,:,sid);
            udata_t{1}=indata.ucp(:,:,sid);
            vdata_t{1}=indata.vcp(:,:,sid);
            coordN.XX = (coord.XX - stmx(sid)*1000)./1000;
            coordN.YY= (coord.YY - stmy(sid)*1000)./1000;
            
            svdata.(varn){id}=data_t;
             if iv==2
                svdata.ucp{id}=udata_t{1};
                svdata.vcp{id}=vdata_t{1};
                svdata.wndmag{id}=wnddata_t;
            end
            
         end
        
        %%% modify coordinate:
      if isNom
          coordN.XX=coordN.XX./Rmax;
          coordN.YY=coordN.YY./Rmax;
      end
        
           svdata.XX{id}=coordN.XX;
           svdata.YY{id}=coordN.YY;
           
        % make plots at different depths:
        if 1==0
        ip=baseid+ncol*(rowid-1);
        if ip~=1+ncol
            labelstr.title=['depth=' num2str(d) 'm'];
            hfig{iv}=figure(fignum+iv);
          %  if strcmp(varn,'hs') || strcmp(varn,'lm')
          if vec_flag==1
              generate_SpaceTime_map_at_deps(hfig{iv},ncol,ip, coordN,data_t,wnddata_t,data_conlev,conlev.wnd_mag,...
                  cbar_ppt,Rmax,labelstr,fontsz,dmsz,tickinc,pos, udata_t,vdata_t);
              % if add vector, also add legend of the vector.
          else
              generate_SpaceTime_map_at_deps(hfig{iv},ncol,ip, coordN,data_t,wnddata_t,data_conlev,conlev.wnd_mag,...
                  cbar_ppt,Rmax,labelstr,fontsz,dmsz,tickinc,pos);
          end
        end
        end
        
    end   % loop for different depths.
    
    
    
    %%
    if vec_flag==1
        svfigname=[upper(varn) '_at_different_deps_withvec_2slopes'];
    else
        svfigname=[upper(varn) '_at_different_deps_2slopes'];
    end
        % save figure:
       % for i=1:length(hfig);
       if 1==0
       if rowid==2
            figure(hfig{iv})
            set(gcf,'paperunits','inches','paperposition',[0 0 12 6]);
            print(gcf,'-djpeg','-r350',[figsvdir filesep svfigname '_' domainstr '_' date '.jpg']);
            savefig(gcf,[figsvdir filesep svfigname '_' domainstr '_' date '.fig']);
      %  end
       end
   
       end

end     % loop for different variables

    %% save the dataset to a location;
    svdatadir='/Users/xychen/Desktop/research_work/paper_writing/ssd_shallowater/figuresInPaper/Figure_PartII/Fig02_03/data';
    save([svdatadir filesep 'Hs_LM_maps_at_selected_depths_slp' num2str(slp) '.mat'],'svdata','deps_q');



return