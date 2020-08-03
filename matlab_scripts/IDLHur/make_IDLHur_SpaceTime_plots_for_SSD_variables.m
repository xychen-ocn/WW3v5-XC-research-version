function make_IDLHur_SpaceTime_plots_for_SSD_variables(DataDir, MatFname, VarNames, unitstr, B, DepParms,conlev)
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
dw_thres=DepParms.dw_thres;
sw_thres=DepParms.sw_thres;
rho_a=1.225;
meanCdFname='DW_Cd.mat';
bulkchoice='mean';              %or median.

%@ 2. define parameters for figure plotting.
nrow=3;
ncol=4;
deps_q=DepParms.deps_q;
nfig=ceil(length(deps_q)/4);
pos=customize_subplot_size(nrow,ncol,0.1,0.05);
jet_sub=jet(64);
jet_sub=jet_sub(17:2:52,:);

%jet_sub=jet(length(conlev.Cd)*2);

Rmax=70;
if ~isNom
    labelstr.x='km';
    labelstr.y='km';
else
    labelstr.x='Normalized Dist.';
    labelstr.y='Normalized Dist.';

end

fontsz.cbar=12;
fontsz.axis=11;
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



switch B
   case 2
       rowid=1;
   case 6
       rowid=2;
   case 12
       rowid=3;
end


figsvdir=[DataDir.SW.main filesep 'fig'];

if ~exist(figsvdir,'dir') 
   feval('mkdir', figsvdir);
end


% load colormap here:
load('/Users/xychen/Desktop/research_work/matlab_projects/color_schemes/homemade_colormap/cool_warm.mat');
idx=[1:10,12,14:23];
coolwarm_sub=coolwarm(idx,:);

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
      if strcmp(varn,'fp')
         cbar_ppt.tag=[ 'Tp (' unitstr{iv} ')'];
      else
        cbar_ppt.tag=[varn ' (' unitstr{iv} ')'];
      end
        cbar_ppt.cmap=jet_sub;
        cbar_ppt.cmax=max(conlev.(varn));
        cbar_ppt.cmin=min(conlev.(varn));
        inc=round((cbar_ppt.cmax-cbar_ppt.cmin)/5);

        cbar_ppt.ytick=[cbar_ppt.cmin: inc :cbar_ppt.cmax];

       

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

       
       ww3_time_sw=SWC.ww3_timenum;
       ww3_time_dw=DW.ww3_timenum;

      %% ------ generate the sudo-space--space 2D matrix: ------ %%
       
         if d<dw_thres
             [ny,nx,nt]=size(indata.dpt);
            % [tmp,depID]=min(abs(d-indata.dpt(round(ny/2),:,round(nt/2))));
            % this part needs to be rewritten.
           dpt_vec=dpt_maxtrix(round(ny/2),:,1);
              
             if strcmp(varn,'fp');
                 indata.(varn)=1./indata.(varn);
             end
                
             [data_t, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(indata.(varn), ...
                 coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);
             
             [wnddata_t, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(indata.wnd_mag, ...
                coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000, stmtime);
                
           % if strcmp(varn,'Cd') ;
               % compute ratio:
               % interpolation from the DW grid to the grid of the shoaling domain:
               %Cd_bulk_interp=interp2(coordN_DW.XX,coordN_DW.YY,Cd_bulk_t, coordN.XX,coordN.YY);
              data_deep_interp=interp2(coordN_DW.XX,coordN_DW.YY,data_deep_t, coordN.XX,coordN.YY);
                
              % Cd2bulk_ratio= data_t./Cd_bulk_interp;
               data2deep_ratio= data_t./data_deep_interp;
           % end

         else  
            % It would be better to sync the deep water results with the
            % shallow water results. 
              if strcmp(varn,'fp');
                 indata.(varn)=1./indata.(varn);
             end
            [tmp,sid]=min(abs(ww3_time_dw-ww3_time_sw(1)));
            data_t=indata.(varn)(:,:,sid);
            wnddata_t=indata.wnd_mag(:,:,sid);
            
            coordN.XX = (coord.XX - stmx(sid)*1000)./1000;
            coordN.YY= (coord.YY - stmy(sid)*1000)./1000;
            
            
            
            %if strcmp(varn,'Cd');
                %@ 3. additional computation
               % [Cd_bulk_t]=compute_Cd_bulk_DW(DataDir.DW, meanCdFname,bulkchoice,wnddata_t);
                data_deep_t= data_t;
                coordN_DW=coordN;
                % compute the ratio also in the deep water domain.
               % Cd2bulk_ratio= data_t./Cd_bulk_t;
                data2deep_ratio= data_t./data_deep_t;
           % end
            
           
        end
      
      if vec_flag==1
        % wind vector:
        [udata_t{1}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(udata(1), ...
            coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime);            
        [vdata_t{1}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(vdata(1), ...
            coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime);
            
        % wind stress vector:
        [udata_t{2}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(udata(2), ...
            coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime);           
        [vdata_t{2}, TT,YT,XX_sudo,coordN]=convert_to_time_and_space_coordinate_updated(vdata(2), ...
            coord.XX,coord.YY,d, dpt_vec, ww3_timenum,stmx*1000,stmy*1000,stmtime); 
      end
        
      %%% modify coordinate:
      if isNom
          coordN.XX=coordN.XX./Rmax;
          coordN.YY=coordN.YY./Rmax;
      end
    %%% create time-space map of data:    
    
        
        labelstr.title=['depth=' num2str(d) 'm'];
       
      %%% make plot: 
   
        hfig{ifig}=figure(iv+length(VarNames)+ifig);
        set(hfig{ifig},'name',varn);
        
        ip=baseid+ncol*(rowid-1);

        
        if vec_flag == 1 
           generate_SpaceTime_map_at_deps(hfig{ifig},ncol,ip, coordN,data_t,wnddata_t,data_conlev,conlev.wnd_mag,...
                                     cbar_ppt,Rmax,labelstr,fontsz,dmsz,tickinc,pos,udata_t,vdata_t)
        else
           generate_SpaceTime_map_at_deps(hfig{ifig},ncol,ip, coordN,data_t,wnddata_t,data_conlev,conlev.wnd_mag,...
                                     cbar_ppt,Rmax,labelstr,fontsz,dmsz,tickinc,pos)
        end
        
        if 1==1
        % if (strcmp(varn,'Cd')) % plot ratio:
           cbar_ppt2.cmap=coolwarm_sub;      %'redblue';
           cbar_ppt2.cmax=1.5; %max(conlev.Cd_ratio);
           cbar_ppt2.cmin=0.5; %min(conlev.Cd_ratio);
           cbar_ppt2.tag='Ratio';%'Cd/Cd_{DW bulk}';
           cbar_ppt2.ytick=[0.5:0.1:1.5];
           hexf{1,ifig}=figure(iv*10+ifig);
           set(hexf{1,ifig},'name','ratio');
%            generate_SpaceTime_map_at_deps(hexf{1,ifig}, ncol,ip, coordN,Cd2bulk_ratio,wnddata_t,conlev.Cd_ratio, conlev.wnd_mag,...
%                                      cbar_ppt2,Rmax,labelstr,fontsz,dmsz,tickinc,pos)
%            
%           cbar_ppt2.tag='Cd/Cd_{DW ssd}';
%            hexf{2,ifig}=figure(iv*20+ifig);   
%           set(hexf{2,ifig},'name','Cd2deep_ratio');
           generate_SpaceTime_map_at_deps(hexf{1,ifig}, ncol,ip, coordN,data2deep_ratio,wnddata_t,conlev.Cd_ratio,conlev.wnd_mag,...
                                     cbar_ppt2,Rmax,labelstr,fontsz,dmsz,tickinc,pos)
%          %end
        end
          
    end    % depth loop
    
   % if rowid==2
        % save figure:
        for i=1:length(hfig);
            svfigname=[varn '_at_different_deps_with_3Blevs'];
            figure(hfig{i});
            set(gcf,'paperunits','inches','paperposition',[0 0 13 7.8]);
            print(gcf,'-djpeg','-r350',[figsvdir filesep svfigname '_' domainstr '_fig' num2str(i) '_' date '.jpg']);
            close gcf
        end
        
      %  if 1==1
       % if strcmp(varn,'Cd');
            for i=1:size(hexf,2)
                svfigname=[varn '2deep_ratio_at_different_deps_with_3Blevs'];
                figure(hexf{1,i});
                set(gcf,'paperunits','inches','paperposition',[0 0 13 7.8]);
                print(gcf,'-djpeg','-r350',[figsvdir filesep svfigname '_' domainstr '_fig' num2str(i) '_' date '.jpg']);
                close gcf
              
               
%                 svfigname='Cd2deep_ratio_at_different_deps_with_3Blevs';
%                 figure(hexf{2,i});
%                 set(gcf,'paperunits','inches','paperposition',[0 0 13 7.8]);
%                 print(gcf,'-djpeg','-r350',[figsvdir filesep svfigname '_' domainstr '_fig' num2str(i) '_' date '.jpg']);
%                 close gcf
            end
       % end
       % end
       
   % end
  
  
end        % variable loop

return