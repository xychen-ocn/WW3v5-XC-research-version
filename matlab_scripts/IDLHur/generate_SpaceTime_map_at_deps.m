function generate_SpaceTime_map_at_deps(hfig,ncol,ip, coordN,data_t,wnddata_t,...
                                   data_conlev,wnd_conlev,cbar_ppt,Rmax,...
                                   labelstr,fontsz,dmsz,tickinc,pos, varargin)
 switch nargin
     case 15
         vec_flag=0;
     case 17
         vec_flag=1;
         uvec=varargin{1};
         vvec=varargin{2};
     otherwise
         disp('too much input argument!');
         return
 end
        
 % predefine parameters:
 nrow=2;
 cmap=cbar_ppt.cmap;
 cmin=cbar_ppt.cmin;
 cmax=cbar_ppt.cmax;
 
 incxy=1;  % km;
 dx=coordN.XX(1,2)-coordN.XX(1,1);
 dy=coordN.YY(2,1)-coordN.YY(1,1);
 
 yinc=floor(incxy/dy);
 xinc=floor(incxy/dx);
 
 [tmp, xstid]= min(abs(coordN.XX(1,:)-(-3+0.25)));
 [tmp,xedid]= min(abs(coordN.XX(1,:)-(3-0.25)));
 
 [tmp, ystid]= min(abs(coordN.YY(:,1)-(-3+0.25)));
 [tmp,yedid]= min(abs(coordN.YY(:,1)-(3-0.25)));
 

 
%  if mod(ip,ncol)==1
%      yinc=8;
%      xinc=8;
%  else
%      if ip/ncol<=1
%          yinc=floor(0.06.*size(coordN.XX,1));
%          xinc=floor(0.04.*size(coordN.XX,2));
%      else
%          yinc=floor(0.08.*size(coordN.XX,1));
%          xinc=floor(0.06.*size(coordN.XX,2));
%      end
%      
%  end
 
 if mod(ip,ncol)~=0
     labelspc=140;
 else
     labelspc=80;
 end
 
 
    
 figure(hfig);

    hsub{ip}=subplot(nrow,ncol,ip)  ;
    % map out input main data:
    pcolor(coordN.XX,coordN.YY,data_t); shading interp;
    hold on;
    
  % [C,hc]=contour(coordN.XX, coordN.YY, wnddata_t, wnd_conlev,'color',[0.38 0.38 0.38],'linestyle','--','linewidth',1.2);

    plot([0 dmsz],[0 0],'--','color','k','linewidth',1.0);
    % add contours of the main data:
   % if mod(ip,ncol)==1
   
   %% plot wind contour:
   if 1==0
       [C_wnd,hc_wnd]=contour(coordN.XX,coordN.YY,wnddata_t,wnd_conlev,'linestyle','--','linewidth',1.2, ...
           'color',[0.65 0.65 0.65],'Visible','off');
       
       [bs, es, levs]=clines(C_wnd);
       for ilev=1:length(wnd_conlev)
           sids=find(levs==wnd_conlev(ilev));
           if length(sids)>1
               % choose the one with longer circumference.
               CL=[];
               for il=1:length(sids)
                   l=sids(il);
                   tmp=bs(l):es(l);
                   CL(il)=length(tmp);
               end
               [CLmax, lid]=max(CL);
               l=sids(lid);
           else
               l=sids;
           end
           plot(C_wnd(1, bs(l):es(l)), C_wnd(2, bs(l):es(l)),'color',[0.75 0.75 0.75],'linestyle','-','linewidth',0.9);
           text(C_wnd(1, bs(l)), C_wnd(2, bs(l)),num2str(wnd_conlev(ilev)),'color',[0.75 0.75 0.75],'linestyle','-','linewidth',0.9);
           
       end
       
   end

       [C,hc]=contour(coordN.XX,coordN.YY,data_t,[data_conlev],'k','linewidth',1.2);
       clabel(C,hc,data_conlev,'labelspacing',labelspc);
       
       
       
   %    [C,hc]=contour(coordN.XX,coordN.YY,data_t,[2.0, 2.5,3.0],'k','linewidth',1.2);
  %     clabel(C,hc,[2.0, 2.5,3.0]);
      
     % clabel(C,hc,[1:0.2:4, 5],'labelspacing',160);
   % else
   %     [C,hc]=contour(coordN.XX,coordN.YY,data_t,data_conlev,'k','linewidth',1.2);
   %     clabel(C,hc,data_conlev,'labelspacing',50);
   % end
   %  clabel(C,hc,'manual');
    % add contours of wind speed:
    %figure(1);
    
    

   % clabel(C_wnd,hc_wnd,wnd_conlev,'labelspacing',300, 'color',[0.65 0.65 0.65]);
    % clabel(C,hc,'manual');
     
    % add vectors at request
    

    % add circle denoting Radius of maximum wind
  %  circle(0,0,1,[0.40 0.40 0.40],1.1);
  circle(0,0,1,'k',2);
    % circle(0,0,2,[0.40 0.40 0.40],0.9);
    %  circle(0,0,3,[0.40 0.40 0.40],0.9);
    plot(0,0,'+k','linewidth',1.0);
   % plot([0 0],[-3 3],'--','color',[0.40 0.40 0.40],'linewidth',1.0);
    

    % add Colorbar
    colormap(hsub{ip},cmap);
   % jet_sub=jet(64);
   % jet_sub=jet_sub(7:3:64,:);
   %jet_sub=jet_sub(17:2:52,:);
   %colormap(jet_sub);
    hb=colorbar(hsub{ip});
    if mod(ip,ncol) ==0
        hb.XLabel.String=cbar_ppt.tag;
        hb.XLabel.FontSize=fontsz.cbar;
%     elseif ip==1
%         hb.XLabel.String='Cd (x1000)';
%         hb.XLabel.FontSize=fontsz.cbar;
    else
        hb.Visible='off';
    end
    
    caxis([cmin, cmax]);   
    
   % if mod(ip, ncol)~=1
      set(hb,'visible','on','ytick',cbar_ppt.ytick); 
   % end

    % label axis:
    set(gca,'fontsize',fontsz.axis);
   % if mod(ip,ncol)==1
        ylabel(labelstr.y);
  %  end
    
    %if (ip<=ncol)
       title(labelstr.title,'fontsize',fontsz.title,'fontweight','bold');
   % end
               
    
  %  if (ip>2*ncol)
       xlabel(labelstr.x);
  %  end
    %hold off

    xlim([-dmsz dmsz]);
    ylim([-dmsz dmsz]);
    set(gca,'xtick',-dmsz:tickinc:dmsz, ...
        'ytick',-dmsz:tickinc:dmsz);
    
    if vec_flag==1
        % wind vector
        refvec=false;
        refcol{1}=[0.1 0.1 0.1]; %[0.635, 0.078, 0.1840];
        refcol{2}=[0.4940, 0.1840, 0.5560];
        quiver(coordN.XX(ystid:yinc:yedid,xstid:xinc:xedid),coordN.YY(ystid:yinc:yedid,xstid:xinc:xedid), ...
            double(uvec{1}(ystid:yinc:yedid,xstid:xinc:xedid)*0.025),double(vvec{1}(ystid:yinc:yedid,xstid:xinc:xedid)*0.025),'color',refcol{1},'AutoScale','off');
        %
        %         for i=length(uvec)
        %             ncquiverref(double(coordN.XX(1:yinc:end,1:xinc:end)),double(coordN.YY(1:yinc:end,1:xinc:end)), ...
        %                 double(uvec{i}(1:yinc:end,1:xinc:end)),double(vvec{i}(1:yinc:end,1:xinc:end)), ...
        %                 'm/s','max',refvec,refcol{i})
        %             hold on;
        %         end
        
        if 1==0
            % wind stress vector
            quiver(coordN.XX(1:yinc:end,1:xinc:end),coordN.YY(1:yinc:end,1:xinc:end), ...
                double(uvec{2}(1:yinc:end,1:xinc:end)*0.15),double(vvec{2}(1:yinc:end,1:xinc:end)*0.15),'color',refcol{2},'AutoScale','off');
        end
    end
    axis('square');
    d_ax1=daspect;
    set(gca,'pos',pos{ip});
%     if vec_flag==1
%        % add legend represents the scale of the vector:
%        %%% 1. build another axis:
%        a1=gca;
%        a2pos=[a1.Position(1), a1.Position(2)+a1.Position(4) a1.Position(3)*0.4, a1.Position(4)*0.2];
%        a2=axes('pos',a2pos);
%        a2.XLim=[-2 2];
%        a2.YLim=[-2 2];
%        a2.XTick='';
%        a2.YTick='';
%        a2.XTickLabel='';
%        a2.YTickLabel='';
%        a2.XLim=[-2 2];
%        a2.Box='on';
%        
%        %%% 2. build loaction of the vector in the new axis.
%        xref=[0];
%        yref=[0];
%        
%        tmp=sqrt( uvec{1}.^2 +vvec{1}.^2);
%        uref_wind=-round(max(tmp(:))/10)*10;
%        vref_wind=0;
%        
%        %%% 3. plot the vector:
%        hq1=quiver(xref,yref, uref_wind*0.05, vref_wind*0.05,'AutoScale','off','color','k','linewidth',1.2);
%        hq1.MaxHeadSize=5;
%        
%       % text(xref(1), yref(2),'wind','color','m','horizontalalignment','center','verticalalignment','bottom');
%        text(xref,yref,[num2str(abs(uref_wind)) ' m/s'],'verticalalignment','bottom','color','k','horizontalalignment','center');
%        
%        daspect(a2,d_ax1);
%     end

    %% plot extra legend:


    
      
return
