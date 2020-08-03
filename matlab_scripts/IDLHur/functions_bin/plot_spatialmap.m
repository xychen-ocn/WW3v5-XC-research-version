function [ handles ] = plot_spatialmap( x,y,var,scaler_flag,size_var,fieldN,...
                                     casename,cvec,axismode, xylabel) 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% do multiple subplots.

% define necessary parameters to make seamless subplots:
plotheight=20;
plotwidth=18;
leftedge=1.2;
rightedge=0.4;   
topedge=1;
bottomedge=1;
spacex=0.2;
spacey=0.2;
fontsize=5;   
inc=15;   % increment of vector.

% get dimension of the input structure variable.
D=[];
num_of_dim=length(size_var);
if num_of_dim ==2
    D(1)=1;
    for i=1:num_of_dim;
        D(i+1)=size_var(i);
    end
elseif num_of_dim==3
    for i=1:num_of_dim;
        D(i)=size_var(i);
    end
else
    disp('input variable can only have 3 dimensions in max.')
    return
end

ncol=D(3);
nrow=D(2); 

% determine cminval and cmaxval.
% need to loop through all the fields..

for i=1:D(1);        % first dimension of the
    hf{i}=figure(i); clf(hf{i});
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure.
    set(gcf,'units','inches');
    figpos=get(gcf,'position');
    set(gcf, 'PaperUnits', 'inches');
    %set(gcf, 'PaperSize', [plotwidth plotheight]);
    %set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', figpos);
    jj=0; 
    for j=D(2):-1:1;
        jj=jj+1;
        max_resv=[]; min_resv=[];
        for k=1:D(3);
            ip=k+(j-1)*D(3);
            
            if num_of_dim==3
                eval(['z=var(i).' fieldN{end-1} '(j).' fieldN{end} '(k).mag(:,:,end);']);
            elseif num_of_dim==2
                eval(['z=var(j).' fieldN{end} '(k).mag(:,:,end);']);
            end
            
            if scaler_flag~=1
                if num_of_dim==3
                    eval(['zdir=var(i).' fieldN{2} '(j).' fieldN{3} '(k).dir(:,:,end);']);
                elseif num_of_dim==2
                    eval(['zdir=var(j).' fieldN{2} '(k).dir(:,:,end);']);
                end
                %eval(['uz=var(i).' fieldN{2} '(j).' fieldN{3}
                %'(k).u(:,:,end);']);  % for testing purpose
                %eval(['vz=var(i).' fieldN{2} '(j).' fieldN{3} '(k).v(:,:,end);']);
                uz=z.*cosd(zdir);
                vz=z.*sind(zdir);
            end 
                
            max_resv(k)=max(max(z(:,5:end)));
            min_resv(k)=min(min(z(:,5:end)));
            
            hsb{ip}=subplot(D(2),D(3),ip);
            
            % get subplot position:
            pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,...
                topedge,ncol,nrow,spacex,spacey);
            
            set(hsb{ip},'pos',pos{jj,k},'Fontsize',fontsize) %,'XGrid','off','XMinorGrid','off','Box','on','Layer','Top');
            
            % use pcolor to make plot.
           % hf{id}=subplot(nrow,ncol,idx);
           % set(hf{id},'pos',pos{idx});
            pcolor(x,y,z); shading interp
            hold on
            [C,h]=contour(x,y,z,cvec,'color','k');
            if scaler_flag ~=1 
                quiver(x(1:inc:end,1:inc:end),y(1:inc:end,1:inc:end),uz(1:inc:end,1:inc:end),vz(1:inc:end,1:inc:end),0.25,'k');
            end
                
            clabel(C,h,cvec,'Fontsize',8,'LabelSpacing',200);
            colormap(hsb{ip},jet(40));
           % hbar{ip}=colorbar;
          %  set(hbar{ip},'Visible','off');
            
            % colorbar visibility:
            if k==D(3) 
                hbar{ip}=colorbar;
                % reinforce subplot size:
                %set(hsb{ip},'Position',pos{jj,k});
                set(hbar{ip},'Visible','on');
               % special_cbid=ip;
            end
            
            if j==1
                title(casename{k})
            end
            
            if k==1 %&& j == round(D(2)/2);
                ylabel(xylabel{2},'fontsize',12);
            end
            if j==D(2) %&& k== round(D(3)/2);
                xlabel(xylabel{1},'fontsize',12);
            end
            
            if j~=D(2) 
                set(gca,'xticklabel',[]);
            end
            
            if k~=1
                set(gca,'yticklabel',[]);
            else
                set(gca,'ytick',[min(y(:,1)):10:max(y(:,1))]);
            end
                
            axis(axismode)
           % caxis([cminval cmaxval]) 
            
        end
        cminval_raw=min(min_resv);
        cmaxval_raw=max(max_resv);
        
        dif=cmaxval_raw-cminval_raw
        
        if dif <=1
            cminval=cminval_raw
            cmaxval=cmaxval_raw
            
        elseif cmaxval_raw<=3 && cmaxval_raw>=0.5
            cmaxval=round(cmaxval_raw);
            cminval=round(cminval_raw);
        else
            if cminval_raw~=0
                scl=get_scale(abs(cminval_raw));
                cminval=round(cminval_raw/scl*10)*scl/10
            else
                cminval=0;
            end
            scl=get_scale(abs(cmaxval_raw));
            cmaxval=round(cmaxval_raw/scl*10)*scl/10
        end
        
        for ip=1+(j-1)*D(3):j*D(3)
            axes(hsb{ip})
            caxis([cminval cmaxval]);
        end
    end
    
    
    % save handles:
    handles(i).fig=hf;
    handles(i).subplot=hsb;
    handles(i).cbar=hbar;
end




end







