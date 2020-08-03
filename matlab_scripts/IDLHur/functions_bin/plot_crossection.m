function [ hf ] = plot_crossection( x,crsec_yid, var,size_var,...
                 fieldN, titlestr, casename,color,lw, marker,...
                           xylabel,svfigname)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

D=[];
num_of_dim=length(size_var);
% if num_of_dim ==2
%     D(1)=1;
%     for i=1:num_of_dim;
%         D(i+1)=size_var(i);
%     end
% elseif num_of_dim==3
    for i=1:num_of_dim;
        D(i)=size_var(i);
    end
% else
%     disp('input variable can only have 3 dimensions in max.')
%     return
% end




for i=1:D(1)   % 1st dimension of input var (usually tailev)
    hf{i}=figure(i); clf(hf{i});
        set(gcf, 'units','normalized','outerposition',[0 0 0.8 0.8]); % Maximize figure.
        set(gcf,'units','inches');
        figpos=get(gcf,'position');
        set(gcf, 'PaperUnits', 'inches');
        %set(gcf, 'PaperSize', [plotwidth plotheight]);
        %set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', figpos);

    for j=1:D(end-1);   % 2nd dimension  (wspd usually)
        
        for k=1:D(end);
            % extract data:
            if num_of_dim==3
                eval(['z=var(i).' fieldN{end-1} '(j).' fieldN{end} '(k).mag(crsec_yid,:,end);']);
            else
                eval(['z=var(j).' fieldN{end} '(k).mag(crsec_yid,:,end);']);
            end
            
            z=squeeze(z);
            z(1:5)=NaN;
            hl(k)=plot(x,z,'color',color{j},'linewidth',lw,'linestyle','none',...
                'marker',marker{k},'MarkerFaceColor',color{j},'MarkerEdgeColor','w',...
                'Markersize',8);
            hold on
            
            
        end
   

    end
 
    % text(double(0.98*x(depthID)),0.99*hs{1}(depthID,mid_ny,end),...
%     num2str(depth{1}(depthID,mid_ny,1),'%5.2f'));
%xlabel('Cross-shore Distance (km)');
xlabel(xylabel{1});
ylabel(xylabel{2});
title(titlestr,'fontsize',14);
set(gca,'xtick',[0:5:100]);   % depth
%set(gca,'xdir','reverse')
set(gca,'fontsize',14);
hld=legend(hl,casename);
set(hld,'loc','southoutside','orientation','horiz');
grid on
% figname='testST2_crossshores_Hs_at_3incoming_waves';
% saveas(gcf,[figname '.fig']);
% set(gcf,'paperunits','inches','paperposition',[0 0 8 6]);
% print(gcf,'-djpeg','-r350',[figname '.jpg']);

end

end

