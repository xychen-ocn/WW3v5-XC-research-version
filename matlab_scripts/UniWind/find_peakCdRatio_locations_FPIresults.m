% This script is used to find the maximum Cd ratio in shallow water and the
% corresponding kp.
% Code Structure:
%    0. apply mask to data (limited water depth and kpD trend (shoaling)).
%    1. Quality Control on Cd.
%    2. find the maximum Cd, find the local maximum (peak) Cd and their occurring location in the nondimensional depth space.
%
% Author: Xuanyu Chen;
% Date  : Dec 4, 2018;  Improved at Dec 22, 2018; Corrected for FPI setup
%         on Apr.6 2019. modified for computing maximum Cd ratio instead.
%

%%% define global parameters:
global TRS_abspath LOC_abspath local_root scptdir trackdir
global wspd fetch FLD_arry tailev_arry expID TRID slp

dwsamemesh_flg=false;
%%% define local parameters:
infile0.SW=[expID '_standard_wvparm.mat'];
infile1.SW=[expID '_WW3spec_at_outpnts.mat'];

infile1.DW='DW_WW3spec_at_outpnts.mat';
%DWpath=[LOC_abspath(1:end-18) filesep 'DW_4km'];
DWpath=[LOC_abspath(1:95) filesep 'DW_4km'];
if dwsamemesh_flg
    DWpath1=[LOC_abspath(1:end) '_dw'];
end
%DWpath=[LOC_abspath '_dw'];

%DWpath=[LOC_abspath(1:103) filesep 'DW']  ; %For swell

if length(tailev_arry)==1
    if dwsamemesh_flg
        svfilename='PeakCdRatio_and_masked_crosshore_variables_vFPI_BGFDLcor_10percentThshold_new.mat';
    else
        svfilename='PeakCdRatio_and_masked_crosshore_variables_vFPI_BGFDLcor_10percentThshold.mat';
    end
else
    if dwsamemesh_flg
        svfilename='PeakCdRatio_and_masked_crosshore_variables_vFPI_new.mat';
    else
        svfilename='PeakCdRatio_and_masked_crosshore_variables_vFPI.mat';
    end
         
end

svdir=LOC_abspath;

varname={'Hs','dpt','MWL','kpD','cp','HpD'};

span=3;
% wspdcrt=15;   % critical wind speed used to separate data preprocessing procedure.
% if wspd>wspdcrt
%     span=5;
% else
%     span=3;   % window length for moving average.
% end

kpD_thres_higher=1.0;
kpD_thres_lower=0;

% coordinate adjustment/lineup:
if slp==2000 || (slp==1000 && fetch==200)
    offset=800-fetch;
else
    offset=0;
end
    
%% load shallow water data:
load([LOC_abspath filesep infile0.SW],'lm');
load([LOC_abspath filesep infile1.SW]);    % check the path here, I need to specify fetch directory.

%% load DW data:
if dwsamemesh_flg
    DWdata=load([DWpath1 filesep infile1.SW]);
else
    DWdata=load([DWpath filesep infile1.DW]);
end
cp_dw=DWdata.spec.cp(:,end);
% cp_dw_intp shares the same coordinate as cp_sw (cp_box);
% instead of interpolation, I can now have the absolute value:
if dwsamemesh_flg
    cp_dw_intp=cp_dw;
else
    cp_dw_intp=interp1(DWdata.pntloc.x(:,1), cp_dw, pntloc.x(:,1)+offset);
end
%

%% 
%%% Read all parameters out:
%%% a. lump standard parameters and spectrum (cross-shore variation)
Hs_box=wvpar.Hs(:,end) ;%mean(wvpar.Hs,2);
dpt_box=wvpar.dpt(:,end); %mean(wvpar.dpt,2);
PWL_box=wvpar.pkwvlen(:,end); %mean(wvpar.pkwvlen,2);
kp_box=spec.kp(:,end); %mean(spec.kp,2);
kpD_box=kp_box.*dpt_box;
cp_box=spec.cp(:,end); %mean(spec.cp,2);

% compute peak wave height:
%%% Initialize data:
Hp_box=zeros(size(Hs_box));
Hs2_tmp=zeros(size(Hs_box));
%%%
for ip=1:length(Hs_box)
    E2d=spec.Efth(:,:,ip,end); %mean(spec.Efth(:,:,ip,:),4);
    freq_bins=spec.freq.cen;
    fp=spec.fp(ip);
    U10dir=wvpar.wnddir(ip);    % double check, also check why subtracted by 180.
    [Hp_box(ip),Hs2_tmp(ip)]=get_windsea_Hp_from_2DSpectrum(E2d,freq_bins,spec.theta,fp,wspd,U10dir);
end
HpD_box=Hp_box./dpt_box;

% need to get mean wave length information from the field output.
[ny,nx,nt]=size(lm);
MWL_box=lm(round(ny/2)+1,end:-1:2,nt)';
% DirSpr_box(:,ic)=spr(round(ny/2)+1,end:-1:2,nt);


%%% develop mask:
% mask 1: depth < 2.5m:
msk1=ones(size(dpt_box));
msk1(dpt_box<=0.5)=NaN;

%%%% This one might not be needed if I restrict the section of the data when finding the maximum. 
% if wspd> wspdcrt  
%     msk1(end-1:end)=NaN;
% end

msk2=ones(size(dpt_box));
if wspd>5
% mask Cd where kd increases at shallow water depth.
mask_criterion=diff(kpD_box)./diff(dpt_box);
tmpid=find(mask_criterion<0);
if ~isempty(tmpid)
    valid= find( dpt_box(tmpid)<30);
    if ~isempty(valid)
        msk2(min(tmpid(valid))+1:end)=NaN;
    else
        msk2(tmpid)=NaN;
    end
    
end
end
mask=msk1.*msk2;

% make plot to check:
hfig=figure(1); clf;
set(hfig,'name','check mask');
plot(dpt_box, kpD_box,'ob');
hold on;
plot(dpt_box.*msk1,kpD_box.*msk1,'+r');
plot(dpt_box.*msk2,kpD_box.*msk2,'+m');
%pause
plot(dpt_box.*mask,kpD_box.*mask,'*k');
hold off;

% apply mask to data:
for iv=1:length(varname)
    eval([varname{iv} '_masked=' varname{iv} '_box.*mask;']);
    %  eval([varname{iv} '_masked(mask,:)=NaN;']);    
end
Cp_dw_masked=cp_dw_intp.*mask;

clear mss_box Cd_box mss_masked Cd_masked
for FLD=FLD_arry
    FLDstr=num2str(FLD);
    if length(tailev_arry)==1
        infile2.SW=[expID '_FLD' FLDstr '_BGFDLcor_diagout.mat'];
    else
        infile2.SW=[expID '_FLD' FLDstr '_diagout.mat'];
    end
       
    load([LOC_abspath filesep infile2.SW]);
    for itlev=1:length(tailev_arry)
        % examine at the itlev th tail level.

      %  tailev=tailev_arry(itlev);
        mss_box{FLD}(:,itlev)=SPECX{itlev}.mss.mag(:,end); %mean(SPECX{itlev}.mss.mag,2);
       % Cd_box{FLD}(:,itlev)=mean(FLD_drv{itlev}.Cd,2);
        Cd_box{FLD}(:,itlev)=FLD_drv{itlev}.Cd(:,end);

        if length(Cd_box{FLD}) ~= length(mask)
            mask=mask(1:end-1);
        end
        mss_masked{FLD}(:,itlev)=mss_box{FLD}(:,itlev).*mask;
        Cd_masked{FLD}(:,itlev)=Cd_box{FLD}(:,itlev).*mask;
        
    end
    
    %% for Deep Water:
    if length(tailev_arry)==1
        infile2.DW=['DW_FLD' FLDstr '_BGFDLcor_diagout.mat'];
    else
        infile2.DW=['DW_FLD' FLDstr '_diagout.mat'];
    end

    if dwsamemesh_flg
        load([DWpath1 filesep infile2.SW]);
    else
        load([DWpath filesep infile2.DW]);
    end
    
    Cd_dw=[]; Cd_dw_intp=[];
    
    
    for itlev=1:length(tailev_arry)
        %if dwsamemesh_flg
        Cd_dw(:,itlev)=FLD_drv{itlev}.Cd(:,end);
       % else
       % Cd_dw0(:,itlev)=tmp_deep.FLD_drv{itlev}.Cd(:,end);
         Cd_dw_intp(:,itlev)=interp1(DWdata.pntloc.x(:,1), Cd_dw(:,itlev), ...
                             pntloc.x(:,1)+offset,'linear');    
       % end
    end
    
    if dwsamemesh_flg
        Cd_ratio{FLD}=Cd_masked{FLD}./Cd_dw;
        Cd_ratio_unmask{FLD}=Cd_box{FLD}./Cd_dw;
    else
        Cd_ratio{FLD}=Cd_masked{FLD}./Cd_dw_intp;
        Cd_ratio_unmask{FLD}=Cd_box{FLD}./Cd_dw_intp;
    end
   
    
    figure(10); hold on;
    plot(DWdata.pntloc.x(:,1)+offset,Cd_dw,'-*k');
    hold on;
    %plot(DWdata0.pntloc.x(:,1),Cd_dw0,'-b');
    %plot(pntloc.x(:,1)+offset,Cd_dw_intp,'-b');
    plot(pntloc.x(:,1)+offset,Cd_box{FLD},'-m');  %+offset

    
    %hold off;
    %pause
    
  
end

hfig=figure(2); clf;
set(hfig,'name','check mask');
for itlev=1:length(tailev_arry);
    plot(dpt_box, Cd_box{FLD}(:,itlev),'ob');
    hold on;
    plot(dpt_box.*mask,Cd_masked{FLD}(:,itlev),'*k');
    plot(dpt_box.*mask,Cd_ratio{FLD}(:,itlev),'-');
    
end
hold off;



%% find the first peak of Cd:
for FLD=FLD_arry
 %   Cd_max{FLD}=NaN(1,3);
 %   kpD_max{FLD}=NaN(1,3);
 %   Cp_max{FLD}=NaN(1,3);
    
%    Cd_pk{FLD}=NaN(3,3);
%    kpD_pk{FLD}=NaN(3,3);
%    Cp_pk{FLD}=NaN(3,3);
    
    hf10=figure(10); clf;
    set(hf10,'WindowStyle','Docked');
    ax10=gca;
    for itlev=1:length(tailev_arry);
        
        %% use running mean to smooth data first:
        %
        % Cd_data_in is acutally Cd_ratio. 
        Cd_data_in=double(flip(Cd_ratio{FLD}(mask==1,itlev)));
        
        
        kpD_data_in=flip(kpD_masked(mask==1));
        Cp_data_in=flip(cp_masked(mask==1));    
        Cp_dw_in = flip(Cp_dw_masked(mask==1));
        HpD_data_in=flip(HpD_masked(mask==1));
        dpt_mask=dpt_box.*mask;
        dep_in=flip(dpt_mask(mask==1));
        
        plot(ax10,kpD_data_in,Cd_data_in,'*','color',[0.5 0.5 0.5],'linewidth',0.8);
        hold on
        xlim([0 1.2]);
        
%         crit=(kpD_data_in-kpD_thres_higher).*(kpD_data_in-kpD_thres_lower);
%         sectionID=find(crit<=0);      % within kpD_thres(lower--> high)
%         

        
        plot(ax10,kpD_data_in,Cd_data_in,'*b','linewidth',0.8);
        
       
         %% smooth data:
         if wspd>15
         n=floor(span/2);
            % Cd_data_sm=NaN(size(Cd_data));
            Cd_data_sm=[];kpD_data_sm=[];
            for i=1+n:length(Cd_data_in)-n
                Cd_data_sm(i-n)=mean(Cd_data_in(i-n:i+n));
                kpD_data_sm(i-n)=mean(kpD_data_in(i-n:i+n));   %kpD data is not moving-averaged.
            end
            Cp_data_sub=Cp_data_in(1+n:length(Cd_data_in)-n);      % Cp data is not moving-averaged. why??
            HpD_data_sub=HpD_data_in(1+n:length(Cd_data_in)-n);
            Cp_dw_sub=Cp_dw_in(1+n:length(Cd_data_in)-n); 
            dpt_data_sub=dep_in(1+n:length(Cd_data_in)-n); 
            
            Cd_data_ready=Cd_data_sm;
            kpD_data_ready=kpD_data_sm;
            Cp_data_ready=Cp_data_sub;
            Cp_dw_ready=Cp_dw_sub;
            HpD_data_ready=HpD_data_sub;
            dpt_data_ready=dpt_data_sub;
            
            %%% ploted out smoothed results to check:
            plot(ax10,kpD_data_sm,Cd_data_sm,'-k','linewidth',1.2);
           % ylim([0 1.1*max(Cd_data_sm)]);
           ylim([0 1.6]);
            
         else
             Cd_data_ready=Cd_data_in;
             kpD_data_ready=kpD_data_in;
             Cp_data_ready=Cp_data_in;
              Cp_dw_ready=Cp_dw_in;
             HpD_data_ready=HpD_data_in;
             dpt_data_ready=dep_in;
             
         end
            
        
        
        %% simply find the maximum:
        crit=(kpD_data_ready-kpD_thres_higher).*(kpD_data_ready-kpD_thres_lower);
        sectionID=find(crit<=0);      % within kpD_thres(lower--> high)
        whos Cd_data_ready
        
        [Cd_max, maxid]=max(Cd_data_ready(sectionID));      
        kpD_max=mean( kpD_data_ready( sectionID(maxid) ) );
        Cp_max=mean( Cp_data_ready( sectionID(maxid) ) );
        HpD_max=mean( HpD_data_ready( sectionID(maxid) ));

        
        % figure(10)
        plot(ax10,kpD_data_ready, Cd_data_ready,'-k');
        plot(ax10,kpD_data_ready( sectionID(maxid) ),Cd_max,'oc');       % check if there are multiple maximums
        plot(ax10,kpD_max, mean(Cd_max,'omitnan'),'+c');
        
        pause(0.1);
        
        
        CdRatio_pk{FLD}(itlev)=mean(Cd_max,'omitnan');
        kpD_pk{FLD}(itlev)=kpD_max;
        Cp_pk{FLD}(itlev)=Cp_max;
        Cp_dw_pk{FLD}(itlev)=mean( Cp_dw_ready( sectionID(maxid) ) );
        CpU10_dw_pk{FLD}(itlev)=mean( Cp_dw_ready( sectionID(maxid) ) )/wspd;
        CpU10_pk{FLD}(itlev)=Cp_max/wspd;
        HpD_pk{FLD}(itlev)=HpD_max;
        dpt_pk{FLD}(itlev)=mean(dpt_data_ready( sectionID(maxid)),'omitnan');
        
        
        %% find range of kpD where Cd_ratio > 1;     
        % only the first peak is recognized here.
        % the potential overshoot in Cd is masked.
        tmp=Cd_ratio{FLD}(:,itlev);      
        % build mask:
        range_mask=zeros(size(kpD_box));
        %% % set threshold for significant increase
       % threshold=1.05;      % 5% increase
        threshold=1.10;      % 10% increase
        %%
        range_mask(tmp>threshold)=1;
        range_mask(kpD_box>2)=0;
       % range_mask(dpt_box<10)=0;
        
        affected_range.kpD{FLD,itlev}=kpD_box(range_mask==1);
        affected_range.dpt{FLD,itlev}=dpt_box(range_mask==1);
        affected_range.cp_sw{FLD,itlev}=cp_box(range_mask==1);
        affected_range.cp_dw{FLD,itlev}=cp_dw_intp(range_mask==1);
        affected_range.CdRatio{FLD,itlev}=tmp(range_mask==1);
        
%         figure(11);clf;
%         plot(dpt_box,tmp,'*b');
%         hold on;
%         plot(dpt_box(range_mask==1), tmp(range_mask==1),'or');
%         pause(0.1)
        
        % now, find the range for overshoot Cd:
        tmp2=Cd_ratio_unmask{FLD}(:,itlev);  
        dumped_mask=isnan(tmp);
        
        overshoot_mask=zeros(size(tmp2));
        overshoot_mask(tmp2>1.1)=1;
        overshoot_mask=overshoot_mask.*dumped_mask;
        
        overshoot_range.kpD{FLD,itlev}=kpD_box(overshoot_mask==1);
        overshoot_range.dpt{FLD,itlev}=dpt_box(overshoot_mask==1);
        overshoot_range.cp_sw{FLD,itlev}=cp_box(overshoot_mask==1);
        overshoot_range.cp_dw{FLD,itlev}=cp_dw_intp(overshoot_mask==1);
        overshoot_range.CdRatio{FLD,itlev}=tmp2(overshoot_mask==1);
        
%         figure(11);clf;
%         plot(dpt_box,tmp2,'*b');
%         hold on;
%         plot(dpt_box(overshoot_mask==1), tmp2(overshoot_mask==1),'or');
%         pause(0.1)
%       
        
    end
end





for iv=1:length(varname);
    vn=varname{iv};
    eval(['masked_qtty.(vn)=' vn '_masked;']);
end

peak_qtty.CdRatio=CdRatio_pk;
peak_qtty.kpD=kpD_pk;
peak_qtty.Cp_sw=Cp_pk;
peak_qtty.Cp_dw=Cp_dw_pk;
peak_qtty.CpU10_sw=CpU10_pk;
peak_qtty.CpU10_dw=CpU10_dw_pk;
peak_qtty.HpinvD=HpD_pk;
peak_qtty.dpt=dpt_pk;

%%% save the QCed data, peak information.
save([svdir filesep svfilename],'Cd_masked','Cd_ratio','Cd_ratio_unmask','masked_qtty', 'mss_masked', ...
                                'peak_qtty','affected_range','overshoot_range');

% pause;
%close all;


