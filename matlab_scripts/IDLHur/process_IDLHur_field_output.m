function process_IDLHur_field_output(absFN, process_time_range, process_freq, ...
                                     wind_fac, DepParms, ...
                                     MatFname, local_dir, TrackInfo, expID,tailev)
% Description:
%
%    This script performs Quality Controls on dataset and also separate data into different wave types.
%    Finish debugging on Jan 21, 2019
%    corrected the fpi-related parameteris. (fpi is radian frequency(omega)
%    already)

%% -------------- Section I --------------------- %%
%@@ define local variables:
grav=9.8;
rad2dgr=180/pi;
hs_threshold=0.5;           % hs <0.5 will be NaN
wnd_threshold=10;           % wnd<10 will be NaN

%trackFN=['stmcenter_cartesian_and_time.mat'];
trackfile=TrackInfo.trackFN;
track_dir=TrackInfo.dir;
dateformat='yyyymmdd HHMMSS';

% wind speed:
xbinwidth=2;
xrange=[10 60];  
if strcmp(expID, 'SW')
    SW_flag=true;
else
    SW_flag=false;
end

% cut-off wavenumber in model(FLD) and in satellite observation:
k_max=366.8313;
k_cutoff=9.686;                     % mean for GPS L-band.;

% accuracy threshold in computing phase speed:
eps=0.0001;


% wave types:
wave_types={'WS','CSP','CSN','OS'};

%@@ make local directory for data saving if exsit:
if ~exist(local_dir,'dir'); 
    feval('mkdir', local_dir);
end


%@@ compute some necessary local variables:
deps_top=DepParms.deps_q.*(1+DepParms.dep_binwidth);
deps_btm=DepParms.deps_q.*(1-DepParms.dep_binwidth);


time_start=process_time_range{1};
time_end=process_time_range{2};

xinc=process_freq.x;
yinc=process_freq.y;
tinc=process_freq.t;

%% -------------- Section II --------------------- %%
disp('--> Preparing ....');
%@@ 1. read data from THREDDS data server:
%%%   data to be extracted is predefined via InqVar:
% InqVar={'hs','lm','dpt','fp','dp','u01','u02','u03','u04','wnd','ust'}; %tws','spr',
% OutVar={'hs','lm','dpt','fp','pwdir','QC_flag','fpi','mssx','mssy','wnd','ust'};%'dirspr', 

InqVar={'hs','dpt','u01','u02','wnd','ust'}; %tws','spr',
OutVar={'hs','dpt','QC_flag','fpi','wnd','ust'};%'dirspr', 

isvec=false(size(InqVar));
isvec(end-1:end)=true;

starts=[1 1 1];
counts=[Inf Inf 144];
FuncName='check_qtty_in_netCDFfile';
    for iv=1:length(InqVar)
        
        if isvec(iv)
        
          varnx=['u' InqVar{iv}];
          varny=['v' InqVar{iv}];
          [hf, odata,coord, ww3_timenum]=feval(FuncName,absFN,0, varnx,varny);
          ODATA.(OutVar{iv})=odata;
          
        else
         
            [hf, odata,coord, ww3_timenum]=feval(FuncName,absFN,0, InqVar{iv});
            ODATA.(OutVar{iv})=odata;
            
        end
    end

     clear hf
%@@ 2. time information management: 
%%%  obtain the indices of the start time and end time 
%%%  output data between selected start time and end time will be processed 

    [trash,tsid]=min(abs(ww3_timenum-datenum(time_start,dateformat)));
    [trash,teid]=min(abs(ww3_timenum-datenum(time_end,dateformat)));
    ww3_tres=(ww3_timenum(2)-ww3_timenum(1))*24;   %unit: in hour;
    clear trash
%%% load IDLHur track data (center location and time)

  %  [stmx,stmy,stm_timenum]=read_ASCII_trackInfo(track_dir, trackfile);
   load([track_dir filesep trackfile]);
    TCtrack_tres=(stm_timenum(2)-stm_timenum(1))*24;                  % unit: hour;   (a.k.a. 15min)
    tres_inc=round(ww3_tres/TCtrack_tres);

%%% align the first entry in the storm time to the ww3 output time:
    [tdif,id1]=min(abs(stm_timenum-ww3_timenum(1)));
        if tdif > TCtrack_tres/24.0;
           disp('Error in matching the TCtrack output time and ww3 output time!');
           return
        end

        if tres_inc ~= 0
            % TCloc information is now increased in the same freq. as the ww3 output data
            stmx=stmx(id1:tres_inc:end);     
            stmy=stmy(id1:tres_inc:end); 
            stmtime=stm_timenum(id1:tres_inc:end);
        else
            disp({'Warning: ww3_data output freq. higher than TCtrack info output freq!';...
                'interpolation is used'});
            % do interpolation to the storm track:
            stmx_org=stmx;
            stmy_org=stmy;
            stmtime_org=stm_timenum;
            
            stmtime=stmtime_org(id1):ww3_tres/24:stmtime_org(end);
            stmx=interp1(stmtime_org,stmx_org,stmtime);
            stmy=interp1(stmtime_org,stmy_org,stmtime);
            %return
        end


%% -------------- Section III --------------------- %%
disp('--> Processing ....');                              

% use only the data between processing time range:
    for iv=1:length(InqVar);
        if isvec(iv)
           ODATA.(OutVar{iv}).u=ODATA.(OutVar{iv}).u(:,:,tsid:tinc:teid);
          ODATA.(OutVar{iv}).v=ODATA.(OutVar{iv}).v(:,:,tsid:tinc:teid);
          
        else
          ODATA.(OutVar{iv})=ODATA.(OutVar{iv})(:,:,tsid:tinc:teid);
        end
    end
    
    % extract information within selected time range.
    ww3_timenum=ww3_timenum(tsid:tinc:teid);
    stmx=stmx(tsid:tinc:teid);
    stmy=stmy(tsid:tinc:teid);
    stmtime=stmtime(tsid:tinc:teid);

%@@ a. define and derive all needed quantites:
%   touch on wind
    ODATA.wnd.v=(ODATA.wnd.v)./wind_fac;
    ODATA.wnd.u=(ODATA.wnd.u)./wind_fac;
    wnd_mag=sqrt((ODATA.wnd.u).^2+(ODATA.wnd.v).^2);
    wndir=atan2((ODATA.wnd.v), (ODATA.wnd.u));

%   touch on wind stress
    ust_mag=sqrt((ODATA.ust.u).^2+(ODATA.ust.v).^2);
    ustdir=atan2((ODATA.ust.v),(ODATA.ust.u));
    if 1==0
%   touch on peak wave direction
    pwdir_cart=get_cartesian_direction(ODATA.pwdir,'Meteo');    % already in degree.
    pwdir_cart(pwdir_cart>180)= pwdir_cart(pwdir_cart>180)-360;
    pwdir_cart=pwdir_cart/rad2dgr;         % in radian
    ODATA.pwdir=pwdir_cart;
    
%   compute mean square slope in satelite observation     
%%@ when attachment wavenumber is smaller than cutoff wavenumber:
    delta_mss=tailev*10^(-3)*log(k_max./k_cutoff);
    mss_tot= ODATA.mssx + ODATA.mssy;
    mss_measured=mss_tot-delta_mss;
    end
    
    %%%  misalignment angle between wind and wind stress:   (in degree)
    misang_wust=(wndir-ustdir)*rad2dgr;
    misang_wust(misang_wust>180)=misang_wust(misang_wust>180)-360;
    misang_wust(misang_wust<-180)=misang_wust(misang_wust<-180)+360;
        
    if 1==0
    %%%  misalignment angle between wind stress and peak wave:
    misang_ustwv=(pwdir_cart-ustdir)*rad2dgr;
    misang_ustwv(misang_ustwv>180)=misang_ustwv(misang_ustwv>180)-360;
    misang_ustwv(misang_ustwv<-180)=misang_ustwv(misang_ustwv<-180)+360;
    
     %%%  misalignment angle between wind  and peak wave:
    misang_wndwv=(pwdir_cart-wndir)*rad2dgr;
    misang_wndwv(misang_wndwv>180)=misang_wndwv(misang_wndwv>180)-360;
    misang_wndwv(misang_wndwv<-180)=misang_wndwv(misang_wndwv<-180)+360;
    end

    %%%  kp and phase speed of dominant waves:
    [ny,nx,nt]=size(ODATA.hs);
    % pre-allocate kp to increase the speed.
    if 1==0
    kp=zeros(ny,nx,nt);kpi=kp;
    if ~SW_flag               % deep water
        cp=grav./(2*pi*ODATA.fp);   %phase speed of the peak waves in deep water
        cpi=grav./(ODATA.fpi); %phase speed of the waves at peak input freq. in deep water

        kp = (2*pi*ODATA.fp).^2./grav;
        kpi= (ODATA.fpi).^2./grav;           % Did I correct this? (Yes).
    else                      % shallow water
        % call function;     % will take a long time.... 
        for it=1:nt
            for ix=1:nx
                for iy=1:ny
                    kp(iy,ix,it)=sig2wn_newton_iteration(2*pi*ODATA.fp(iy,ix,it),ODATA.dpt(iy,ix,it),eps); 
                    kpi(iy,ix,it)=sig2wn_newton_iteration(ODATA.fpi(iy,ix,it),ODATA.dpt(iy,ix,it),eps);
                end
            end
        end
        cp = 2*pi.*ODATA.fp./kp;
        cpi=ODATA.fpi./kpi;
    end


    %%%  wave age: (cp*cos<cp,ust> / ust)    (input wave age here.)
    inpwvage=cpi./ust_mag;      % wave age.
    wvage=cp./ust_mag;
    % including the angle is more appropriate:
%     wvage=cp.*cosd(misang_ustwv)./ust_mag;
%     wvage_T18=wnd_mag.* cosd(misang_wndwv)./cp;

    %%%  misalignment angle between wind and peak waves.
    dtheta=pwdir_cart - wndir;   % misalignment angle between wind and peak waves.  (in radian)
    end
    
    % put derived variables into ODATA structure:
    ODATA.wnd.mag=wnd_mag;
    ODATA.wnd.dir=wndir;
    
    ODATA.ust.mag=ust_mag;
    ODATA.ust.dir=ustdir;
    
    ODATA.misang_wust=misang_wust;
if 1==0
    ODATA.misang_ustwv=misang_ustwv;
    ODATA.misang_wndwv=misang_wndwv;
    
    ODATA.dtheta=dtheta;
    ODATA.inpwvage=inpwvage;
    ODATA.wvage=wvage;
    ODATA.cp = cp;
    ODATA.kp = kp;
    ODATA.mss_tot=mss_tot;
    ODATA.mss_measured=mss_measured;
end   
    %ODATA.swell_indx=(1-ODATA.tws)./ODATA.tws;
    
% save everything:
%save([local_dir filesep 'DWdata_unmasked.mat'], 'ODATA','-v7.3');

%%------------- Quality Control Start.. -----------------%%
disp('--> Quality Controling.. ');

    %@@ generate QC mask for data:
    QC_mask=true(ny,nx,nt);
    QC_mask(ODATA.QC_flag~=1)=false;

    hs_mask=true(ny,nx,nt);
    wnd_mask=true(ny,nx,nt);

    hs_mask(ODATA.hs<hs_threshold)=false;
    wnd_mask(wnd_mag<wnd_threshold)=false;

    hs_mask=hs_mask.*(~isnan(ODATA.hs));     % non-nan=1 and will be remained 1 where hs_mask==1.
    wnd_mask=wnd_mask.*(~isnan(wnd_mag));

    % take subset of wnd and hs mask:
    QC_mask=wnd_mask.*hs_mask.*QC_mask;

    data_mask=ones(size(QC_mask));
    data_mask(QC_mask==false)=NaN;

    %% !!!***************************************************!!!   %%
    
   fieldN=fieldnames(ODATA);
    %%%@   apply mask to data extracted and derived above: 
    for iv=1:length(fieldN)
        varn=fieldN{iv};
        
        if strcmp(varn,'wnd') || strcmp(varn,'ust');
            QCed_fieldVar.([varn '_mag'])= data_mask .* ODATA.(varn).mag;
            QCed_fieldVar.(['u' varn])= data_mask .* ODATA.(varn).u;
            QCed_fieldVar.(['v' varn])= data_mask .* ODATA.(varn).v;
            QCed_fieldVar.([varn '_dir'])= data_mask .* ODATA.(varn).dir;
        else
            QCed_fieldVar.(varn)= data_mask .* ODATA.(varn);
        end
           
    end
    
    %% !!!***************************************************!!!   %%

    % computed extra variables 
    QCed_fieldVar.Cd=( (QCed_fieldVar.ust_mag).^2) ./ ((QCed_fieldVar.wnd_mag).^2).* 1000;
%     QCed_fieldVar.ucp=QCed_fieldVar.cp .* cos(QCed_fieldVar.pwdir);
%     QCed_fieldVar.vcp=QCed_fieldVar.cp .* sin(QCed_fieldVar.pwdir);
    
%@  Save QC controlled data:
disp({'--> Finished Quality Control.. Saving data to '; ...
      local_dir});
      
ReadMe_String={['Note:dtheta is the misalignment angle (in radian) between wind and peak waves'];...
               ['pwdir is in radian and the direction is interpret from cartesian coord. .']}; 
       

% save everything:
save([local_dir filesep MatFname{1}], 'ODATA','QCed_fieldVar','ww3_timenum','coord', ...
                                    'stmtime','stmx','stmy','ReadMe_String','-v7.3');

% update stmtime information:
%save([local_dir filesep MatFname{1}],'stmtime','stmx','stmy');
if 1==0                           
%% ------------ Section IV  ------------------------- %%  
disp('--> Sorting data based on wave type ....'); 

%% 2.2 preprocessing on selected data(tsid:teid) and perform data sorting:
% Log:
% Jul 2, 2018: added the mean square slope, and misalignment angle, 
%              fp/U10 (this fp should also only denotes the incoming deep water fp, in the shoaling zone, it is fine, because fp~constant, but not fine for the surf zone.) 
% Jul 3, 2018: added the shallowness parameter (kh) into the preprocessing:
% fp/U10, kh should be added with a more complex data referncing method.
% (storm following and in deep water.)
% ----------------------------------------------------------------------- %

% add relative location information to each grid point:
for it=1:length(stmx)
    coord_rel.XX(:,:,it)=coord.XX/1E3 - stmx(it);
    coord_rel.YY(:,:,it)=coord.YY/1E3 - stmy(it);
end


%  i. select data during certain period: 24-hr as default.
% time dependent variables:
fieldN = fieldnames(QCed_fieldVar);

for i=1:length(fieldN)
    tmp=QCed_fieldVar.(fieldN{i});
    if SW_flag
        DataSub.(fieldN{i}) = tmp(1:yinc:end,1:xinc:end,:);   % How many data in time should I take?
    else
        DataSub.(fieldN{i}) = tmp(1:yinc:end,1:xinc:end,1:2:48);
    end
end
clear tmp

if SW_flag
    DataSub.XX=coord_rel.XX(1:yinc:end,1:xinc:end,:);
    DataSub.YY=coord_rel.YY(1:yinc:end,1:xinc:end,:);
else
    DataSub.XX=coord_rel.XX(1:yinc:end,1:xinc:end,1:2:48);
    DataSub.YY=coord_rel.YY(1:yinc:end,1:xinc:end,1:2:48);
end

fieldN = fieldnames(DataSub);
%@@  reshape all the quantities into vector form.
for i=1:length(fieldN)
    DataSubVec.(fieldN{i})=reshape(DataSub.(fieldN{i}),1,[]);
end


%@@ extract only non-NaN data:
QC_mask=~isnan(DataSubVec.(fieldN{1}));            %QC_mask==1: non NAN;

for i=1:length(fieldN)
    DataSubVec.(fieldN{i})=DataSubVec.(fieldN{i})(QC_mask==1);
end

%------- Doubting if this block is necessary... ------- %
% compute mean Cd in deep water:
if ~SW_flag
    DW_Cd_vec=DataSubVec.Cd;
    
    %%% 1. sort data by wind speed bin:
    [DW_Cd_binned, DW_U10_binned, center_wspd, DW_idx_binned]=sort_data_by_windspeed(DW_Cd_vec,DataSubVec.wnd_mag,xbinwidth,xrange);
    num_of_wspdbins=length(DW_Cd_binned);
    
    %%% 2. compute the mean Cd in deep water for each wind speed bin:
    DW_mean_Cd=[]; DW_median_Cd=[];
    for i=1:num_of_wspdbins
        DW_mean_Cd(i)=mean(DW_Cd_binned{i},'omitnan');
        DW_median_Cd(i)=median(DW_Cd_binned{i},'omitnan');
    end
    DW_Cd.mean=DW_mean_Cd;
    DW_Cd.median=DW_median_Cd;
    DW_Cd.U10 = center_wspd;
    clear DW_mean_Cd DW_median_Cd center_wspd 
    save([local_dir filesep 'DW_Cd.mat'],'DW_Cd');
  
end
% ---------------------------------------------------------- %

%% iv. further separate ssd stress based on different quantities
%% -- a. wind-pkwave misalignment angle (dth_QCed) -----  %%%
%%% ----> 1) build mask for marking wave characteristics:
% Note: Apr 29, the "wind sea" is incorrectly technically speaking.
% following swell. % will stick to this version for now.
%wave_mask=NaN(size(DataSubVec.(fieldN{1})));  % temperary mask.

dth_QCed=DataSubVec.dtheta;    % vector form and only contains valid data.

% category criteria
% Wind Sea or following swell:
wave_mask=NaN(size( dth_QCed ));  % temperary mask.
wave_mask( cos(dth_QCed) >= sqrt(2)/2 ) =1;
mask.WS=wave_mask;

% Cross-Swell Positive:
wave_mask=NaN(size(dth_QCed));  % temperary mask.
wave_mask( cos(dth_QCed) < sqrt(2)/2) =1;
wave_mask( cos(dth_QCed) < 0 ) =NaN;
mask.CSP = wave_mask;

% Cross-Swell Negative:
wave_mask=NaN(size(dth_QCed));  % temperary mask.
wave_mask( cos(dth_QCed) < 0) =1;
wave_mask( cos(dth_QCed) <= -sqrt(2)/2 ) =NaN;
mask.CSN = wave_mask;

% Opposing-Swell:
wave_mask=NaN(size(dth_QCed));
wave_mask( cos(dth_QCed) <= -sqrt(2)/2) =1;
mask.OS = wave_mask;


%%% ----> 2) apply masks and separate data: <-----  %%%

for ic=1:length(wave_types);
    wvcn=wave_types{ic};
    
    for iv = 1:length(fieldN)
    
        data_tmp=DataSubVec.(fieldN{iv}) .* mask.(wvcn);
        WvType_sorted.(wvcn).(fieldN{iv})= data_tmp(mask.(wvcn)==1);
        
    end
    
    % compute the percentage of each category :   
    masksize=length(mask.(wvcn));
    n=length(mask.(wvcn)==1);
    WvType_sorted.(wvcn).percent= n/masksize;
  
end

save([local_dir filesep MatFname{2}],'DataSubVec','WvType_sorted','-v7.3');

end




return