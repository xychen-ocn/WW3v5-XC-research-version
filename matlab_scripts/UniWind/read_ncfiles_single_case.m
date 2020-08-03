% Purpose: This script is used to aggregate data from both DW and SW experiments.

%%% 0a. define global variables:
global TRS_abspath LOC_abspath local_root scptdir trackdir
global wspd fetch FLD_arry tailev_arry expID fieldFN specFN wind_fac
%%%

qtty_pool={'hs','lm','fp','dp','dir','spr','dpt','wnd','ust'};   % qtties to be read from the above specified file.
qtty_longname={'SWH','MWL','Freq_{pk}','MWD','DirSpr','Depth','10-m Wind','Wind Stress'};
qtty_units={'m','m','s^(-1 )','degree','degree','m','m/s','N/m^2'};
scaler_flag=[1 1 1 1 1 1 1 0 0];                          % 1=scaler; used to identify whether the qtty is scaler or vector.

grav=9.806;
srcterm='ST4';
%% reading and preprocessing data:

% ------------------------------------------------------------------------%
% 1. set path and filenames for input and output files:
% ------------------------------------------------------------------------%
% data is read from opendap:
% data will be processed and stored in LOC_abspath

absfieldFN=[TRS_abspath filesep fieldFN];
absspecFN=[TRS_abspath filesep specFN];

svfname1=[expID '_standard_wvparm.mat'];
svfname2=[expID '_WW3spec_at_outpnts.mat'];



% ------------------------------------------------------------------------%
% 2. read in data and perform necessary preprocessing
% ------------------------------------------------------------------------%
if 1==0
%%% ------ o_o: start reading data from field output file -------------- %%
%--quantity pool/list: (use lower case for reference.)
%   hs, lm, fp, dp, dir, dpt, wnd, ust
for iqtty=1:length(qtty_pool)
    qttyn=qtty_pool{iqtty};
    if scaler_flag(iqtty)==1
        eval(['[hfig,' qttyn ', coord, time_field]=check_qtty_in_netCDFfile(absfieldFN,0,qttyn);']);
        
    else  % vector field (wnd, ust)
        uqttyn=['u' qttyn];
        vqttyn=['v' qttyn];
        eval(['[hfig,' qttyn ', coord, time_field]=check_qtty_in_netCDFfile(absfieldFN,0,uqttyn,vqttyn);']);
    end
end

%--Derived Quantity:
%  retrieve dimensions of qtties:
NX = size(coord.XX,2);
NY = size(coord.YY,1);
NT = length(time_field);

% convert coordinate to a relative coordinate: (in km)
% also, @ central transect, YY=0;
coord_rel.XX = (coord.XX-coord.XX(1,1))/1000;
coord_rel.YY = (coord.YY - 0.5*(coord.YY(1,1)+coord.YY(end,1)))/1000;

%  wind_dir, ust_dir, misalignment angle
%  Cp (peak wave phase speed), wave age;
%  ust/sqrt(gHs): another parameter similar to wave age, used to quantify
%  wave breaking.

% apply wind factor to correct wind:

wnd.u=wnd.u./wind_fac;
wnd.v=wnd.v./wind_fac;
wnd.mag=sqrt(wnd.u.^2+ wnd.v.^2);

winddir=atan2(wnd.v,wnd.u).*(180/pi);       % wind direction in cartesian coordinate and in degree.
winddir(winddir<0)=winddir(winddir<0)+360;  % and range from [0 360];
wnd.dir=winddir;
        
% drag coefficient:
ust.mag=sqrt(ust.u.^2 + ust.v.^2);
Cd=ust.mag.^2./wnd.mag.^2 .*10^3;      %Cd1000;

% convert wave direction to cartesian convention:
dp=get_cartesian_direction(dp,'Meteo');      % in degree.
dir=get_cartesian_direction(dir,'Meteo');       % mean wave direction.


save([LOC_abspath filesep svfname1],'hs','lm','fp','dp','dir','spr','dpt','wnd','ust','Cd','coord_rel','time_field','coord');
end

if 1==0
%%%------------- Finish reading data from field ouptut file ------------%%%
%-------------------------------------------------------------------------%
%%% -----o_o: start reading data from spectral point output file -------%%%
% get data from point output:
clear spec wvpar pntloc time_pnt
[spec, wvpar,pntloc,time_pnt]=extract_qtty_ww3pntNC(absspecFN);

% Cp, wave age, ust/sqrt(gHs)
% @ last time step:

%NPNT=length(pntloc.x);
NPNT=size(pntloc.x,1);
NT=length(time_pnt);

% use pntloc information to find peak frequency from the field WW3
% output instead of using fp from the output.   % This is not necessary for
% ST4.
% for ip=1:NPNT
%     [trash,yid]=min(abs(coord.YY(:,1)/1000-pntloc.y(ip)));
%     [trash,xid]=min(abs(coord.XX(1,:)/1000-pntloc.x(ip)));
%     pntloc_fp(ip) = fp(yid, xid,end);
%     %pause
% end
% spec.fp=pntloc_fp;
if strcmp(expID,'DW')
    kp=(2*pi*spec.fp).^2./9.8;         % corrected on Apr 2, 2019....
    kpi=(2*pi*spec.fpi).^2./9.8;         % corrected on Apr 2, 2019....

else
    kp=[];kpi=[];
    for it=1:NT
        for ip=1:NPNT
            kp(ip,it)=sig2wn_newton_iteration(2*pi*spec.fp(ip,it),wvpar.dpt(ip,it),0.0001);
            kpi(ip,it)=sig2wn_newton_iteration(2*pi*spec.fpi(ip,it),wvpar.dpt(ip,it),0.0001);
        end
    end
end
spec.kp = kp;
spec.kpi= kpi;
wvpar.pkwvlen=2*pi./kp;
spec.cp = 2*pi.*spec.fp./kp;
spec.cpi = 2*pi.*spec.fpi./kpi;
spec.cpU10 = spec.cp ./(wvpar.wnd./wind_fac);

% convert point location to relative location as in the field data:
pntloc_rel.x=(pntloc.x-coord.XX(1,1)/1000);       % correct, b/c pntloc.x is in km.
pntloc_rel.y=(pntloc.y - 0.5*(coord.YY(1,1)+coord.YY(end,1))/1000);

save([LOC_abspath filesep svfname2],'spec','wvpar','pntloc','pntloc_rel','time_pnt'); %
%%%------------- Finish reading data from point ouptut file ------------%%%
end
%-------------------------------------------------------------------------%
%%% The following is only performed on SW Exps.
%%% -----   o_o: start reading data from FLD output file         -------%%%
%fnames={'no_breaking';'with_Pbreaking'};
%if strcmp(expID, 'SW')
if 1==1
    load([LOC_abspath filesep  svfname2]);
    load([LOC_abspath filesep  svfname1]);
    NPNT=size(pntloc.x,1);
    
%  NT=length(time_pnt);
%for ii=1:length(fnames)
   % FLDfolder =fnames{ii};
   FLDfolder='';
% read diag_FLD for stress:
for FLDID=FLD_arry           % save as 2 separate files.
    FLDstr=num2str(FLDID);
    itlev=0;
    clear FLD_ust SPECX PROFS FLD_drv
    for tailev=tailev_arry;
        itlev=itlev+1;
        tailevstr=num2str(tailev,'%2.2i');
        tailev_value=tailev*10^(-3);
        
        % read out fldXX_diagout_BXX.nc for the
        % sea-state-dependent stress.
        if length(tailev_arry)==1 && tailev==1
            header='T';
        else
            header='B';
        end
        diagFN=[TRS_abspath filesep  'fld' FLDstr '_diagout_' header tailevstr '.nc'];
        % disp(diagFN);
        % pause
        % write a function to read diagout_BXX.nc file:
        [fld_ust,spec2Dx,profs]=get_fld_diagnc(diagFN, FLDID);
        
        FLD_ust{itlev}=fld_ust;
        SPECX{itlev}=spec2Dx;
        if FLDID == 1
            PROFS{itlev}=profs;
        else
            PROFS{itlev}=NaN;
        end
        
        NT1=size(fld_ust.magqced,2);
        NT2=size(wvpar.wnd,2);
        NTcom=min(NT1,NT2);
        % calculate Drag:
        if size(fld_ust.magqced,1) ~= size(wvpar.wnd,1)
            NPT2=min(size(fld_ust.magqced,1),size(wvpar.wnd,1));
            Cd= (fld_ust.magqced.^2./wspd.^2).*10^3;
            MAA_wndust=wdir-fld_ust.dirqced;
        else           
            Cd= (fld_ust.magqced(:,1:NTcom).^2./(wvpar.wnd(:,1:NTcom)./wind_fac).^2).*10^3;
            % misaglinment angle:
            MAA_wndust=wvpar.wnddir(:,1:NTcom)-fld_ust.dirqced(:,1:NTcom);     % > 0: ust to the right of wind vector; <0 ust to the left of wind vector;
            NPT2 = NPNT;
        end
       id=find(MAA_wndust>350);
        MAA_wndust(id)=MAA_wndust(id) - 360;
        id=find(MAA_wndust<-350);
        MAA_wndust(id)=MAA_wndust(id) + 360;
        
        
        
        %  if size(spec.cp,1)== size(fld_ust.magqced,1)
       % disp(['NT2=' num2str(NPT2)]);
            
        wvage = spec.cp(1:NPT2,1:NTcom)./ fld_ust.magqced(1:NPT2,1:NTcom);
        wvagei = spec.cpi(1:NPT2,1:NTcom)./ fld_ust.magqced(1:NPT2,1:NTcom);

        %  else
        %      wvage = spec.cp./ fld_ust.magqced';
        %      if iscolumn(wvage)==0
        %          wvage=wvage';
        %      end
        %  end
        
        eqiwvage = fld_ust.magqced(1:NPT2,1:NTcom)./sqrt(grav.*wvpar.Hs(1:NPT2,1:NTcom));    % equvilant inverse waveage
        
        FLD_drv{itlev}.Cd=Cd;
        FLD_drv{itlev}.MAA_wndust=MAA_wndust;
        FLD_drv{itlev}.Cp=spec.cp;
        FLD_drv{itlev}.CpU10=spec.cpU10;
        FLD_drv{itlev}.wvage=wvage;
        FLD_drv{itlev}.wvagei=wvagei;
        FLD_drv{itlev}.eqiwvage=eqiwvage;
        
    end
    
    if length(tailev_arry)==1
        svfname3=[expID '_FLD' FLDstr '_BGFDLcor_diagout.mat'];
    else 
        svfname3=[expID '_FLD' FLDstr '_diagout.mat'];
    end

    %  iii. save data:
    % save diagnositc FLD output:
    %save([LOC_abspath filesep 'test_tail_attachment/Fmean_2.5x_3.5x/' svfname3], 'FLD_ust','SPECX','PROFS','FLD_drv');
    save([LOC_abspath filesep svfname3], 'FLD_ust','SPECX','PROFS','FLD_drv');
    
end

end
%end