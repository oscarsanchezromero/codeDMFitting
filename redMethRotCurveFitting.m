1;
clear();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Data adressing and variable declaration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Smoothed hybrid Halpha-HI rotation curves
%% Data files from High-resolution rotation curves of low surface brightness galaxies: Data
%% McGaugh, S.S., Rubin, V.C., & de Blok, W.J.G. 2001, AJ, 122, 2381
%% http://astroweb.case.edu/ssm/data/
%%
%% The code reads data from isolated ascii files. 
%% Example: the "data_format_F563_1.txt" file should contain
%%
%%      1.40000      3.30000      17.2000      0.00000      45.1000      24.7000
%%      2.80000      8.00000      20.6000      0.00000      81.6000      17.5000
%%      6.40000      17.5000      26.4000      0.00000      90.0000      24.8000
%%      7.90000      ...          ...          ...          ...          ...
%% 
%% where columns are:  radius (kpc), rotation caused by the observed gas (km/s),
%% rotation caused by the observed disk stars (km/s), rotation caused by the observed 
%% bulge stars (km/s), the observed [smoothed] velocity (km/s) and uncertainty in V (km/s).
%%
%% File name data definition and location
catalog= ['DDO185';'DDO189';'DDO47';'DDO52';'ESO0140040';'ESO0840411';'ESO1200211';'ESO1870510';'ESO2060140';'ESO3020120';...
'ESO3050090';'ESO4250180';'ESO4880049';'F563_1';'F568_3';'F571_8';'F579_V1';'F583_1';'F583_4';'F730_V1';'IC2233';'LSBCF563';'NGC100';...
'NGC1560';'NGC2366';'NGC3274';'NGC4395';'NGC4455';'NGC5023';'U11454';'U11557';'U11583';'U11616';'U11648';'U11748';'U11819';...
'UGC10310';'UGC1230';'UGC1281';'UGC3137';'UGC3371';'UGC4173';'UGC4325';'UGC5005';'UGC5750';'UGC711'];
sourcepath='./data_lsb';    % Folder with data files
resultspath='./fittings/';  % Folder created in advance for results saving. 
prefixfilename='/data_format_';
extensionfilename='.txt';


global radii; 
global vrot;
global weights;
global CteDim=10/((3.0856776^2)*4.51697); %
% This constant changes if mass, spaces and velocity units of data are different
% from solar masses, kpc and km/s respectively.
global vbary;
global somenullvbary;
global totalnullvbary;

tol=10^(-2); % tolerance for determining interval of search of varphi minimizers


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code flux keys (select as true the variables if the corresponding action should be realized) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Halo model fitting to be computed (all of them are allowed at the same time)
%
fittingISO=true;   % Isothermal model fitting   
%fittingISO=false;
fittingBUR=true;  % Burkert model fitting
%fittingBUR=false;
fittingNFW=true;  % NFW model fitting
%fittingNFW=false;
%
% Graphics to be shown during the script execution (optional)
%
%varphiGraph=true;  % varphi(s) function
varphiGraph=false;  
%bestFitGraph=true; % best fitted rotation curve
bestFitGraph=false;
%
% Fitting results printing
% 
%printFittings=false; 
printFittings=true;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bucle for multiple galaxies fitting   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for galaxyfile=1:length(catalog);

  %%%%%%%%%%%%%%%%%%
  % data reading  %
  %%%%%%%%%%%%%%%%%%
  data=load(strcat(sourcepath,prefixfilename,catalog(galaxyfile,:),extensionfilename));
  mkdir(strcat(resultspath,catalog(galaxyfile,:)));
  %% Source data units:
  %%Column 1: r (kpc); Column 2: vgas (km/s); Column 3: disk (km/s)
  %%Column 4: vbul (km/s); Column 5: vrot (km/s); Column 6: e_vrot (km/s) : error
  radii=data(:,1);
  vrot=abs(data(:,5));
  errors=data(:,6);
  vbary=sqrt(data(:,2).^2+data(:,3).^2+data(:,4).^2);
  n=length(radii);
  % Detection of total or partial null rotational velocity due to barionic matter
  somenullvbary=false;
  totalnullvbary=false; 
  if (sum(vbary)==0) totalnullvbary=true; elseif(round(prod(vbary))==0) somenullvbary=true; endif   
  % suitable for \nu-parametric fittings. Changue nu value for more free parameter fits.
  nu=2; 
  if(n <=nu) error("More than two data are needed for the fitting process \n"); endif
  weights=1./((n-nu)*errors.^2);
  % Calculation of two first terms on the right hand side of equation 19, that
  % will be complemented by alphaMV***.m functions in order to define \varphi.
  vv=WeighProd(vrot,vrot,weights);
  vvbary=WeighProd(vbary,vbary,weights); 
  vones=ones(size(radii));  % vector of ones appearing in several expressions 
  savefitting=[]; % Variable for saving the galaxy results. 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Isothermal fitting
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (fittingISO==true && n-nu > 0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Computation of the limit values for varphi(s)  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % varphi limit value when s tends to \infty
    ginf=inline("ones(size(r))","r"); 
    % Function defining equation 37
    equationVLimInf=inline("WeighProd(ginf(radii),vones,weights) - WeighProd(vrot,ginf(radii)./sqrt(t*ginf(radii)+(vbary.^2)),weights)","t");
    % The solution to 37 with vbary=0, which is an upper bound for the solutions to 37
    XVbaryNull=(WeighProd(vrot,sqrt(ginf(radii)),weights)./WeighProd(ginf(radii),vones,weights)).^2;
    % When some rotational velocity due to baryonic matter is null equationVLimInf is not defined at t=0
    if (totalnullvbary==true) X=XVbaryNull;
    elseif(somenullvbary==true)  
      j=-3;
      while (equationVLimInf(10^(j)*XVbaryNull) > 0) j--; endwhile 
      X = fzero(equationVLimInf,[10^(j)*XVbaryNull,XVbaryNull]); % Solves equation 37
    elseif (equationVLimInf(0)>=0) % Checking condition 39
      X=0;
    else
      X = fzero(equationVLimInf,[0,XVbaryNull]); % Solves equation 37
    endif
    % Calculation of the limit value by using Lemma 1.1 and development 19
    varphiLimInf=vv+vvbary-2.*WeighProd(vrot,sqrt(X*ginf(radii)+(vbary.^2)),weights)+X.*WeighProd(ginf(radii),vones,weights); 
    %
    % varphi limit value when s tends to 0
    g0=inline("r.^2","r"); 
    % Function defining equation 41
    equationVLim0=inline("WeighProd(g0(radii),vones,weights) - WeighProd(vrot,g0(radii)./sqrt(t.*g0(radii)+(vbary.^2)),weights)","t");
    % The solution to 41 with vbary=0 is an upper bound for the solutions to 41
    XVbaryNull=(WeighProd(vrot,sqrt(g0(radii)),weights)./WeighProd(g0(radii),vones,weights)).^2;
    % When some rotational velocity due to baryonic matter is null equationVLim0 is not defined at t=0
    if (totalnullvbary==true) X=XVbaryNull;
    elseif (somenullvbary==true)   % In this case condition 42 never holds.
      j=-3;
      while (equationVLim0(10^(j)*XVbaryNull) > 0) j--; endwhile 
      X = fzero(equationVLim0,[10^(j)*XVbaryNull,XVbaryNull]); % Solves equation 41
    elseif (equationVLim0(0)>=0) % Checking condition 42
      X=0;
    else
      X = fzero(equationVLim0,[0,XVbaryNull]); % Solves equation 41
    endif
    % Calculation of the limit value by using Lemma 1.2 and development 19
    varphiLim0=vv+vvbary+X.*WeighProd(g0(radii),vones,weights)-2.*WeighProd(vrot,sqrt(X.*g0(radii)+(vbary.^2)),weights);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Minimization interval determination    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    interval=[-3,3];
    direction=-1;
    maxiter=0;
    while (direction!=0 && maxiter < 50)
      maxiter++;
      eval= abs(vv+vvbary+alphaMVISO(10.^(interval(1)+(-0.2:0.1:0.2)))-varphiLim0)./varphiLim0;
      test1= (sum(eval(3:5)) < tol);
      test2= (sum(eval(1:3)) < tol);
      if(test1==0 && test2==1 ) direction=0;
      elseif (test1==0 && test2==0) direction = -1; interval(1)=interval(1)+0.3*direction;
      elseif (test1==1 && test2==1) direction = 1; interval(1)=interval(1)+0.3*direction;
      else direction = -1; interval(1)=interval(1)+0.3*direction;
      endif
    endwhile
    direction=1;
    maxiter=0;
    while (direction!=0 && maxiter < 50)
      maxiter++;
      eval= abs(vv+vvbary+alphaMVISO(10.^(interval(2)+(-0.2:0.1:0.2)))-varphiLimInf)./varphiLimInf;
      test1= (sum(eval(1:3)) < tol);
      test2= (sum(eval(3:5)) < tol);
      if(test1==0 && test2==1) direction = 0;
      elseif (test1==0 && test2==0) direction = 1; interval(2)=interval(2)+0.3*direction;
      elseif (test1==1 && test2==1) direction = -1; interval(2)=interval(2)+0.3*direction;
      else direction = 1; interval(2)=interval(2)+0.3*direction;
      endif
    endwhile 
    if (interval(1) > interval(2)) 
      warning("The interval selection module has not worked properly");
      interval=[-3,5];
    endif
    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Varphi minimization:     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % The minimizer is aproximated in the constrained interval.
    %[sol,alphaSol,info2]=fminbnd('alphaMVISO',10.^interval(1),10.^interval(2)); 
    %info1 = min([varphiLim0,varphiLimInf,vv+vvbary+alphaSol])< min([varphiLim0,varphiLimInf]);
    % When info1= 1 the existence of best fitting can be assured, on the other hand
    % info2 recognizes the convergence of the minimization method.
    % It would be also possible to aproximate the minimizer in a bounded interval by using:
    [minAprox,iminAprox]=min(alphaMVISO(10.^[interval(1):0.1:interval(2)]));
    xminAprox=10^(interval(1)+(iminAprox-1)*0.1);
    [sol,alphaSol]=fminbnd('alphaMVISO',xminAprox/(10^0.2),xminAprox*(10^0.2));
    info1 = min([varphiLim0,varphiLimInf,vv+vvbary+alphaSol])< min([varphiLim0,varphiLimInf]);
      
   
    % save best fitting parameters in savefitting variable
    rh=1/sol;
    rhocero=rhoISO(sol);
    redChiSquare= vv+vvbary+alphaSol;
    savefitting=[savefitting; rhocero, rh, redChiSquare, varphiLim0, varphiLimInf, info1];
    %
    % Optional plot  of the varphi function and limit values in logaritmic scale 
    if (varphiGraph==true)
      s=10.^(interval(1):0.1:interval(2));
      loglog(s,vv+vvbary+alphaMVISO(s),s,ones(size(s))*varphiLimInf,s,ones(size(s))*varphiLim0);
      title("varphi Isothermal fitting");
      pause();
      clf;
    endif
    % Optional plot  of the varphi function in logaritmic scale    
    if(bestFitGraph==true)     
      hold on;
      errorbar(radii,vrot,errors,'*');
      plot(radii,sqrt(((vISO(radii,sol)).^2).*(ones(length(radii),1)*(rhocero./(CteDim*sol.^3)))+(vbary.^2)*ones(size(sol))));
      title("Isothermal best fitted rotation curve");
      print(strcat(resultspath,catalog(galaxyfile,:),'/fittingISOMV.pdf'),'-dpdf')
      hold off;
      pause(5);
      clf;
    endif 
  endif % end of Isothermal fitting

  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Burkert fitting
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 if (fittingBUR==true && n-nu > 0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Computation of the limit values for varphi(s)   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % varphi limit value when s tends to \infty
    ginf=inline("r.^(-1)","r"); 
    % Function defining equation 37
    equationVLimInf=inline("WeighProd(ginf(radii),vones,weights) - WeighProd(vrot,ginf(radii)./sqrt(t*ginf(radii)+(vbary.^2)),weights)","t");
    % The solution to 37 with vbary=0, which is an upper bound for the solutions to 37
    XVbaryNull=(WeighProd(vrot,sqrt(ginf(radii)),weights)./WeighProd(ginf(radii),vones,weights)).^2;
    % When some rotational velocity due to baryonic matter is null equationVLimInf is not defined at t=0
    if (totalnullvbary==true) X=XVbaryNull;
    elseif(somenullvbary==true)  
      j=-3;
      while (equationVLimInf(10^(j)*XVbaryNull) > 0) j--; endwhile 
      X = fzero(equationVLimInf,[10^(j)*XVbaryNull,XVbaryNull]); % Solves equation 37
    elseif (equationVLimInf(0)>=0) % Checking condition 39
      X=0;
    else
      X = fzero(equationVLimInf,[0,XVbaryNull]); % Solves equation 37
    endif
    % Calculation of the limit value by using Lemma 1.1 and development 19
    varphiLimInf=vv+vvbary-2.*WeighProd(vrot,sqrt(X*ginf(radii)+(vbary.^2)),weights)+X.*WeighProd(ginf(radii),vones,weights); 
    %
    % varphi limit value when s tends to 0
    g0=inline("r.^2","r"); 
    % Function defining equation 41
    equationVLim0=inline("WeighProd(g0(radii),vones,weights) - WeighProd(vrot,g0(radii)./sqrt(t.*g0(radii)+(vbary.^2)),weights)","t");
    % The solution to 41 with vbary=0 is an upper bound for the solutions to 41
    XVbaryNull=(WeighProd(vrot,sqrt(g0(radii)),weights)./WeighProd(g0(radii),vones,weights)).^2;
    % When some rotational velocity due to baryonic matter is null equationVLim0 is not defined at t=0
    if (totalnullvbary==true) X=XVbaryNull;
    elseif (somenullvbary==true)   % In this case condition 42 never holds.
      j=-3;
      while (equationVLim0(10^(j)*XVbaryNull) > 0) j--; endwhile 
      X = fzero(equationVLim0,[10^(j)*XVbaryNull,XVbaryNull]); % Solves equation 41
    elseif (equationVLim0(0)>=0) % Checking condition 42
      X=0;
    else
      X = fzero(equationVLim0,[0,XVbaryNull]); % Solves equation 41
    endif
    % Calculation of the limit value by using Lemma 1.2 and development 19
    varphiLim0=vv+vvbary+X.*WeighProd(g0(radii),vones,weights)-2.*WeighProd(vrot,sqrt(X.*g0(radii)+(vbary.^2)),weights);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Minimization interval determination    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    interval=[-1,1];
    direction=-1;
    maxiter=0;
    while (direction!=0 && maxiter < 50)
      maxiter++;
      eval= abs(vv+vvbary+alphaMVBUR(10.^(interval(1)+(-0.2:0.1:0.2)))-varphiLim0)./varphiLim0;
      test1= (sum(eval(3:5)) < tol);
      test2= (sum(eval(1:3)) < tol);
      if(test1==0 && test2==1 ) direction=0;
      elseif (test1==0 && test2==0) direction = -1; interval(1)=interval(1)+0.3*direction;
      elseif (test1==1 && test2==1) direction = 1; interval(1)=interval(1)+0.3*direction;
      else direction = -1; interval(1)=interval(1)+0.3*direction;
      endif
    endwhile
    direction=1;
    maxiter=0;
    while (direction!=0 && maxiter < 50)
      maxiter++;
      eval= abs(vv+vvbary+alphaMVBUR(10.^(interval(2)+(-0.2:0.1:0.2)))-varphiLimInf)./varphiLimInf;
      test1= (sum(eval(1:3)) < tol);
      test2= (sum(eval(3:5)) < tol);
      if(test1==0 && test2==1) direction = 0;
      elseif (test1==0 && test2==0) direction = 1; interval(2)=interval(2)+0.3*direction;
      elseif (test1==1 && test2==1) direction = -1; interval(2)=interval(2)+0.3*direction;
      else direction = 1; interval(2)=interval(2)+0.3*direction;
      endif
    endwhile
    if (interval(1) > interval(2)) 
      warning("The interval selection module has not worked properly");
      interval=[-3,5];
    endif

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Varphi minimization:     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
     % The minimizer is aproximated in the constrained interval.
    %[sol,alphaSol,info2]=fminbnd('alphaMVBUR',10.^interval(1),10.^interval(2)); 
    %info1 = min([varphiLim0,varphiLimInf,vv+vvbary+alphaSol])< min([varphiLim0,varphiLimInf]);
    % When info1= 1 the existence of best fitting can be assured, on the other hand
    % info2 recognizes the convergence of the minimization method.
    % It would be also possible to aproximate the minimizer in a bounded interval by using:
    [minAprox,iminAprox]=min(alphaMVBUR(10.^[interval(1):0.1:interval(2)]));
    xminAprox=10^(interval(1)+(iminAprox-1)*0.1);
    [sol,alphaSol]=fminbnd('alphaMVBUR',xminAprox/(10^0.2),xminAprox*(10^0.2));
    info1 = min([varphiLim0,varphiLimInf,vv+vvbary+alphaSol])< min([varphiLim0,varphiLimInf]);
    %
    % save best fitting parameters in savefitting variable
    rh=1/sol;
    rhocero=rhoBUR(sol);
    redChiSquare= vv+vvbary+alphaSol;
    savefitting=[savefitting; rhocero, rh, redChiSquare, varphiLim0, varphiLimInf, info1];
    %
    % Optional plot  of the varphi function and limit values in logaritmic scale 
    if (varphiGraph==true)
      s=10.^(interval(1):0.1:interval(2));
      loglog(s,vv+vvbary+alphaMVBUR(s),s,ones(size(s))*varphiLimInf,s,ones(size(s))*varphiLim0);
      title("varphi Burkert fitting");
      pause();
      clf;
    endif
    if(bestFitGraph==true)
      % Optional plot  of the varphi function in logaritmic scale
      hold on;
      errorbar(radii,vrot,errors,'*');
      plot(radii,sqrt(((vBUR(radii,sol)).^2).*(ones(length(radii),1)*(rhocero./(CteDim*sol.^3)))+(vbary.^2)*ones(size(sol))));
      title("Burkert best fitted rotation curve");
      print(strcat(resultspath,catalog(galaxyfile,:),'/fittingBURMV.pdf'),'-dpdf')
      hold off;
      pause(5);
      clf;
    endif 
  endif % end of Burkert fitting

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % NFW fitting
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (fittingNFW==true && n-nu > 0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Computation of the limit values for varphi(s)  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % varphi limit value when s tends to \infty
    ginf=inline("r.^(-1)","r"); 
    % Function defining equation 37
    equationVLimInf=inline("WeighProd(ginf(radii),vones,weights) - WeighProd(vrot,ginf(radii)./sqrt(t*ginf(radii)+(vbary.^2)),weights)","t");
    % The solution to 37 with vbary=0, which is an upper bound for the solutions to 37
    XVbaryNull=(WeighProd(vrot,sqrt(ginf(radii)),weights)./WeighProd(ginf(radii),vones,weights)).^2;
    % When some rotational velocity due to baryonic matter is null equationVLimInf is not defined at t=0
    if (totalnullvbary==true) X=XVbaryNull;
    elseif(somenullvbary==true)  
      j=-3;
      while (equationVLimInf(10^(j)*XVbaryNull) > 0) j--; endwhile 
      X = fzero(equationVLimInf,[10^(j)*XVbaryNull,XVbaryNull]); % Solves equation 37
    elseif (equationVLimInf(0)>=0) % Checking condition 39
      X=0;
    else
      X = fzero(equationVLimInf,[0,XVbaryNull]); % Solves equation 37
    endif
    % Calculation of the limit value by using Lemma 1.1 and development 19
    varphiLimInf=vv+vvbary-2.*WeighProd(vrot,sqrt(X*ginf(radii)+(vbary.^2)),weights)+X.*WeighProd(ginf(radii),vones,weights); 
    %
    % varphi limit value when s tends to 0
    g0=inline("r","r"); 
    % Function defining equation 41
    equationVLim0=inline("WeighProd(g0(radii),vones,weights) - WeighProd(vrot,g0(radii)./sqrt(t.*g0(radii)+(vbary.^2)),weights)","t");
    % The solution to 41 with vbary=0 is an upper bound for the solutions to 41
    XVbaryNull=(WeighProd(vrot,sqrt(g0(radii)),weights)./WeighProd(g0(radii),vones,weights)).^2;
    % When some rotational velocity due to baryonic matter is null equationVLim0 is not defined at t=0
    if (totalnullvbary==true) X=XVbaryNull;
    elseif (somenullvbary==true)   % In this case condition 42 never holds.
      j=-3;
      while (equationVLim0(10^(j)*XVbaryNull) > 0) j--; endwhile 
      X = fzero(equationVLim0,[10^(j)*XVbaryNull,XVbaryNull]); % Solves equation 41
    elseif (equationVLim0(0)>=0) % Checking condition 42
      X=0;
    else
      X = fzero(equationVLim0,[0,XVbaryNull]); % Solves equation 41
    endif
    % Calculation of the limit value by using Lemma 1.2 and development 19
    varphiLim0=vv+vvbary+X.*WeighProd(g0(radii),vones,weights)-2.*WeighProd(vrot,sqrt(X.*g0(radii)+(vbary.^2)),weights);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Minimization interval determination    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    interval=[-1,1];
    direction=-1;
    maxiter=0;
    while (direction!=0 && maxiter < 50)
      maxiter++;
      eval= abs(vv+vvbary+alphaMVNFW(10.^(interval(1)+(-0.2:0.1:0.2)))-varphiLim0)./varphiLim0;
      test1= (sum(eval(3:5)) < tol);
      test2= (sum(eval(1:3)) < tol);
      if(test1==0 && test2==1 ) direction=0;
      elseif (test1==0 && test2==0) direction = -1; interval(1)=interval(1)+0.3*direction;
      elseif (test1==1 && test2==1) direction = 1; interval(1)=interval(1)+0.3*direction;
      else direction = -1; interval(1)=interval(1)+0.3*direction;
      endif
    endwhile
    direction=1;
    maxiter=0;
    while (direction!=0 && maxiter < 50)
      maxiter++;
      eval= abs(vv+vvbary+alphaMVNFW(10.^(interval(2)+(-0.2:0.1:0.2)))-varphiLimInf)./varphiLimInf;
      test1= (sum(eval(1:3)) < tol);
      test2= (sum(eval(3:5)) < tol);
      if(test1==0 && test2==1) direction = 0;
      elseif (test1==0 && test2==0) direction = 1; interval(2)=interval(2)+0.3*direction;
      elseif (test1==1 && test2==1) direction = -1; interval(2)=interval(2)+0.3*direction;
      else direction = 1; interval(2)=interval(2)+0.3*direction;
      endif
    endwhile
    if (interval(1) > interval(2)) 
      warning("The interval selection module has not worked properly");
      interval=[-3,5];
    endif
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Varphi minimization:     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % The minimizer is aproximated in the constrained interval.
    %[sol,alphaSol,info2]=fminbnd('alphaMVNFW',10.^interval(1),10.^interval(2)); 
    %info1 = min([varphiLim0,varphiLimInf,vv+vvbary+alphaSol])< min([varphiLim0,varphiLimInf]);
    % When info1= 1 the existence of best fitting can be assured, on the other hand
    % info2 recognizes the convergence of the minimization method.
    % It would be also possible to aproximate the minimizer in a bounded interval by using:
    [minAprox,iminAprox]=min(alphaMVNFW(10.^[interval(1):0.1:interval(2)]));
    xminAprox=10^(interval(1)+(iminAprox-1)*0.1);
    [sol,alphaSol]=fminbnd('alphaMVNFW',xminAprox/(10^0.2),xminAprox*(10^0.2));
    info1 = min([varphiLim0,varphiLimInf,vv+vvbary+alphaSol])< min([varphiLim0,varphiLimInf]);
      
   
    % save best fitting parameters in savefitting variable
    rh=1/sol;
    rhocero=rhoNFW(sol);
    redChiSquare= vv+vvbary+alphaSol;
    savefitting=[savefitting; rhocero, rh, redChiSquare, varphiLim0, varphiLimInf, info1];
    %
    % Optional plot  of the varphi function and limit values in logaritmic scale 
    if (varphiGraph==true)
      s=10.^(interval(1):0.1:interval(2));
      loglog(s,vv+vvbary+alphaMVNFW(s),s,ones(size(s))*varphiLimInf,s,ones(size(s))*varphiLim0);
      title("varphi NFW fitting");
      pause();
      clf;
    endif
    % Optional plot  of the varphi function in logaritmic scale    
    if(bestFitGraph==true)     
      hold on;
      errorbar(radii,vrot,errors,'*');
      plot(radii,sqrt(((vNFW(radii,sol)).^2).*(ones(length(radii),1)*(rhocero./(CteDim*sol.^3)))+(vbary.^2)*ones(size(sol))));
      title("NFW best fitted rotation curve");
      print(strcat(resultspath,catalog(galaxyfile,:),'/fittingNFWMV.pdf'),'-dpdf')
      hold off;
      pause(5);
      clf;
    endif 
  endif % end of NFW fitting
  




if(printFittings==true)
  mkdir(strcat(resultspath,catalog(galaxyfile,:)));
  save('-ascii',strcat(resultspath,catalog(galaxyfile,:),'/fittingCLASSICS.txt'),'savefitting')
endif

endfor







