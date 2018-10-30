clc;
clear all;
close all;



%% Script to generate data and response
% fn = @(x1,x2) [1 + (x1 + x2 + 1).^2.*(19 - 14*x1 + 3*x1.^2 - 14*x2 + 6*x1.*x2 + 3*x2.^2)]...
%                 .*[30 + (2*x1-3*x2).^2.*(18-32*x1+12*x1.^2+48*x2-36*x1.*x2+27*x2.^2)];
% 
% YE = @(x) (6.*x-2).^2.*sin(12.*x-4);
% YC = @(x) 0.5.*YE(x)+10.*(x-0.5)+5;
% %YC = @(x) (6.*x-2).^2.*sin(12.*x-4);        

load('FTX_det_2D');
load('FTX_det_3D');
MM = [3];
for i = 1:1
M = MM(i);
xStep = 0.875;
tauP  = 28;
posX = [[0:xStep:tauP]./tauP]';
Nx   = length(posX);

%x1en = [1 [sort(ceil(lhsdesign(33,1).*33))]'];
%x1en = [1  3     8    11    17    19    24    26    31];
x1en = [1 1:32];
x1en = [1 8 17 26];

x1cn = [1 1:2:32];
X1E = posX(x1en,1);
X1C = posX(x1cn,1);


YE =  [FTX_det_3D(x1en,M)]; 
YC =  [FTX_det_2D(x1cn,151+M) ; FTX_det_2D(x1en,151+M)];


%% ------------------ initializing parameters----------------------
saveFlag = 1;   
Xe_P_XcFlag = 0;
gpr_F_RMS_2D  = gpr_mdl_Krig;
gpr_F_RMS_SUR = gpr_mdl_coKrig;
gpr_TRUE      = gpr_mdl_Krig;

gpr_mdl_2D    = gpr_F_RMS_2D;
gpr_mdl_SUR   = gpr_F_RMS_SUR;
gpr_mdl_YE    = gpr_TRUE;


lbX = [0];
ubX = [1];

X1En = normalizeX(X1E,lbX,ubX);
X1Cn = normalizeX(X1C,lbX,ubX);

% Define the correlation function
% 1 : Gaussian
% 2 : Gaussian x Periodic
% 3 : GEK
% 4 : Gaussian symmertic
% 5 : Periodic
RcorrType = 5;

ntheta = 1;
nP = 1;
nT = 0;
k = size(X1Cn,2);
rsb = 1; % the index for subplot. rsb :rows, csb:colums
csb = 1;
effectType = 'LIN' ; % other can be 'RAND'
predictorStr = {'x','h_p','l_{C1}','h_{C1}','w_{C1}','h_N','a'};
numXnew = 50; % number of random points for plotting
intCon = [];

%Limits for hyperparamsfor krig mdl
lbTheta_P = [1e0.*ones(1,ntheta)  1.991.*ones(1,nP)   ] %         ];     % Lower Bound of Variables
ubTheta_P = [  5e1*ones(1,ntheta)   1.999.*ones(1,nP) ] %        ];    % Upper bounds

%Limits for hyperparams for coKrig mdl
lbTheta_P_rho = [1e0.*ones(1,ntheta) 1.991.*ones(1,nP) -10] %         ];     % Lower Bound of Variables
ubTheta_P_rho = [5e1*ones(1,ntheta)   1.999.*ones(1,nP)  10] %        ];    % Upper bounds


numTest   = 1;
numSmpl = size(X1Cn,1) - numTest;
numSmpl3DTest = 1;
numSmpl3D =size(X1En,1) - numSmpl3DTest;

extractMmLHD;

%% ---------------------- same for any kind of data----------------
dummy = 0;
if (not(exist('gpr_2D.mat')))
        save gpr_2D.mat dummy
end
if (not(exist('gpr_3D.mat')))
        save gpr_3D.mat dummy
end
%% load the MmLHD 
load('index.mat')
ind3DTest = index.ind3DTest;
ind3D     = index.ind3D;
ind2DTest = index.ind2DTest;
ind2D     = index.ind2D;

%% Script to assign data and response
objVar2DStr = 'YC';
objVar3DStr = 'YE';
objVar2D = eval(objVar2DStr);
objVar3D = eval(objVar3DStr);


if (Xe_P_XcFlag == 1)
x1c     = [X1Cn(ind2D,:) ; X1En(ind3D,:) ];
x1cTest = [X1Cn(ind2DTest,:) ; X1En(ind3DTest,:)];
x1e     = X1En(ind3D,:);
x1eTest = X1En(ind3DTest,:);

YC_C     = [objVar2D(ind2D,:) ; objVar2D(numSmpl+numTest+ind3D,:)] ; 
YC_CTest = [objVar2D(ind2DTest,:); objVar2D(numSmpl+numTest+ind3DTest,:)]; 

YE_ETest = objVar3D(ind3DTest,1);
YE_E = objVar3D(ind3D,:); 
YC_E = objVar2D(numSmpl+numTest+ind3D,1);
YC_ETest = objVar2D(numSmpl+numTest+ind3DTest,1);

else
    x1c     = [X1Cn(ind2D,:)];
    x1cTest = [X1Cn(ind2DTest,:)];
    x1e     = X1En(ind3D,:);
    x1eTest = X1En(ind3DTest,:);

    YC_C     = [objVar2D(ind2D,:)] ; 
    YC_CTest = [objVar2D(ind2DTest,:)]; 

    YE_ETest = objVar3D(ind3DTest,1);
    YE_E = [objVar3D(ind3D,:)]; 
    YC_E = objVar2D(numSmpl+numTest+ind3D,1);
    YC_ETest = objVar2D(numSmpl+numTest+ind3DTest,1);

end


if (size(x1c,2)==1)
    figure(987),hold all;
    plot(x1c,YC_C,'-o');
    plot(x1e,YE_E,'-x');
    plot(x1e,YC_E,'-^');grid on,xlabel('x'),ylabel('Response'),
    legend('show'),lgnd = legend('sample YC_C','sample YE_E','sample YC_E');
end
    
gpr_mdl_2D = initKrig(gpr_mdl_2D,x1c,x1cTest,YC_C,YC_CTest,lbX,ubX,RcorrType);
gpr_mdl_SUR = initCoKrig(gpr_mdl_SUR,x1e,x1eTest,YE_E,YE_ETest,lbX,ubX,gpr_mdl_2D,RcorrType);    
gpr_mdl_YE   =  initKrig(gpr_mdl_YE,x1e,x1eTest,YE_E,YE_ETest,lbX,ubX,RcorrType);

%% -----------------Optimization of log-likehood for low fidelity--------------------------
nvars = ntheta + nT + nP;
nPop = 1000;
MaxIt = 75;
%%Always set the bounds first at large values to see the trend, and then
%%reduce it
lb = [lbTheta_P];      % Lower Bound of Variables
ub = [ubTheta_P];      % Upper bounds

gpr_mdl_2D = OPT_ConLL(gpr_mdl_2D,nvars,lb,ub,MaxIt,nPop,intCon);
gpr_mdl_YE  = OPT_ConLL(gpr_mdl_YE,nvars,lb,ub,MaxIt,nPop,intCon);

if saveFlag == 1    
     save gpr_2D.mat gpr_F_RMS_2D gpr_mdl_YE -append;
     %save gpr_mdl_YE.mat gpr_mdl_YE
else
end
%% ----------------- Evaluation of low fidelity model---------------
if saveFlag == 1 
    load('gpr_2D');
    gpr_mdl_2D = gpr_F_RMS_2D;
else
end
gpr_mdl_2D.crossValidate;
gpr_mdl_2D.paramEffects(rsb,csb,effectType,predictorStr,numXnew);
gpr_mdl_YE.paramEffects(rsb,csb,effectType,predictorStr,numXnew);

%% ---------------------- loading the low fidelity GP------------------------------
if saveFlag == 1 
    load('gpr_2D');
    gpr_mdl_2D = gpr_F_RMS_2D;
else
end

%% -----------------Optimization of log-likehood for d = ye - rho*yc--------------------------
nvars =  ntheta + nT + nP + 1;
A = [];
b = [];
Aeq = [];
beq = [];
nPop = 5000;
MaxIt = 75;

lb = [lbTheta_P_rho];      % Lower Bound of Variables
ub = [ubTheta_P_rho]; 

[gpr_mdl_SUR, gpr_mdl_ERROR]= OPT_ConLLD(YC_E,gpr_mdl_SUR,gpr_mdl_2D,nvars,lb,ub,MaxIt,nPop,intCon);


if saveFlag == 1    
     save gpr_3D.mat gpr_F_RMS_SUR -append
else
end
%% ------------------ Prediction of coKriging model-------------------
if saveFlag == 1 
    load('gpr_3D');
    gpr_mdl_SUR = gpr_F_RMS_SUR;
    gpr_mdl_2D = gpr_F_RMS_2D;
else
end
gpr_mdl_SUR.crossValidate;
gpr_mdl_SUR.paramEffects(rsb,csb,effectType,predictorStr,numXnew);



%%
plottingGPs






end


