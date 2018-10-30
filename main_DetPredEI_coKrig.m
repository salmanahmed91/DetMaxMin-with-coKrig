%clear all;
close all;
clc;

colors10Class = [lines(5)];
newDefaultColors = colors10Class;
set(gca, 'ColorOrder',newDefaultColors,'NextPlot','replacechildren')

%% load the necessary variables and models
load('gpr_3D.mat')
load('gpr_mdl_YE.mat')
gpr_mdl_SUR = gpr_F_RMS_SUR;
gpr_mdl_YE   = gpr_mdl_YE;

%initialization
rsb = 2; % the index for subplot. rsb :rows, csb:colums
csb = 1;
effectType = 'LIN' ; % other can be 'RAND'
predictorStr = {'x','h_p','l_{C1}','h_{C1}','w_{C1}','h_N','a'};
numXnew = 500; % number of random points for plotting
ntheta = size(gpr_mdl_SUR.X,2);
nP     = size(gpr_mdl_SUR.X,2);

%Limits for hyperparams for coKrig mdl
lbTheta_P_rho = [1e0.*ones(1,ntheta) 1.991.*ones(1,nP) -10] ;%         ];     % Lower Bound of Variables
ubTheta_P_rho = [15e1*ones(1,ntheta)   1.999.*ones(1,nP)  10] ;%        ];    % Upper bounds

minMaxFlag = 0;

iter = [3];
xNew = [gpr_mdl_SUR.lbX:0.001:gpr_mdl_SUR.ubX]';
ymaxTrue = [max(gpr_mdl_YE.Eval(xNew))];
yminTrue = [min(gpr_mdl_YE.Eval(xNew))];
ymaxIter = [max(gpr_mdl_SUR.Eval(xNew))];
yminIter = [min(gpr_mdl_SUR.Eval(xNew))];

j = [ 1 2 3 4 9 10 11 12];



for i = 1:10
    
MaxIt = 20;
nPop = 50;
intCon = [];
ymin = min(gpr_mdl_SUR.Y);
ymax = max(gpr_mdl_SUR.Y);
[maxEI_max, xMax] =  maxEI(ymax,gpr_mdl_SUR,'max');
[maxEI_min, xMin] =  maxEI(ymin,gpr_mdl_SUR,'min');
%[maxEI_min, xMin] =  maxProdEI(ymin,ymax,gpr_mdl_coYE);
% if maxEI_max(end) < maxEI_min(end)
%     minMaxFlag = 1;
% else
%     minMaxFlag = 0;
% end
%% evaluating the high fidelity model
%this section shall be replaced by results from 3D FEM
ye_xMax = gpr_mdl_YE.Eval(xMax); 
ye_xMin = gpr_mdl_YE.Eval(xMin);

%% update the HF coKrig model
if (minMaxFlag == 0)
        gpr_mdl_SUR.X   =  [gpr_mdl_SUR.X ; xMin    ];
        gpr_mdl_SUR.Y   =  [gpr_mdl_SUR.Y ; ye_xMin ];
        minMaxFlag = 1;
else
    if minMaxFlag == 1
        gpr_mdl_SUR.X   =  [gpr_mdl_SUR.X ; xMax   ];
        gpr_mdl_SUR.Y   =  [gpr_mdl_SUR.Y ; ye_xMax];
        minMaxFlag = 0;
    end
end

%% re calculate the hyperparameters for coKrig model
% the krig model of 2D sims remains the same

%% -----------------Optimization of log-likehood for d = ye - rho*yc--------------------------
nvars =  size(gpr_mdl_SUR.theta_P,2)+1 ;
A = [];
b = [];
Aeq = [];
beq = [];
nPop = 2000;
MaxIt = 50;
intCon = [];

lb = [lbTheta_P_rho];      % Lower Bound of Variables
ub = [ubTheta_P_rho]; 
YC_E = gpr_mdl_SUR.gpr_mdl_2D.Eval(gpr_mdl_SUR.X);

[gpr_mdl_SUR, gpr_mdl_ERROR]= OPT_ConLLD(YC_E,...
                                         gpr_mdl_SUR,...
                                         gpr_mdl_SUR.gpr_mdl_2D,...
                                         nvars,lb,ub,MaxIt,nPop,intCon);

%% ------------------ Prediction of coKriging model-------------------
% 
% gpr_mdl_SUR.crossValidate;
% gpr_mdl_SUR.paramEffects(rsb,csb,effectType,predictorStr,numXnew);


ymin = min(gpr_mdl_SUR.Y);
ymax = max(gpr_mdl_SUR.Y);
figure(100);hold all;
[Fpred,RMSE]= gpr_mdl_SUR.Eval(xNew);
subplot(4,4,j(i)),plot(xNew,gpr_mdl_YE.Eval(xNew),'--','linewidth',1);hold all;
subplot(4,4,j(i)),plot(gpr_mdl_SUR.X,gpr_mdl_SUR.Y,'^','markersize',8);
subplot(4,4,j(i)),plot(xNew,[Fpred],'linewidth',2);
subplot(4,4,j(i)),plot(xNew,[Fpred+2.*RMSE Fpred-2*RMSE],'--','linewidth',1);
xlabel('x'),ylabel('Response'),hold all,grid on
% legend('show'),hldng = legend('EI max','EI min');

EI_xMax  = EI(ymax,xNew,gpr_mdl_SUR,'max');
EI_xMin  = EI(ymin,xNew,gpr_mdl_SUR,'min');
subplot(4,4,j(i)+4),plot(xNew,[EI_xMax EI_xMin],'linewidth',2);
xlabel('x'),ylabel('EI(x)'),hold all,grid on;
% legend('show'),hldng = legend('EI max','EI min');

figure(232)
iter = [iter ; size(gpr_mdl_SUR.X,1)];
ymaxTrue = [ymaxTrue ; max(gpr_mdl_YE.Eval(xNew))];
yminTrue = abs([yminTrue ; min(gpr_mdl_YE.Eval(xNew))]);
ymaxIter = [ymaxIter ; max(gpr_mdl_SUR.Eval(xNew))];
yminIter = abs([yminIter ;min(gpr_mdl_SUR.Eval(xNew))]);
subplot(2,1,1),plot(iter,[ymaxTrue ymaxIter],'--o','linewidth',2),grid on,
ylabel('Response')
xlabel('Number of HF samples')
subplot(2,1,2),plot(iter,[yminTrue yminIter],'--o','linewidth',2),grid on,
ylabel('-Response')
xlabel('Number of HF samples')


plottingGPs
end