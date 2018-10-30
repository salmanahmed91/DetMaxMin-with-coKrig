%close all; 
f1 = figure(767)
%f1.Name = num2str(MM(i)); 
 
xStep = 0.875;
tauP  = 28;
posX = [[0:xStep:tauP]./tauP]';
plot(posX,FTX_det_3D(:,M),'linewidth',2),hold all; 
plot(posX,FTX_det_2D(1:33,151+M),'linewidth',2),hold all; 
plot(gpr_mdl_2D.X,gpr_mdl_2D.Y(1:end),'o','linewidth',1.5,'MarkerSize',7);
plot(gpr_mdl_SUR.X,gpr_mdl_SUR.Y(1:end),'^','linewidth',1.5,'MarkerSize',7);

posXfine = [0:0.001:1]';

gprYC = gpr_mdl_2D.Eval(posXfine);
plot(posXfine, [gprYC ],'linewidth',2);

gprYE = gpr_mdl_YE.Eval(posXfine);
plot(posXfine, [gprYE ],'--','linewidth',1);


gprYE_coK = gpr_mdl_SUR.Eval(posXfine);
plot(posXfine, [gprYE_coK ],'linewidth',2);

gprE  = gpr_mdl_ERROR.Eval(posXfine);
gprRecon = gpr_mdl_SUR.rho.*gpr_mdl_2D.Eval(posXfine)+gpr_mdl_ERROR.Eval(posXfine);
plot(posXfine, [gprE gprRecon],'linewidth',2);


grid on;xlabel('x'),ylabel('Response'),
legend('show')
legend('HF True','LF True','LF samples','HF samples','gpr Yc','gpr Ye',...
        'coKrig','gpr error','gpr Recon','location','eastoutside','orientation','vertical')
hold off;

% figName = strcat('case',num2str(MM(i)),'.fig');
% saveas(gcf,figName);