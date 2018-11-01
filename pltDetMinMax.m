


ymin = min(gpr_mdl_SUR.Y);
f1 = figure(450);
xNew = [gpr_mdl_SUR.lbX:0.01:gpr_mdl_SUR.ubX]';
[Fpred,RMSE] = gpr_mdl_SUR.Eval(xNew);
EI_xMax  = EI(ymax,xNew,gpr_mdl_SUR,'max');
EI_xMin  = EI(ymin,xNew,gpr_mdl_SUR,'min');
subplot(211),plot(xNew,[gpr_mdl_YE.Eval(xNew)],'-','linewidth',1,'DisplayName','HF True');hold all
subplot(211),plot(xNew,[Fpred] ,'linewidth',2,'DisplayName','$\hat{y}_e(\mathbf{x})$');grid on
subplot(211),plot(xNew,[Fpred+2.*RMSE],'-.','linewidth',1,'DisplayName','$\hat{y}_e(\mathbf{x})+2\hat{s}(\mathbf{x})$')
subplot(211),plot(xNew,[Fpred-2.*RMSE],'-.','linewidth',1,'DisplayName','$\hat{y}_e(\mathbf{x})-2\hat{s}(\mathbf{x})$')
subplot(211),
plot(gpr_mdl_SUR.X,gpr_mdl_SUR.Y(1:end),...
    's','linewidth',1,...
    'MarkerSize',6,...
    'DisplayName','$\mathbf{{y_e}}$',...
    'MarkerFaceColor','y',...
    'MarkerEdgeColor','r');%xlabel('x'),
ylabel('Response');

%MarkerTxt(gpr_mdl_SUR.X,gpr_mdl_SUR.Y,1:length(gpr_mdl_SUR.X));
%legend('show'),

hold off;


subplot(212),plot(xNew,[EI_xMax],'linewidth',2,'DisplayName','EI max'),grid on,xlabel('x'),ylabel('EI(x)');hold all
subplot(212),plot(xNew,[EI_xMin],'linewidth',2,'DisplayName','EI min'),grid on,xlabel('x'),ylabel('EI(x)');
ylim([-1e-9 inf]);
%legend('show');hold all
%hldng = legend('EI max','EI min');
hold off;
reScalingPlots(f1,14);
