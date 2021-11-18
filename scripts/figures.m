%diagrams

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#D95319';
color1 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#0072BD';
color2 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#77AC30';
color3 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
color4 = [164 167 233]/255;
color5 = [182 180 200]/255;
color6 = [140 200 147]/255; 

%Basic system
f1 = figure('Name', 'Basic system');
f1.Position(1:2) = [120 550];
plot(doo, Hdi_sim, '.', doo, Hdi, 'Color', color1, 'LineWidth', 1.5, 'MarkerSize', 13); 
hold on;
plot(doo, Qdi_sim, '.', doo, Qdi, 'Color', color2, 'LineWidth', 2, 'MarkerSize', 13);
plot(doo, Hdi_sim.*Qdi_sim, '.', doo, Hdi.*Qdi, 'Color', color3, 'LineWidth', 2, 'MarkerSize', 13);
grid on;
L1(1) = plot(nan, nan, 'Color', color1);
L1(2) = plot(nan, nan, 'Color', color2);
L1(3) = plot(nan, nan, 'Color', color3);
legend(L1,{'Connection probability', 'Capture probability', 'Coverage probability'}, 'Location', 'Southwest')
set(gca,'FontSize',12)
xlabel('Distance [km]','fontsize',16);
ylabel('Probabilities','fontsize',16)
set(gca,'YLim',[0 1])
set(gca,'YTick',(0:0.2:1))
legend('boxoff')

%Capture probability
f2 = figure('Name', 'Capture probability with SIC');
f2.Position(1:2) = [710 550];
plot(doo, Qdi_sim, '.', doo, Qdi, 'Color', color2, 'LineWidth', 2, 'MarkerSize', 13);
hold on;
plot(doo, Qdk_1sim+Qdk_2sim.*Qsic_sim, '.', doo, Qdk_1.*Pr1+Pr_bigger_than_one.^2.*Qdi.*Qdk_2, 'Color', color4, 'LineWidth', 1.5, 'MarkerSize', 13); 
plot(doo, Qdi_sim+Qdk_1sim+Qdk_2sim.*Qsic_sim, '.', doo, Qdi+Qdk_1.*Pr1+Pr_bigger_than_one.^2.*Qdi.*Qdk_2, 'Color', color5, 'LineWidth', 2, 'MarkerSize', 13);
grid on;
L(1) = plot(nan, nan, 'Color', color2);
L(2) = plot(nan, nan, 'Color', color4);
L(3) = plot(nan, nan, 'Color', color5);
legend(L,{'Capture probability', 'Successive SIC Probability', 'Total Capture Probability'}, 'Location', 'Northeast')
set(gca,'FontSize',12)
xlabel('Distance [km]','fontsize',16);
ylabel('Probabilities','fontsize',16)
legend('boxoff')

%Coverage probabilities
f3 = figure('Name', 'Coverage probabilities');
f3.Position(1:2) = [1300 550];
plot(doo, Hdi_sim.*Qdi_sim+Hdi_sim.*(Qdk_1sim+Qdk_2sim.*Qsic_sim), '.', doo, Hdi.*Qdi+Hdi.*(Qdk_2.*Pr_bigger_than_one.^2.*Qdi+Qdk_1.*Pr1), 'Color', color6, 'LineWidth', 2, 'MarkerSize', 13)
hold on;
plot(doo, Hdi_sim.*Qdi_sim, '.', doo, Hdi.*Qdi, 'Color', color2, 'LineWidth', 2, 'MarkerSize', 13);

grid on;
L2(1) = plot(nan, nan, 'Color', color6);
L2(2) = plot(nan, nan, 'Color', color3);

legend(L2,{'Coverage Probability after SIC', 'Coverage Probability before SIC'}, 'Location', 'Northeast')
set(gca,'FontSize',12)
xlabel('Distance [km]','fontsize',16);
ylabel('Probabilities','fontsize',16)
set(gca,'YLim',[0 1])
set(gca,'YTick',(0:0.2:1))
legend('boxoff')

