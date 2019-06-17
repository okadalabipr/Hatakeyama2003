figure('Renderer', 'painters', 'Position', [0 0 900 600]);

subplot(3,3,1);
hold on;
plot(t,sim_RP(:,1),'-','LineWidth',2);
plot(t,sim_RP(:,2),'-.','LineWidth',2);
plot(t,sim_RP(:,3),'--','LineWidth',2);
set(gca, 'LineWidth', 1, 'FontSize', 8, 'FontName', 'Arial');
set(gca, 'Xtick', [0 200 400 600 800 1000 1200 1400 1600 1800],...
    'Xticklabel', [0 200 400 600 800 1000 1200 1400 1600 1800]);
xlim([0 1800]);
xlabel('TIME (s)');
ylim([0 90]);
ylabel('RP (%)');
legend('HRG=10nM','HRG=1nM','HRG=0.1nM');
hold off;
box on;

subplot(3,3,2);
hold on;
plot(t,sim_ShP(:,1),'-','LineWidth',2);
plot(t,sim_ShP(:,2),'-.','LineWidth',2);
plot(t,sim_ShP(:,3),'--','LineWidth',2);
set(gca, 'LineWidth', 1, 'FontSize', 8, 'FontName', 'Arial');
set(gca, 'Xtick', [0 200 400 600 800 1000 1200 1400 1600 1800],...
    'Xticklabel', [0 200 400 600 800 1000 1200 1400 1600 1800]);
xlim([0 1800]);
xlabel('TIME (s)');
ylim([0 100]);
ylabel('ShP (%)');
hold off;
box on;

subplot(3,3,3);
hold on;
plot(t,sim_PI3K_act(:,1),'-','LineWidth',2);
plot(t,sim_PI3K_act(:,2),'-.','LineWidth',2);
plot(t,sim_PI3K_act(:,3),'--','LineWidth',2);
set(gca, 'LineWidth', 1, 'FontSize', 8, 'FontName', 'Arial');
set(gca, 'Xtick', [0 200 400 600 800 1000 1200 1400 1600 1800],...
    'Xticklabel', [0 200 400 600 800 1000 1200 1400 1600 1800]);
xlim([0 1800]);
xlabel('TIME (s)');
ylim([0 70]);
ylabel('PI3K* (%)');
hold off;
box on;

subplot(3,3,4);
hold on;
plot(t,sim_Raf_act(:,1),'-','LineWidth',2);
plot(t,sim_Raf_act(:,2),'-.','LineWidth',2);
plot(t,sim_Raf_act(:,3),'--','LineWidth',2);
set(gca, 'LineWidth', 1, 'FontSize', 8, 'FontName', 'Arial');
set(gca, 'Xtick', [0 200 400 600 800 1000 1200 1400 1600 1800],...
    'Xticklabel', [0 200 400 600 800 1000 1200 1400 1600 1800]);
xlim([0 1800]);
xlabel('TIME (s)');
ylim([0 4]);
ylabel('Raf* (%)');
hold off;
box on;

subplot(3,3,5);
hold on;
plot(t,sim_MEKPP(:,1),'-','LineWidth',2);
plot(t,sim_MEKPP(:,2),'-.','LineWidth',2);
plot(t,sim_MEKPP(:,3),'--','LineWidth',2);
set(gca, 'LineWidth', 1, 'FontSize', 8, 'FontName', 'Arial');
set(gca, 'Xtick', [0 200 400 600 800 1000 1200 1400 1600 1800],...
    'Xticklabel', [0 200 400 600 800 1000 1200 1400 1600 1800]);
xlim([0 1800]);
xlabel('TIME (s)');
ylim([0 90]);
ylabel('MEKPP (%)');
hold off;
box on;

subplot(3,3,6);
hold on;
plot(t,sim_ERKPP(:,1),'-','LineWidth',2);
plot(t,sim_ERKPP(:,2),'-.','LineWidth',2);
plot(t,sim_ERKPP(:,3),'--','LineWidth',2);
set(gca, 'LineWidth', 1, 'FontSize', 8, 'FontName', 'Arial');
set(gca, 'Xtick', [0 200 400 600 800 1000 1200 1400 1600 1800],...
    'Xticklabel', [0 200 400 600 800 1000 1200 1400 1600 1800]);
xlim([0 1800]);
xlabel('TIME (s)');
ylim([0 90]);
ylabel('ERKPP (%)');
hold off;
box on;

subplot(3,3,8);
hold on;
plot(t,sim_Akt_PI_PP(:,1),'-','LineWidth',2);
plot(t,sim_Akt_PI_PP(:,2),'-.','LineWidth',2);
plot(t,sim_Akt_PI_PP(:,3),'--','LineWidth',2);
set(gca, 'LineWidth', 1, 'FontSize', 8, 'FontName', 'Arial');
set(gca, 'Xtick', [0 200 400 600 800 1000 1200 1400 1600 1800],...
    'Xticklabel', [0 200 400 600 800 1000 1200 1400 1600 1800]);
xlim([0 1800]);

xlabel('TIME (s)');
ylim([0 50]);
ylabel('AktPP (%)');
hold off;
box on;