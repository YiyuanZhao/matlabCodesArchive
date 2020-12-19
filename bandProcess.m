load('D:\OneDrive - tongji.edu.cn\vscode_workspace\DESKTOP-DQVLUVG\VScode_WorkSpace\model\data\spectrum_orig.dat');
fermiLevel = -2.4474445;
plot(kPointsMesh(4, :), Hkin - fermiLevel, 'r');
hold on;
plot(kPointsMesh(4, end)/spectrum_orig(end, 1)*spectrum_orig(:, 1), spectrum_orig(:, 2) - fermiLevel, 'b');
hold off;
legend([{'Model'}, {'Fortran Code'}]);