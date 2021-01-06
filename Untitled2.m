for i = 2: 11
    hoppingMatrixTemp{i} = hoppingMatrixOrigin{i};
    len = length(hoppingMatrixTemp{i}(:, 1));
    hoppingOrder = i - 1;
    hoppingMatrixTemp{i} = [hoppingOrder * ones(len, 1), hoppingMatrixTemp{i}, zeros(len, 1), ones(len, 2)];
    if i == 2
        out = hoppingMatrixTemp{i};
    else
        out = [out ;hoppingMatrixTemp{i}];
    end
%     writecell({[hoppingOrder, len]}, 'hoppingPara.txt', 'WriteMode', 'append', 'Delimiter', 'tab');
%     writemcell(hoppingMatrixTemp(i), 'hoppingPara.txt', 'WriteMode', 'append', 'Delimiter', 'tab');
end

load('D:\OneDrive - tongji.edu.cn\vscode_workspace\DESKTOP-DQVLUVG\VScode_WorkSpace\model\data\spectrum_orig.dat');
figure;
plot(spectrum_orig(:, 1)*kPointsMesh(4, end)/spectrum_orig(end, 1), spectrum_orig(:, 2))
load('D:\OneDrive - tongji.edu.cn\vscode_workspace\DESKTOP-DQVLUVG\VScode_WorkSpace\model\data\wanner90_band.dat');
hold on
plot(wanner90_band(:, 1), wanner90_band(:, 2))
hold off
kPointsMesh(2, 2)
legend