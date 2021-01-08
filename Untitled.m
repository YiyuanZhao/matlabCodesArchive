clear variables
%% Data Preprocessing
% hoppingInfo = load('D:\OneDrive - tongji.edu.cn\vscode_workspace\DESKTOP-DQVLUVG\VScode_WorkSpace\model\data\grapProcessed.dat');
% degeneracy = load('degeneracy.mat');
hoppingInfo = load('C:\Users\SCES\AppData\Roaming\MathWorks\MATLAB Add-Ons\Apps\VASP_Plotter - 1\hoppingParameterProcessed.dat');
degeneracy = load('C:\Users\SCES\AppData\Roaming\MathWorks\MATLAB Add-Ons\Apps\VASP_Plotter - 1\degeneracy.dat');
hoppingInfo(:, 8) = degeneracy';

% Danger zone degeneracy correction
hoppingInfo(:, 6) = hoppingInfo(:, 6) ./ hoppingInfo(:, 8);
% Danger zone ended

hoppingParameter = hoppingInfo(:, 6);
a1 = [3.3298712493556994    0.0000000000000000    0.0000000000000000];
a2 = [-1.6649356246778497    2.8837530932734619   -0.0000000000000000];
a3 = [-0.0000000000000001   -0.0000000000000000   23.1180000394582805];
transMatA = [a1;a2;a3];
V = dot(a1, cross(a2, a3));
b1 = 2*pi*cross(a2,a3)/V;
b2 = 2*pi*cross(a3,a1)/V;
b3 = 2*pi*cross(a1,a2)/V;
transMatB = [b1;b2;b3];
%% Divides the hopping terms by hopping parameter difference
[~ , indexSorted] = sort(hoppingParameter);
hoppingInfoSorted = hoppingInfo(indexSorted, :);
[~, indexUnique] = uniquetol(hoppingInfoSorted(:, 6), 1e-6);

hoppingParameterOrderLoop = cell(1, length(indexUnique));
hoppingMatrixOrderLoop = cell(1, length(indexUnique));
for loopIndex = 1: length(indexUnique)
    try
        hoppingParameterOrderLoop{loopIndex} = hoppingInfoSorted(indexUnique(loopIndex), 6);
        hoppingMatrixOrderLoop{loopIndex} = hoppingInfoSorted(indexUnique(loopIndex): indexUnique(loopIndex + 1) - 1, 1: 3);
    catch
        hoppingParameterOrderLoop{loopIndex} = hoppingInfoSorted(end, 6);
        hoppingMatrixOrderLoop{loopIndex} = hoppingInfoSorted(indexUnique(loopIndex): end, 1: 3);
    end
end

%% Read K-Mesh
%  load('savedFigureData.mat');
%  selectionIndex = ~isnan(exportData.bandStructure_3D.bandeigenvalueSelected);
%  kPointsMesh = exportData.hubbardModel.kpointsRecipMeshSelected(1:3, selectionIndex);

%% OR Generate K-Mesh \Gamma --> M --> K --> \Gamma
kpoints{1} = [0 0 0]* transMatB;
kpoints{2} = [0 0.5 0]* transMatB;
kpoints{3} = [-0.3333 0.6667 0.0]* transMatB;
kpoints{4} = [0 0 0]* transMatB;
Ntot = 51;
for loopIndex = 1: length(kpoints) - 1
    kPointsMesh(Ntot*(loopIndex - 1) + 1: Ntot* loopIndex, 1: 3) = [linspace(kpoints{loopIndex}(1), kpoints{loopIndex + 1}(1), Ntot)', linspace(kpoints{loopIndex}(2), kpoints{loopIndex + 1}(2), Ntot)', ...
        linspace(kpoints{loopIndex}(3), kpoints{loopIndex + 1}(3), Ntot)'];
end
numIdx = Ntot: Ntot:(length(kpoints) - 2)*Ntot;
kPointsMesh(numIdx, :) = [];
kPathDistance = 0;
kPointsMesh(1, 4) = 0;
for loopIndex = 2: length(kPointsMesh)
    sqrtPart(1) = (kPointsMesh(loopIndex, 1) - kPointsMesh(loopIndex - 1, 1));
    sqrtPart(2) = (kPointsMesh(loopIndex, 2) - kPointsMesh(loopIndex - 1, 2));
    sqrtPart(3) = (kPointsMesh(loopIndex, 3) - kPointsMesh(loopIndex - 1, 3));
    kPathDistance = kPathDistance + sqrt(sqrtPart(1)^2 + sqrtPart(2)^2 + sqrtPart(3)^2 );
    kPointsMesh(loopIndex, 4) = kPathDistance;
end
kPointsMesh = kPointsMesh';
%% Calculate H_kin
Hkin = zeros(1, length(kPointsMesh));
for loopIndex = 1:length(kPointsMesh)
    Hkin(loopIndex) = calculateKineticHamiltonian(transMatA, kPointsMesh(1:3, loopIndex)', hoppingParameterOrderLoop, hoppingMatrixOrderLoop);
end
%  scatter3(kPointsMesh(1, :), kPointsMesh(2, :), Hkin, 1);
scatter(kPointsMesh(4, :), Hkin, 1);

%% Compare the calculation result to the wannier fitting
load('D:\OneDrive - tongji.edu.cn\vscode_workspace\DESKTOP-DQVLUVG\VScode_WorkSpace\model\data\wanner90_band.dat');
hold on;
plot(wanner90_band(:, 1), wanner90_band(:, 2));
hold off;
atomPosition = hoppingInfoSorted(:, 1).* a1 + hoppingInfoSorted(:, 2).* a2 + hoppingInfoSorted(:, 3).* a3;
generacyOne = hoppingInfoSorted(:, 8) == 1;
generacyTwo = hoppingInfoSorted(:, 8) == 2;
fig = figure();
scatter(atomPosition(generacyOne, 1), atomPosition(generacyOne, 2), 'filled', 'r');
hold on
scatter(atomPosition(generacyTwo, 1), atomPosition(generacyTwo, 2), 'filled', 'b');
xlim([-40 40]);
%% H_{kin} = exp{i \vec{k} \cdot \vec{r}}
function Hkin = calculateKineticHamiltonian(transMat, kpointsRecip, hoppingParameter, hoppingMatrixProcessed)
Hkin = 0;
for hoppingOrder = 1: length(hoppingParameter)
    sum = 0;
    hoppingLength = length(hoppingMatrixProcessed{hoppingOrder}(:, 1));
    a1part = cell(1, hoppingLength);
    a2part = cell(1, hoppingLength);
    a3part = cell(1, hoppingLength);
    dotpart = zeros(1, hoppingLength);
    for numIdx = 1: hoppingLength
        a1part{numIdx} = hoppingMatrixProcessed{hoppingOrder}(numIdx, 1).*transMat(1, :);
        a2part{numIdx} = hoppingMatrixProcessed{hoppingOrder}(numIdx, 2).*transMat(2, :);
        a3part{numIdx} = hoppingMatrixProcessed{hoppingOrder}(numIdx, 3).*transMat(3, :);
        dotpart(numIdx) = dot(kpointsRecip, a1part{numIdx} + a2part{numIdx} + a3part{numIdx});
        sum = sum + exp(1i*dotpart(numIdx));
    end
    Hkin = Hkin + hoppingParameter{hoppingOrder}*sum;
end
end