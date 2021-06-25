% This is a collection of small scripts used in daily life.
% clear variables;
a1 = [3.3291   -0.0000    0.0000];
a2 = [-1.6645    2.8831   -0.0000];
a3 = [0.0000   -0.0000   23.1180];
transMatA = [a1;a2;a3];
V = dot(a1, cross(a2, a3));
b1 = 2*pi*cross(a2,a3)/V;
b2 = 2*pi*cross(a3,a1)/V;
b3 = 2*pi*cross(a1,a2)/V;
transMatB = [b1;b2;b3];
% load('savedFigureData.mat');
highSymmetryPoints = [0, 0, 0; 0, 0.5, 0; -0.33333, 0.66667, 0.0; 0, 0, 0];
selectionIndex = ~isnan(exportData.bandStructure_3D.bandeigenvalueSelected);
kPointsMesh = exportData.hubbardModel.kpointsRecipMeshSelected(1:3, selectionIndex);
targetHkin = exportData.hubbardModel.bandeigenvalueSelected(selectionIndex);
highSymmetryPoints = highSymmetryPoints*transMatB;
pathSelection = cell(1, length(highSymmetryPoints) - 1);
resolution = 1.6e-2;
for numIdx = 1: length(highSymmetryPoints) - 1
    if abs(highSymmetryPoints(numIdx + 1, 1) - highSymmetryPoints(numIdx, 1)) < 1e-3
        pathSelection{numIdx} = abs(kPointsMesh(1, :) - highSymmetryPoints(numIdx, 1)) < resolution ...
            & kPointsMesh(2, :) >= min([highSymmetryPoints(numIdx, 2), highSymmetryPoints(numIdx + 1, 2)]) ...
            & kPointsMesh(2, :) <= max([highSymmetryPoints(numIdx, 2), highSymmetryPoints(numIdx + 1, 2)]);
    elseif abs(highSymmetryPoints(numIdx + 1, 2) - highSymmetryPoints(numIdx, 2)) < 1e-3
        pathSelection{numIdx} = abs(kPointsMesh(2, :) - highSymmetryPoints(numIdx, 2)) < resolution ...
            & kPointsMesh(1, :) >= min([highSymmetryPoints(numIdx, 1), highSymmetryPoints(numIdx + 1, 1)]) ...
            & kPointsMesh(1, :) <= max([highSymmetryPoints(numIdx, 1), highSymmetryPoints(numIdx + 1, 1)]);       
    else
        for loopIndex = 1: length(kPointsMesh(1, :))
            pathSelection{numIdx}(loopIndex) = abs((highSymmetryPoints(numIdx + 1, 1) - highSymmetryPoints(numIdx, 1)) ...
                *(kPointsMesh(2, loopIndex) - highSymmetryPoints(numIdx, 2)) - ...
                (highSymmetryPoints(numIdx + 1, 2) - highSymmetryPoints(numIdx, 2))* ...
                (kPointsMesh(1, loopIndex) - highSymmetryPoints(numIdx, 1))) < resolution;
        end
        pathSelection{numIdx} = pathSelection{numIdx} ...
            & kPointsMesh(2, :) >= min([highSymmetryPoints(numIdx, 2), highSymmetryPoints(numIdx + 1, 2)]) ...
            & kPointsMesh(2, :) <= max([highSymmetryPoints(numIdx, 2), highSymmetryPoints(numIdx + 1, 2)]) ...
            & kPointsMesh(1, :) >= min([highSymmetryPoints(numIdx, 1), highSymmetryPoints(numIdx + 1, 1)]) ...
            & kPointsMesh(1,: ) <= max([highSymmetryPoints(numIdx, 1), highSymmetryPoints(numIdx + 1, 1)]);
    end
end
pathSelectionMix = false(1, length(kPointsMesh(1, :)));
for numIdx = 1: length(pathSelection)
    pathSelectionMix = pathSelectionMix | pathSelection{numIdx};
end
pathSelectionMix = find(pathSelectionMix == 1);
scatter(kPointsMesh(1, pathSelectionMix), kPointsMesh(2, pathSelectionMix), 2);

weight = ones(1, length(kPointsMesh(1, :)));
weight(pathSelectionMix) = 2;
% weightloop = zeros(1, length(kPointsMesh(1, :)));
% kxMax = max(kPointsMesh(1, :));
% for numIdx = 1: length(pathSelectionMix)
%     mu1 = kPointsMesh(1, pathSelectionMix(numIdx));
%     mu2 = kPointsMesh(2, pathSelectionMix(numIdx));
%     sigma = 0.1995;
%     for loopIndex = 1: length(kPointsMesh(1, :))
%         x = kPointsMesh(1, loopIndex);
%         y = kPointsMesh(2, loopIndex);
%         distance = sqrt(kPointsMesh(1, loopIndex)^2 + kPointsMesh(2, loopIndex)^2);
%         weightloop(loopIndex) = 1/(2*pi()*sigma^2)*exp(-1/2*((x-mu1)^2/sigma^2 + (y-mu2)^2/sigma^2));
%     end
%     weight = weight + weightloop;
% end
scatter3(kPointsMesh(1, :), kPointsMesh(2, :), weight, 2);

%% Calculate methods (energy Level)
maxWeight = 2;
lifeWidth = 0.70316;
weight = calculateWeightProportialToEnergy(maxWeight, lifeWidth, targetHkin, 0);
figure;
plot(kPointsMesh(4, :), weight);
%% Calculate lattice vectors
clear variables
a1 = [6.5085482243     -0.0000          0.0000];
a2 = [-3.2542741121    5.6365681040    -0.0000];
a3 = [0.0000   -0.0000   21.118000039458281];
transMatA = [a1;a2;a3];
V = dot(a1, cross(a2, a3));
b1 = 2*pi*cross(a2,a3)/V;
b2 = 2*pi*cross(a3,a1)/V;
b3 = 2*pi*cross(a1,a2)/V;
transMatB = [b1;b2;b3];
%% Read bandstructure
clear variables
load("D:\OneDrive - tongji.edu.cn\BLG\nm\wannier90_band.dat");
load("D:\OneDrive - tongji.edu.cn\vscode_workspace\DESKTOP-DQVLUVG\VScode_WorkSpace\model\data\spectrum_orig_modified.dat");
load("D:\OneDrive - tongji.edu.cn\vscode_workspace\DESKTOP-DQVLUVG\VScode_WorkSpace\model\data\spectrum_orig.dat");
spectrum_orig_modified(:, 1) = spectrum_orig_modified(:, 1)*(1.5227/max(spectrum_orig_modified(:, 1)));
spectrum_orig(:, 1) = spectrum_orig(:, 1)*(1.5227/max(spectrum_orig(:, 1)));
% scatter(wannier90_band(:, 1), wannier90_band(:, 2), 1);
hold on;
% scatter(spectrum_orig(:, 1), spectrum_orig(:, 2), 1);
plot(spectrum_orig_modified(:, 1), spectrum_orig_modified(:, 2) + 0.51633);
plot(spectrum_orig(:, 1), spectrum_orig(:, 2) + 1.842);
hold off
%% Recognize Circle
clear variables;
latticePara = 13.3194849974227978;
% supercell = imread('C:\Users\SCES\Desktop\supercell.jpg');
supercell = imread('C:\Users\SCES\Desktop\supercellProcess2.jpg');
featureOrigin = imread('C:\Users\SCES\Desktop\feature.jpg');
supercellGray = rgb2gray(supercell);
featureOriginGray = rgb2gray(featureOrigin);
% points1 = detectSURFFeatures(featureOriginGray);
% [features1, valid_points1] = extractFeatures(featureOriginGray, points1);
% points2 = detectSURFFeatures(supercellGray);
% [features2, valid_points2] = extractFeatures(supercellGray, points2);
% indexPairs = matchFeatures(features1,features2);
% matchedPoints1 = valid_points1(indexPairs(:,1),:);
% matchedPoints2 = valid_points2(indexPairs(:,2),:);
% figure;
% showMatchedFeatures(featureOriginGray,supercellGray,matchedPoints1,matchedPoints2);

bw = imbinarize(supercellGray, 0.64);
[centers, radii, metric] = imfindcircles(bw, [30 100], 'Sensitivity', 0.85, ...
    'EdgeThreshold', [], 'Method', 'PhaseCode', 'ObjectPolarity', 'dark');
imshow(bw);
viscircles(centers, radii, 'EdgeColor', 'b');
% imshow(featureOrigin);
% hold on
% plot(features.selectStrongest(10));
% hold off;

% tanTheta = (centers(2,1) - centers(9,1))/(centers(2,2) - centers(9,2));
% alpha = sqrt(3)/(3*tanTheta);
% centers(:, 1) = centers(:, 1)*alpha;
% 
% centers = centers - centers(2, :);
% scalingPara = latticePara / norm(centers(6, :));
% centers(6, :) = [];
% centers = centers * scalingPara;
% centers(:, 2) = -centers(: ,2);

tanTheta = (centers(4,1) - centers(8,1))/(centers(4,2) - centers(8,2));
alpha = sqrt(3)/(3*tanTheta);
centers(:, 1) = centers(:, 1)*alpha;

centers = centers - centers(4, :);
scalingPara = latticePara / norm(centers(3, :));
centers(8, :) = [];
centers(3, :) = [];
centers = centers * scalingPara;
centers(:, 2) = -centers(: ,2);
figure;
scatter(centers(:, 1), centers(:, 2));
%% Get magnetization
clear variables
clc;
mag = load('D:\tmp\mag.dat');
workMode = 'fm';
atomNumber = length(mag(:, 1))/3;
switch workMode
    case 'fm'
        magmom_V  = mean(mag(1: atomNumber, 5));
        magmom_Se = mean(mag(atomNumber + 1: end, 5));
        var_V     = var(mag(1: atomNumber, 5));
    case 'afm'
        magmom_V  = mean(abs(mag(1: atomNumber, 5)));
        magmom_Se = mean(abs(mag(atomNumber + 1: end, 5)));
        var_V     = var(abs(mag(1: atomNumber, 5)));
end
% output = [magmom_V, magmom_Se, var_V];
% clipboard('copy', output);
%% Calculate DOS from H_k
clear variables
numberOfKpoints = 17091;
% fermiLevel = -2.08436568;
fermiLevel = 0;
gamma = 0.001;
Hk_origin = load('D:\OneDrive - tongji.edu.cn\vscode_workspace\DESKTOP-DQVLUVG\VScode_WorkSpace\model\data\HK.dat');
index = 1;
for kNumIdx = 1: numberOfKpoints
    for iorbNumidx = 1: 28
        for jorbNumIdx = 1: 28
            Hk(jorbNumIdx, iorbNumidx, kNumIdx) = Hk_origin(index, 1) + 1i*Hk_origin(index, 2);
            index = index + 1;
        end
    end
end

% omega = linspace(-13, 16, 10000);
omega = linspace(-3, 3, 10000);
GreenFunc = zeros(28, 28, 10000);
parfor omegaNumIdx = 1: length(omega)
    for kNumIdx = 1: numberOfKpoints
        GreenFunc(:, :, omegaNumIdx) = GreenFunc(:, :, omegaNumIdx) + ...
            inv(diag(ones(1, length(Hk(:, :, kNumIdx))), 0)*(omega(omegaNumIdx) + fermiLevel + 1i*gamma) - Hk(:, :, kNumIdx));
    end
    GreenFunc(:, :, omegaNumIdx) = GreenFunc(:, :, omegaNumIdx)/numberOfKpoints;
end

for omegaNumIdx = 1: length(omega)
   dos(omegaNumIdx) = -1/pi()* imag(trace(GreenFunc(:, :, omegaNumIdx)));
end
plot(omega, dos);
dosPlot = [omega', dos'];
save('dos.txt', 'dosPlot', '-ascii','-double');
%% Determine Fermi Level of BLG
numIdxMin = 1;
numIdxMax = length(dos);

while true
    numIdxMid = floor((numIdxMin + numIdxMax)/2);
    n1 = sum(dos(1: numIdxMid)*(omega(2) - omega(1)));
    n2 = sum(dos*(omega(2) - omega(1)));
    if n2/n1 > 2
        numIdxMin = numIdxMid;
    else
        numIdxMax = numIdxMid;
    end
    
    if numIdxMax - numIdxMin == 1
        break
    end
end
omega(numIdxMid)
% plot(omega, dos);
%% Validate hopping paramter
% import("wannier90hrmodified.dat");
clear hoppingParameter;
checkOrder = 33;
hoppingParameter = zeros(length(output.hopping{1, checkOrder}), 1);
for numIdx = 1: length(output.hopping{1, checkOrder})
    idx1 = output.hopping{1, checkOrder}(numIdx, 1);
    idx2 = output.hopping{1, checkOrder}(numIdx, 2);
    try 
        hoppingParameter(numIdx, 1) = wannier90hrmodified(wannier90hrmodified(:,1) == idx1 & wannier90hrmodified(:, 2) == idx2, 6);
    catch
        hoppingParameter(numIdx, 1) = NaN;
    end
end

%% Construct primitive cell

%% Math
clear variables
a = linspace(4/sqrt(3), 8/sqrt(3)-0.0001,100);
x_A = (a + sqrt(a.^2 - 4*(a.^2 - 16)))/4;
x_lb = 4/sqrt(3);
x_ub = 10;
% for i = 1: length(a)
%     err = Inf;
%     while (true)
%         if err < 1e-8
%             dist(i) = dist_Orig;
%             break
%         end
%         x_mid = (x_lb + x_ub)/2;
%         [dist_lb, ~] = calcDist(x_A(i), a(i), x_lb);
%         [dist_ub, ~] = calcDist(x_A(i), a(i), x_ub);
%         [dist_mid,dist_Orig] = calcDist(x_A(i), a(i), x_mid);
%         if dist_lb*dist_mid > 0
%             x_lb = x_mid;
%         else
%             x_ub = x_mid;
%         end
%         err = abs(dist_lb - dist_ub);
%         disp(err);
%     end 
% end
fminuncOptions = optimoptions('fminunc', "Algorithm","quasi-Newton", "MaxIterations", 1000,...
    "Display","iter-detailed",'FunctionTolerance',1e-9,...
    "MaxFunctionEvaluations", 40000, "UseParallel",true, ...
    "FiniteDifferenceType", 'central', 'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-7);
% fminsearchOptions = optimset('Display', 'iter', 'FunValCheck', 'on', 'PlotFcns', @optimplotfval);
% [x,fval,exitflag,output] = fminsearch(f, initialHoppingParameter, fminsearchOptions);
for i = 1: length(a)
    x = a(i) + 1;
    f = @(x) calcDist(x_A(i), a(i), x);
    [x_D(i),fval,exitflag,output,~,~] = fminunc(f, x, fminuncOptions);
    y(i) = -(x_A(i) - a(i))*(x_D(i) - x_A(i))/(sqrt(3)*x_A(i)) + sqrt(3)*x_A(i);
    dist_Orig(i) = sqrt(x_D(i)^2 + y(i)^2);
end
plot(a, dist_Orig)
%% Wave test
clear variables
A = 1;
omega = 0.2;
lambda = 20;
t = 0: 0.1: 1000;
x = linspace(0, 100);
% an = animatedline;
y = 2*A*cos(2*pi()*x/lambda)*cos(omega*t(1));
fig = plot(x,y);
axis([0 100 -2 2]);
for i = 2: length(t)
    y = 2*A*cos(2*pi()*x/lambda)*cos(omega*t(i));
    fig.YData = y;
    drawnow 
end
%% Wave demo 2
clear variables;
A = 1;
waveLength = 50;
frequency = 2;
rho = 5000;
F = 500;
% lambda = sqrt(F/rho)/frequency;
lambda = 100;
x = linspace(0, 100, 500);
t = 0: 0.001: 1000;
phi_1 = pi()/2;% - 2*pi()/lambda.*x;
y = 2*A*cos(2*pi()/lambda.*x + phi_1)*cos(2*pi()*frequency*t(1));
fig = plot(x, y);
axis([0 100 -2 2]);
for i= 2: length(t)
    y = 2*A*cos(2*pi()/lambda.*x + phi_1)*cos(2*pi()*frequency*t(i));
    fig.YData = y;
    drawnow;
end
%% Construct primitive cell
% clear variables
vacuumLength = 20;
numberOfLayers = 2;
a = 3.3171310374808005;  % a-axies of the lattice
b = a;      % b-axies of the lattice
gamma = 120;% angle of <a, b>
sizeLattice = 1;  % Scale of the system (should be odd)

% Initialize of the variables
centerOrder = (sizeLattice + 1) ./ 2;
lattice.x = zeros(sizeLattice);
lattice.y = zeros(sizeLattice);
% Difference in x and y in the primitive cell
deltaA = [a, 0];
deltaB = [b*cosd(gamma), b*sind(gamma)];
% Set the zero of the axies
lattice.x(sizeLattice, 1) = 0;
lattice.y(sizeLattice, 1) = 0;
% Initialize the position of the center
lattice.x(centerOrder, centerOrder) = ((centerOrder - 1)*deltaA(1) + (sizeLattice - centerOrder)*deltaB(1));
lattice.y(centerOrder, centerOrder) = ((centerOrder - 1)*deltaA(2) + (sizeLattice - centerOrder)*deltaB(2));
% Calculate the distance from the center
for row = sizeLattice: -1: 1
    for column = 1: sizeLattice
        lattice.x(row, column) = ((column - 1)*deltaA(1) + (sizeLattice - row)*deltaB(1));
        lattice.y(row, column) = ((column - 1)*deltaA(2) + (sizeLattice - row)*deltaB(2));
        lattice.hoppingA(row, column) = column - centerOrder;
        lattice.hoppingB(row, column) = centerOrder - row;
    end
end
% Move the center to the axis zero
lattice.x = lattice.x - lattice.x(centerOrder, centerOrder);
lattice.y = lattice.y - lattice.y(centerOrder, centerOrder);
% Construct layers
d = 1.58106107700766;   %1.58106107700766;   % h-phase:1.59498653898723, t-phase:1.58106107700766
D = 3.52122839340594;   %2.809105554065346;  % h-phase:3.67176692979057, t-phase:3.12122839340594
                        % Origin: 3.12122839340594

if isfield(distance, 'mean_d')
    d = distance.mean_d;
    if isfield(distance, 'D_X2_X3')
        D = distance.D_X2_X3;
    end
end
                        
offsetSe1 = [deltaA', deltaB']*[2/3; 1/3];
% offsetSe1 = [deltaA', deltaB']*[1/3; 2/3];
offsetSe2 = [deltaA', deltaB']*[1/3; 2/3];
layer{1}.Se1.x = lattice.x + offsetSe1(1);
layer{1}.Se1.y = lattice.y + offsetSe1(2);
layer{1}.Se1.z = zeros(size(lattice.x));
layer{1}.V.x = lattice.x;
layer{1}.V.y = lattice.y;
layer{1}.V.z = zeros(size(lattice.x)) + d;
layer{1}.Se2.x = lattice.x + offsetSe2(1);
layer{1}.Se2.y = lattice.y + offsetSe2(2);
layer{1}.Se2.z = zeros(size(lattice.x)) + 2*d;

layer{2} = layer{1};
layer{2}.Se1.z = layer{1}.Se1.z + D + 2*d;
layer{2}.Se2.z = layer{1}.Se2.z + D + 2*d;
layer{2}.V.z = layer{1}.V.z + D + 2*d;

% Calculate POSCAR
switch numberOfLayers
    case 2
        primitiveCellC = vacuumLength + 4*d + D;
    case 1
        primitiveCellC = vacuumLength + 2*d;
end
VatomNumber = 0;
SeatomNumber = 0;
primitiveCellLattice = a;
primitiveCell.V.x = layer{1}.V.x;
primitiveCell.V.y = layer{1}.V.y;
primitiveCell.V.z = layer{1}.V.z;
primitiveCell.Se.x = cat(1, layer{1}.Se1.x, layer{1}.Se2.x);
primitiveCell.Se.y = cat(1, layer{1}.Se1.y, layer{1}.Se2.y);
primitiveCell.Se.z = cat(1, layer{1}.Se1.z, layer{1}.Se2.z);
for i = 1: numberOfLayers
    VatomNumber = VatomNumber + length(layer{i}.V.x);
    SeatomNumber = SeatomNumber + length(layer{i}.Se1.x) + length(layer{i}.Se2.x);
    if i >= 2
        primitiveCell.V.x = cat(1, primitiveCell.V.x, layer{i}.V.x);
        primitiveCell.V.y = cat(1, primitiveCell.V.y, layer{i}.V.y);
        primitiveCell.V.z = cat(1, primitiveCell.V.z, layer{i}.V.z);
        primitiveCell.Se.x = cat(1, primitiveCell.Se.x, layer{i}.Se1.x, layer{i}.Se2.x);
        primitiveCell.Se.y = cat(1, primitiveCell.Se.y, layer{i}.Se1.y, layer{i}.Se2.y);
        primitiveCell.Se.z = cat(1, primitiveCell.Se.z, layer{i}.Se1.z, layer{i}.Se2.z);
    end
end
% Write POSCAR
fileId = fopen("POSCAR", 'w');
fprintf(fileId, "Primitive Cell, D = %5.4f\n", D);
fprintf(fileId, "   1.00000000000000 \n");
fprintf(fileId, "%24.15f %24.15f %24.15f\n", [primitiveCellLattice, 0, 0]);
fprintf(fileId, "%24.15f %24.15f %24.15f\n", [-primitiveCellLattice/2, sqrt(3)/2*primitiveCellLattice, 0]);
fprintf(fileId, "%24.15f %24.15f %24.15f\n", [0, 0, primitiveCellC]);
fprintf(fileId, "   V    Se\n");
fprintf(fileId, "%6d %6d\n", [VatomNumber, SeatomNumber]);
fprintf(fileId, "Cartesian\n");
for numIdx = 1: VatomNumber
    fprintf(fileId, "% 18.16f % 21.16f % 21.16f\n", [primitiveCell.V.x(numIdx), primitiveCell.V.y(numIdx), primitiveCell.V.z(numIdx)]);
end

for numIdx = 1: SeatomNumber
    fprintf(fileId, "% 18.16f % 21.16f % 21.16f\n", [primitiveCell.Se.x(numIdx), primitiveCell.Se.y(numIdx), primitiveCell.Se.z(numIdx)]);
end
fclose(fileId);

%% Read POSCAR, obtain structure information
clear variables
filename = "C:\Users\SCES\Desktop\matlabCodes\POSCAR";
opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [1, Inf];
opts.Delimiter = ["\t", " "];
opts.VariableNames = ["col1", "col2", "col3", "col4"];
opts.VariableTypes = ["string", "string", "string", "string"];
opts = setvaropts(opts, [1, 2, 3, 4], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
POSCARtemp = readtable(filename, opts);
POSCARtemp = table2array(POSCARtemp);
POSCARtemp(ismissing(POSCARtemp)) = "";
for i = 1:length(POSCARtemp)
    tmp.numIdx = find( strlength(POSCARtemp(i, :)) );
    POSCAR.stringMat(i, 1: length(tmp.numIdx)) = POSCARtemp(i, tmp.numIdx);
end
POSCAR.stringMat(ismissing(POSCAR.stringMat)) = "";
POSCAR.element.name = POSCAR.stringMat(6, :);
POSCAR.element.name( POSCAR.element.name=='' )= [];
POSCAR.element.num = str2double (POSCAR.stringMat(7, :) );
POSCAR.element.num(isnan(POSCAR.element.num)) = [];
POSCAR.element.totalAtomNumber = sum(POSCAR.element.num);
tmp.totAtomNumber = sum(POSCAR.element.num);
clear opts i POSCARtemp

SCAL = str2double(POSCAR.stringMat(2,1));
a1 = str2double(POSCAR.stringMat(3, 1: 3))*SCAL;
a2 = str2double(POSCAR.stringMat(4, 1: 3))*SCAL;
a3 = str2double(POSCAR.stringMat(5, 1: 3))*SCAL;
POSCAR.atomPos = str2double(POSCAR.stringMat(9:(8 + POSCAR.element.totalAtomNumber), 1: 3));
% Remove the translation symmetry
POSCAR.atomPos(POSCAR.atomPos(:, 3) > 0.9, 3) = -(1 - POSCAR.atomPos(POSCAR.atomPos(:, 3) > 0.9, 3));
switch POSCAR.element.totalAtomNumber
    case 3
        numberOfLayers = 1;
        distance.latticeConst = a1(1);
        distance.d_V_X1 = dot((POSCAR.atomPos(1, :) - POSCAR.atomPos(2, :)), a3);
        distance.d_V_X2 = dot((POSCAR.atomPos(3, :) - POSCAR.atomPos(1, :)), a3);
        distance.mean_d = sum([distance.d_V_X1, distance.d_V_X2])/2;
    case 6
        numberOfLayers = 2;
        distance.latticeConst = a1(1);
        distance.d_V1_X1 = dot((POSCAR.atomPos(1, :) - POSCAR.atomPos(3, :)), a3);
        distance.d_V1_X2 = dot((POSCAR.atomPos(4, :) - POSCAR.atomPos(1, :)), a3);
        distance.D_X2_X3 = dot((POSCAR.atomPos(5, :) - POSCAR.atomPos(4, :)), a3);
        distance.d_V2_X3 = dot((POSCAR.atomPos(2, :) - POSCAR.atomPos(5, :)), a3);
        distance.d_V2_X4 = dot((POSCAR.atomPos(6, :) - POSCAR.atomPos(2, :)), a3);
        distance.mean_d = sum([distance.d_V1_X1, distance.d_V1_X2, distance.d_V2_X3, distance.d_V2_X4])/4;
end
disp(distance);
%% Temp Function
function weight = calculateWeightProportialToEnergy(maxWeight, lifeWidth, targetHkin, FermiLevel)
% Gauss distribution, peaks at Fermilevel with input variable maxWeight,
% and decays to unity at lifeWidth.
% Distribution Density: $$ y = \frac{1}{\sqrt{2\pi}\sigma} exp{[-\frac{(E_{tar} - E_f)^2}{k \cdot 2\sigma^2}]}$$
% where k is a renormalization parameter which satisfies $$ k = \frac{\pi L^2 N^2}{ln N}$$, others are standard gauss
% distribution probablility density function(pdf).
sigma = 1/(maxWeight * sqrt(2*pi()));
k = pi()*lifeWidth^2*maxWeight^2/log(maxWeight);
weight = zeros(1, length(targetHkin));
for numIdx = 1: length(targetHkin)
    weight(numIdx) = 1/(sqrt(2*pi())*sigma)*exp(-(targetHkin(numIdx)-FermiLevel)^2/(k*2*sigma^2));
    if weight(numIdx) < 1
        weight(numIdx) = 1;
    end
end
end

function dist = calcDist(x_A, a, x)
y = -(x_A - a)*(x - x_A)/(sqrt(3)*x_A) + sqrt(3)*x_A;
dist = abs(sqrt((x-x_A)^2 + (y-sqrt(3)*x_A)^2) - 2);
% dist_Orig = sqrt(x^2 + y^2);
end
