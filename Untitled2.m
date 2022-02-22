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
%% Sphere Sampling
clear variables
r = 10;
u = rand(1, 1000);
v = rand(1, 1000);
theta = 2*pi*u;
phi = acos(2.*v - 1);
% x = r.*cos(theta).*sin(phi);
% y = r.*sin(theta).*sin(phi);
% z = r.*cos(phi);
[x, y, z] = sph2cart(theta, phi - pi/2, r);
scatter3(x, y, z);
%% Construct primitive cell
clear variables
vacuumLength = 20;
numberOfLayers = 2;
% Lattice Constant Database:
% VS2:
a = 3.188283699968635;
% VSe2:
% a = 3.3171310374808005;
% VTe2:
% a = 3.579672793241685;
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
% Struct Database:
% VS2:
d = 1.467435675425866;   % h-phase:N/A, t-phase:1.467435675425866
D = 3.002888428993187;   % h-phase:N/A, t-phase:3.002888428993187
% VSe2:
% d = 1.58106107700766;   % h-phase:1.59498653898723, t-phase:1.58106107700766
% D = 3.42122839340594;   % h-phase:3.67176692979057, t-phase:3.12122839340594
% VTe2:
% d = 1.725763872572049;   % h-phase:N/A, t-phase:1.725763872572049
% D = 1.92122839340594;   % h-phase:N/A, t-phase:3.305618713378540

% if isfield(distance, 'mean_d')
%     d = distance.mean_d;
%     if isfield(distance, 'D_X2_X3')
%         D = distance.D_X2_X3;
%     end
% end
                        
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
fprintf(fileId, "   V    S\n");
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
%% Temp, Data Process for CPS Abstract
clear variables;
load('querySummary_for_abstract.mat');
zeroHorizontalLine = [linspace(lowerBound, upperBound, 100); zeros(1, 100)];
fig = figure();
hold on;
plot(VSe2_pri(:, 1), VSe2_pri(:, 4), '-o', 'Color', '#0072BD', 'LineWidth', 1.25, 'MarkerFaceColor', '#0072BD');
plot(VSe2_sup(:, 1), VSe2_sup(:, 4), '--o', 'Color', '#0072BD', 'LineWidth', 1.25, 'MarkerFaceColor', '#0072BD');
plot(VS2_pri(:, 1), VS2_pri(:, 4), '-o', 'Color', '#D95319', 'LineWidth', 1.25, 'MarkerFaceColor', '#D95319');
plot(VS2_sup(:, 1), VS2_sup(:, 4), '--o', 'Color', '#D95319', 'LineWidth', 1.25, 'MarkerFaceColor', '#D95319');
plot(VTe2_pri(:, 1), VTe2_pri(:, 4), '-o', 'Color', '#7E2F8E', 'LineWidth', 1.25, 'MarkerFaceColor', '#7E2F8E');
plot(VTe2_sup(:, 1), VTe2_sup(:, 4), '--o', 'Color', '#7E2F8E', 'LineWidth', 1.25, 'MarkerFaceColor', '#7E2F8E');
plot(zeroHorizontalLine(1, :), zeroHorizontalLine(2, :), ':', 'Color', 'black', 'LineWidth', 1.5);
% scatter(VS2_pri(:, 1), VS2_pri(:, 4), 48, 'filled', 'o', 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#D95319');
% scatter(VSe2_pri(:, 1), VSe2_pri(:, 4), 48, 'filled', 'o', 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#0072BD');
% scatter(VTe2_pri(:, 1), VTe2_pri(:, 4), 48, 'filled', 'o', 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#7E2F8E');
% scatter(VS2_sup(:, 1), VS2_sup(:, 4), 48, 'filled', 'o', 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#D95319');
% scatter(VSe2_sup(:, 1), VSe2_sup(:, 4), 48, 'filled', 'o', 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#0072BD');
% scatter(VTe2_sup(:, 1), VTe2_sup(:, 4), 48, 'filled', 'o', 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#7E2F8E');
% scatter(VS2_pri(1, 1), VS2_pri(1, 4), 48, 'x', 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#D95319');
% scatter(VSe2_pri(1, 1), VSe2_pri(1, 4), 48, 'x', 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#0072BD');
% scatter(VTe2_pri(1, 1), VTe2_pri(1, 4), 48, 'x', 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#7E2F8E');

hold off;
xlim([lowerBound, upperBound]);
legend({'VSe_2 PrimitiveCell', 'VSe_2 SuperCell', 'VS_2  PrimitiveCell', 'VS_2  SuperCell', 'VTe_2 PrimitiveCell', 'VTe_2 Supercell'});
ylabel("E_{FM} - E_{AFM}   (meV/atom)");
xlabel("D (Å)");
fig.Position = [0 0 800 450];
fig.Children(1).FontSize = 12;
fig.Children(1).Location = 'best';
fig.Children(2).Box = 'on';
fig.Children(2).FontSize = 14;
fig.Children(2).XAxis.MinorTick = 'on';
fig.Children(2).YAxis.MinorTick = 'on';
fig.Children(2).LineWidth = 1;
%% Data Process for CPS Poster
clear variables;
load('querySummary_for_abstract.mat');
% zeroHorizontalLine = [linspace(lowerBound, upperBound, 100); zeros(1, 100)];
fig = figure();
hold on;
plot(VSe2_pri(:, 1), VSe2_pri(:, 5), '-o', 'Color', '#0072BD', 'LineWidth', 1.25, 'MarkerFaceColor', '#0072BD');
plot(VSe2_pri(:, 1), VSe2_pri(:, 6), '--o', 'Color', '#0072BD', 'LineWidth', 1.25, 'MarkerFaceColor', '#0072BD');
plot(VSe2_sup(:, 1), VSe2_sup(:, 5), '-o', 'Color', '#FF8C00', 'LineWidth', 1.25, 'MarkerFaceColor', '#FF8C00');
plot(VSe2_sup(:, 1), VSe2_sup(:, 6), '--o', 'Color', '#FF8C00', 'LineWidth', 1.25, 'MarkerFaceColor', '#FF8C00');
plot(VS2_pri(:, 1), VS2_pri(:, 5), '-o', 'Color', '#A52A2A', 'LineWidth', 1.25, 'MarkerFaceColor', '#A52A2A');
plot(VS2_pri(:, 1), VS2_pri(:, 6), '--o', 'Color', '#A52A2A', 'LineWidth', 1.25, 'MarkerFaceColor', '#A52A2A');
plot(VS2_sup(:, 1), VS2_sup(:, 5), '-o', 'Color', '#7E2F8E', 'LineWidth', 1.25, 'MarkerFaceColor', '#7E2F8E');
plot(VS2_sup(:, 1), VS2_sup(:, 6), '--o', 'Color', '#7E2F8E', 'LineWidth', 1.25, 'MarkerFaceColor', '#7E2F8E');
plot(VTe2_pri(:, 1), VTe2_pri(:, 5), '-o', 'Color', '#77AC30', 'LineWidth', 1.25, 'MarkerFaceColor', '#77AC30');
plot(VTe2_pri(:, 1), VTe2_pri(:, 6), '--o', 'Color', '#77AC30', 'LineWidth', 1.25, 'MarkerFaceColor', '#77AC30');
plot(VTe2_sup(:, 1), VTe2_sup(:, 5), '-o', 'Color', '#696969', 'LineWidth', 1.25, 'MarkerFaceColor', '#696969');
plot(VTe2_sup(:, 1), VTe2_sup(:, 6), '--o', 'Color', '#696969', 'LineWidth', 1.25, 'MarkerFaceColor', '#696969');
% plot(zeroHorizontalLine(1, :), zeroHorizontalLine(2, :), ':', 'Color', 'black', 'LineWidth', 1.5);
hold off;
xlim([lowerBound, upperBound]);
ylim([0, 2.5])
legend({'VSe_2 PrimitiveCell FM', 'VSe_2 PrimitiveCell AFM', 'VSe_2 SuperCell FM', 'VSe_2 SuperCell AFM', 'VS_2  PrimitiveCell FM', 'VS_2  PrimitiveCell AFM', 'VS_2  SuperCell FM', 'VS_2  SuperCell AFM', 'VTe_2 PrimitiveCell FM', 'VTe_2 PrimitiveCell AFM', 'VTe_2 Supercell FM', 'VTe_2 Supercell AFM'});
ylabel("m   (\mu_{B}/atom)");
xlabel("D (Å)");

fig.Position = [0 0 800 450];
fig.Children(1).FontSize = 12;
fig.Children(1).Location = 'best';
fig.Children(2).Box = 'on';
fig.Children(2).FontSize = 14;
fig.Children(2).XAxis.MinorTick = 'on';
fig.Children(2).YAxis.MinorTick = 'on';
fig.Children(2).LineWidth = 1;

fig.Children(1).FontSize = 8;
fig.Children(1).NumColumns = 2;
fig.Children(1).Location = 'northeast';
%% Plot SuperCell/primitiveCell data, determine critical point
clear variables;
load('querySummary.mat');
combinedMat = [VS2_pri; VS2_sup; VSe2_pri; VSe2_sup; VTe2_pri; VTe2_sup];
lowerBound = min(combinedMat(:, 1));
upperBound = max(combinedMat(:, 1));
zeroHorizontalLine = [linspace(lowerBound, upperBound, 100); zeros(1, 100)];
fig = figure();
hold on;
plot(VSe2_pri(:, 1), VSe2_pri(:, 4), '-', 'Color', '#0072BD');
plot(VSe2_sup(:, 1), VSe2_sup(:, 4), '--', 'Color', '#0072BD');
plot(VS2_pri(:, 1), VS2_pri(:, 4), '-', 'Color', '#D95319');
plot(VS2_sup(:, 1), VS2_sup(:, 4), '--', 'Color', '#D95319');
plot(VTe2_pri(:, 1), VTe2_pri(:, 4), '-', 'Color', '#7E2F8E');
plot(VTe2_sup(:, 1), VTe2_sup(:, 4), '--', 'Color', '#7E2F8E');
plot(zeroHorizontalLine(1, :), zeroHorizontalLine(2, :), ':', 'Color', 'black');
scatter(VS2_pri(4, 1), VS2_pri(4, 4), 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#D95319');
scatter(VSe2_pri(1, 1), VSe2_pri(1, 4), 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#0072BD');
scatter(VTe2_pri(3, 1), VTe2_pri(3, 4), 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#7E2F8E');
hold off;
xlim([lowerBound, upperBound]);
legend({'VSe2 PrimitiveCell', 'VSe2 SuperCell', 'VS2 PrimitiveCell', 'VS2 SuperCell', 'VTe2 PrimitiveCell', 'VTe2 Supercell'});
title('Energy Difference(meV/atom)');

fig2 = figure();
hold on;
plot(VSe2_pri(:, 1), VSe2_pri(:, 5), '-', 'Color', '#0072BD');
plot(VSe2_pri(:, 1), VSe2_pri(:, 6), '--', 'Color', '#0072BD');
plot(VSe2_sup(:, 1), VSe2_sup(:, 5), '-', 'Color', '#D95319');
plot(VSe2_sup(:, 1), VSe2_sup(:, 6), '--', 'Color', '#D95319');
plot(VS2_pri(:, 1), VS2_pri(:, 5), '-', 'Color', '#EDB120');
plot(VS2_pri(:, 1), VS2_pri(:, 6), '--', 'Color', '#EDB120');
plot(VS2_sup(:, 1), VS2_sup(:, 5), '-', 'Color', '#7E2F8E');
plot(VS2_sup(:, 1), VS2_sup(:, 6), '--', 'Color', '#7E2F8E');
plot(VTe2_pri(:, 1), VTe2_pri(:, 5), '-', 'Color', '#77AC30');
plot(VTe2_pri(:, 1), VTe2_pri(:, 6), '--', 'Color', '#77AC30');
plot(VTe2_sup(:, 1), VTe2_sup(:, 5), '-', 'Color', '#4DBEEE');
plot(VTe2_sup(:, 1), VTe2_sup(:, 6), '--', 'Color', '#4DBEEE');
scatter(VS2_pri(4, 1), VS2_pri(4, 5), 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#EDB120');
scatter(VSe2_pri(1, 1), VSe2_pri(1, 5), 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#0072BD');
scatter(VTe2_pri(3, 1), VTe2_pri(3, 5), 'MarkerEdgeColor', 'None', 'MarkerFaceColor', '#77AC30');
hold off;
xlim([lowerBound, upperBound]);
legend({'VSe2 PrimitiveCell FM', 'VSe2 PrimitiveCell AFM', 'VSe2 SuperCell FM', 'VSe2 SuperCell AFM', ...
    'VS2 PrimitiveCell FM', 'VS2 PrimitiveCell AFM', 'VS2 SuperCell FM', 'VS2 SuperCell AFM',  ...
    'VTe2 PrimitiveCell FM', 'VTe2 PrimitiveCell AFM', 'VTe2 Supercell FM', 'VTe2 Supercell AFM'});
%% Interpolation
x = flip(VTe2_sup(:, 1));
v = flip(VTe2_sup(:, 4));
F = griddedInterpolant(x, v);
lowerBound = min(x);
upperBound = max(x);
xq = linspace(lowerBound, upperBound, 1000);
vq = F(xq);
[localMin, localMinIdx] = mink(abs(vq), 10);
disp(xq(localMinIdx));
%% GPU benchmark
% mat = rand(1000);
tic;
for i = 1:10
    dev_matA = gpuArray(rand(2000));
    dev_matB = gpuArray(rand(2000));
    dev_matC = mtimes(dev_matA, dev_matB);
    toc
end
%% Distance dependency of DOS, Visualization of DOS/projected DOS from VASP Plotter
clear variables
load('dosData.mat');
% Plot Total DOS
fig = figure();
hold on;
plot(data312.totalDOS.energy, data312.totalDOS.dos, 'LineWidth', 0.8);
plot(data232.totalDOS.energy, data232.totalDOS.dos, 'LineWidth', 0.8);
plot(data212.totalDOS.energy, data212.totalDOS.dos, 'LineWidth', 0.8);
plot(data192.totalDOS.energy, data192.totalDOS.dos, 'LineWidth', 0.8);
hold off;
xlim([-4 4]);
legend({'D = 3.12', 'D = 2.32', 'D = 2.12', 'D = 1.92'});
% Plot projected DOS, V_dxy
fig2 = figure();
hold on;
plot(data312.projectDOS.projectedDOSforEachAtom{1}.energy - data312.fermiLevel, ...
    data312.projectDOS.projectedDOSforEachAtom{1}.dxy_DOS, 'LineWidth', 0.8);
plot(data232.projectDOS.projectedDOSforEachAtom{1}.energy - data232.fermiLevel, ...
    data232.projectDOS.projectedDOSforEachAtom{1}.dxy_DOS, 'LineWidth', 0.8);
plot(data212.projectDOS.projectedDOSforEachAtom{1}.energy - data212.fermiLevel, ...
    data212.projectDOS.projectedDOSforEachAtom{1}.dxy_DOS, 'LineWidth', 0.8);
plot(data192.projectDOS.projectedDOSforEachAtom{1}.energy - data192.fermiLevel, ...
    data192.projectDOS.projectedDOSforEachAtom{1}.dxy_DOS, 'LineWidth', 0.8);
hold off;
xlim([-4 4]);
legend({'D = 3.12', 'D = 2.32', 'D = 2.12', 'D = 1.92'});
title('V dxy projectedDOS')
% Plot projected DOS, V_dyz
fig3 = figure();
hold on;
plot(data312.projectDOS.projectedDOSforEachAtom{1}.energy - data312.fermiLevel, ...
    data312.projectDOS.projectedDOSforEachAtom{1}.dyz_DOS, 'LineWidth', 0.8);
plot(data232.projectDOS.projectedDOSforEachAtom{1}.energy - data232.fermiLevel, ...
    data232.projectDOS.projectedDOSforEachAtom{1}.dyz_DOS, 'LineWidth', 0.8);
plot(data212.projectDOS.projectedDOSforEachAtom{1}.energy - data212.fermiLevel, ...
    data212.projectDOS.projectedDOSforEachAtom{1}.dyz_DOS, 'LineWidth', 0.8);
plot(data192.projectDOS.projectedDOSforEachAtom{1}.energy - data192.fermiLevel, ...
    data192.projectDOS.projectedDOSforEachAtom{1}.dyz_DOS, 'LineWidth', 0.8);
hold off;
xlim([-4 4]);
legend({'D = 3.12', 'D = 2.32', 'D = 2.12', 'D = 1.92'});
title('V dyz projectedDOS')
% Plot projected DOS, V_dz2_r2
fig4 = figure();
hold on;
plot(data312.projectDOS.projectedDOSforEachAtom{1}.energy - data312.fermiLevel, ...
    data312.projectDOS.projectedDOSforEachAtom{1}.dz2_r2__DOS, 'LineWidth', 0.8);
plot(data232.projectDOS.projectedDOSforEachAtom{1}.energy - data232.fermiLevel, ...
    data232.projectDOS.projectedDOSforEachAtom{1}.dz2_r2__DOS, 'LineWidth', 0.8);
plot(data212.projectDOS.projectedDOSforEachAtom{1}.energy - data212.fermiLevel, ...
    data212.projectDOS.projectedDOSforEachAtom{1}.dz2_r2__DOS, 'LineWidth', 0.8);
plot(data192.projectDOS.projectedDOSforEachAtom{1}.energy - data192.fermiLevel, ...
    data192.projectDOS.projectedDOSforEachAtom{1}.dz2_r2__DOS, 'LineWidth', 0.8);
hold off;
xlim([-4 4]);
legend({'D = 3.12', 'D = 2.32', 'D = 2.12', 'D = 1.92'});
title('V dz2-r2 projectedDOS')
% Plot projected DOS, V_dxz
fig5 = figure();
hold on;
plot(data312.projectDOS.projectedDOSforEachAtom{1}.energy - data312.fermiLevel, ...
    data312.projectDOS.projectedDOSforEachAtom{1}.dxz_DOS, 'LineWidth', 0.8);
plot(data232.projectDOS.projectedDOSforEachAtom{1}.energy - data232.fermiLevel, ...
    data232.projectDOS.projectedDOSforEachAtom{1}.dxz_DOS, 'LineWidth', 0.8);
plot(data212.projectDOS.projectedDOSforEachAtom{1}.energy - data212.fermiLevel, ...
    data212.projectDOS.projectedDOSforEachAtom{1}.dxz_DOS, 'LineWidth', 0.8);
plot(data192.projectDOS.projectedDOSforEachAtom{1}.energy - data192.fermiLevel, ...
    data192.projectDOS.projectedDOSforEachAtom{1}.dxz_DOS, 'LineWidth', 0.8);
hold off;
xlim([-4 4]);
legend({'D = 3.12', 'D = 2.32', 'D = 2.12', 'D = 1.92'});
title('V dxz projectedDOS')
% Plot projected DOS, V_dx2-y2
fig6 = figure();
hold on;
plot(data312.projectDOS.projectedDOSforEachAtom{1}.energy - data312.fermiLevel, ...
    data312.projectDOS.projectedDOSforEachAtom{1}.dx2_y2_DOS, 'LineWidth', 0.8);
plot(data232.projectDOS.projectedDOSforEachAtom{1}.energy - data232.fermiLevel, ...
    data232.projectDOS.projectedDOSforEachAtom{1}.dx2_y2_DOS, 'LineWidth', 0.8);
plot(data212.projectDOS.projectedDOSforEachAtom{1}.energy - data212.fermiLevel, ...
    data212.projectDOS.projectedDOSforEachAtom{1}.dx2_y2_DOS, 'LineWidth', 0.8);
plot(data192.projectDOS.projectedDOSforEachAtom{1}.energy - data192.fermiLevel, ...
    data192.projectDOS.projectedDOSforEachAtom{1}.dx2_y2_DOS, 'LineWidth', 0.8);
hold off;
xlim([-4 4]);
legend({'D = 3.12', 'D = 2.32', 'D = 2.12', 'D = 1.92'});
title('V dxz projectedDOS')
%% Aero-Orbit
clear variables;
v1 = 42100;
G = 6.67259e-11;
M = 1.989e30;
r1 = 1.496e11;
r_target = 1.524*r1;
r2 = (2*G*M - sqrt(4*G^2*M^2 + 4*(v1^2 + 2*G*M/r1)*v1^2*r1^2))/(2*(v1^2 + 2*G*M/r1));
v2 = v1*r1/r2;
r_inAU = r2/r1;
epsilon = (r2 - r1)/(r2 + r1);
p = (2*r2*r1)/(r2+r1);
theta = linspace(0, 2*pi, 200);
r = p./(1 + epsilon.*cos(theta));
x = r.*cos(theta);
y = r.*sin(theta);
% scatter(x, y);

%% Read Data from MySQL database
clear variables;
% Parameter Setting: material, cellType, layerDistance
conn = database('nm_dos','sces','');
material = '''VSe2''';
cellType = '''Primary''';
% cellType = '''Super''';
if strcmp(cellType,'''Primary''')
    totAtomNumber = 6;
else
    totAtomNumber = 42;
end
layerDistanceDatabase = [2.12 3.12];
totLayerDistanceNumber = length(layerDistanceDatabase);
% Judge structure type, allocate space
fmdata = cell(2, totAtomNumber);
afmdata = cell(2, totAtomNumber);
nmdata = cell(2, totAtomNumber);
for layerIdx = 1: totLayerDistanceNumber
    layerDistance = layerDistanceDatabase(layerIdx);
    for atomNumber = 1: totAtomNumber
        magConfiguration = '''FM''';
        query = ['select * from DOSCAR ' ...
            'where  material = ', material, ' and cellType = ', cellType, ' and '...
            'magConfiguration = ',magConfiguration, ' and atomIdx = ', num2str(atomNumber), ...
            ' and round(layerDistance, 2) = ', num2str(layerDistance) ,' ;'];
        fmdata{layerIdx, atomNumber} = fetch(conn, query);
        
        magConfiguration = '''AFM''';
        query = ['select * from DOSCAR ' ...
            'where  material = ', material, ' and cellType = ', cellType, ' and '...
            'magConfiguration = ',magConfiguration, ' and atomIdx = ', num2str(atomNumber), ...
            ' and round(layerDistance, 2) = ', num2str(layerDistance) ,' ;'];
        afmdata{layerIdx, atomNumber} = fetch(conn, query);
        
        magConfiguration = '''NM''';
        query = ['select * from nm_DOSCAR ' ...
            'where  material = ', material, ' and cellType = ', cellType, ' and '...
            'magConfiguration = ',magConfiguration, ' and atomIdx = ', num2str(atomNumber), ...
            ' and round(layerDistance, 2) = ', num2str(layerDistance) ,' ;'];
        nmdata{layerIdx, atomNumber} = fetch(conn, query);
    end
end
close(conn);

%% Process total dos
% Get number of V atoms
switch totAtomNumber
    case 6
        totVatomNumber = 2;
        cellType = 'Primary';
    case 42
        totVatomNumber = 14;
        cellType = 'Super';
end
% Calculate average value of V atom from NM, FM, AFM Configuration
totalDosFmAverage = cell(totLayerDistanceNumber, 1);
totalDosAFmAverage = cell(totLayerDistanceNumber, 1);
totalDosNmAverage = cell(totLayerDistanceNumber, 1);

for layerIdx = 1: length(layerDistanceDatabase)
    % Calulate the FM total DOS average data
    totalDosFmAverage{layerIdx, 1}.energy = fmdata{layerIdx, 1}.energy;
    totalDosFmAverage{layerIdx, 1}.tot_up = fmdata{layerIdx, 1}.tot_up;
    totalDosFmAverage{layerIdx, 1}.tot_down = fmdata{layerIdx, 1}.tot_down;
    % Calulate the AFM total DOS average data
    totalDosAFmAverage{layerIdx, 1}.energy = afmdata{layerIdx, 1}.energy;
    totalDosAFmAverage{layerIdx, 1}.tot_up = afmdata{layerIdx, 1}.tot_up;
    totalDosAFmAverage{layerIdx, 1}.tot_down = afmdata{layerIdx, 1}.tot_down;
    % Calculate the NM total DOS average data
    totalDosNmAverage{layerIdx, 1}.energy = nmdata{layerIdx, 1}.energy;
    totalDosNmAverage{layerIdx, 1}.tot = nmdata{layerIdx, 1}.tot;   
    % Data visualization for V atoms
    fig = figure();
    fig.WindowState = 'maximized';
    hold on;
    fmLineUp = plot(totalDosFmAverage{layerIdx, 1}.energy, totalDosFmAverage{layerIdx, 1}.tot_up, 'Color', 'blue');
    plot(totalDosFmAverage{layerIdx, 1}.energy, -totalDosFmAverage{layerIdx, 1}.tot_down, 'Color', 'blue');
    afmLineUp = plot(totalDosAFmAverage{layerIdx, 1}.energy, totalDosAFmAverage{layerIdx, 1}.tot_up, 'Color', 'red');
    plot(totalDosAFmAverage{layerIdx, 1}.energy, -totalDosAFmAverage{layerIdx, 1}.tot_down, 'Color', 'red');
    nmLine = plot(totalDosNmAverage{layerIdx, 1}.energy, totalDosNmAverage{layerIdx, 1}.tot*2/3, 'Color', 'black');
    xlim([-6 6]);
    plot([0 0], [-20 20], '--', 'Color', [0.25 0.25 0.25]);
    hold off;
    legend([fmLineUp afmLineUp nmLine], {'FM', 'AFM', 'NM'});
    title(['total DOS at D = ', num2str(layerDistanceDatabase(layerIdx)), ', Cell Type: ', cellType]);
    imgName = ['./images/', 'tdos', cellType, 'Dis', num2str(layerDistanceDatabase(layerIdx)), '.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');
end

%% Process projected DOS for specific atoms

orbitalDataBase = ["s", "p_y", "p_z", "p_x", "d_xy", "d_yz", "d_x2_r2", "d_xz", "d_x2_y2"];
% Parameter Setting: Orbital Selection
orbitalPickAtomV = [orbitalDataBase(5), orbitalDataBase(7)];
orbitalPickAtomSe = [orbitalDataBase(2), orbitalDataBase(3), orbitalDataBase(4)];
orbitalPickAtomVLength = length(orbitalPickAtomV);
orbitalPickAtomSeLength = length(orbitalPickAtomSe);
% Initialize DOS for V atom
projectedDosUpFmMatAtomV = cell(layerIdx, orbitalPickAtomVLength);
projectedDosDownFmMatAtomV = cell(layerIdx, orbitalPickAtomVLength);
projectedDosFmAverageAtomV = cell(layerIdx, orbitalPickAtomVLength);
projectedDosUpAFmMatAtomV = cell(layerIdx, orbitalPickAtomVLength);
projectedDosDownAFmMatAtomV = cell(layerIdx, orbitalPickAtomVLength);
projectedDosAFmAverageAtomV = cell(layerIdx, orbitalPickAtomVLength);
projectedDosNmMatAtomV = cell(layerIdx, orbitalPickAtomVLength);
projectedDosNmAverageAtomV = cell(layerIdx, orbitalPickAtomVLength);
% Initialize DOS for Se atom
projectedDosUpFmMatAtomSe = cell(layerIdx, orbitalPickAtomSeLength);
projectedDosDownFmMatAtomSe = cell(layerIdx, orbitalPickAtomSeLength);
projectedDosFmAverageAtomSe = cell(layerIdx, orbitalPickAtomSeLength);
projectedDosUpAFmMatAtomSe = cell(layerIdx, orbitalPickAtomSeLength);
projectedDosDownAFmMatAtomSe = cell(layerIdx, orbitalPickAtomSeLength);
projectedDosAFmAverageAtomSe = cell(layerIdx, orbitalPickAtomSeLength);
projectedDosNmMatAtomSe = cell(layerIdx, orbitalPickAtomSeLength);
projectedDosNmAverageAtomSe = cell(layerIdx, orbitalPickAtomSeLength);

for layerIdx = 1: length(layerDistanceDatabase)
    % Calculate DOS for V atom
    for orbitPickAtomVIdx = 1: length(orbitalPickAtomV)
        % Calculate the FM total DOS average data
        for atomIdx = 1: totVatomNumber
            projectedDosUpFmMatAtomV{layerIdx, orbitPickAtomVIdx} = ...
                cat(2, projectedDosUpFmMatAtomV{layerIdx, orbitPickAtomVIdx}, ...
                    fmdata{layerIdx, atomIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_up'));
            projectedDosDownFmMatAtomV{layerIdx, orbitPickAtomVIdx} = ...
                cat(2, projectedDosDownFmMatAtomV{layerIdx, orbitPickAtomVIdx}, ...
                    fmdata{layerIdx, atomIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_down'));
        end
        projectedDosFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.energy = fmdata{layerIdx, 1}.energy;
        projectedDosFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_up') = ...
            mean(projectedDosUpFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, 1: totVatomNumber), 2);
        projectedDosFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.energy = fmdata{layerIdx, 1}.energy;
        projectedDosFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_down') = ...
            mean(projectedDosDownFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, 1: totVatomNumber), 2);
        % Calulate the AFM total DOS average data
        for atomIdx = 1: totVatomNumber
            projectedDosUpAFmMatAtomV{layerIdx, orbitPickAtomVIdx} = ...
                cat(2, projectedDosUpAFmMatAtomV{layerIdx, orbitPickAtomVIdx}, ...
                    afmdata{layerIdx, atomIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_up'));
            projectedDosDownAFmMatAtomV{layerIdx, orbitPickAtomVIdx} = ...
                cat(2, projectedDosDownAFmMatAtomV{layerIdx, orbitPickAtomVIdx}, ...
                    afmdata{layerIdx, atomIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_down'));
        end
        projectedDosAFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.energy = afmdata{layerIdx, 1}.energy;
        projectedDosAFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_up') = ...
            mean(projectedDosUpAFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, 1: totVatomNumber), 2);
        projectedDosAFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.energy = afmdata{layerIdx, 1}.energy;
        projectedDosAFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_down') = ...
            mean(projectedDosDownAFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, 1: totVatomNumber), 2);
        % Calulate the NM total DOS average data
        for atomIdx = 1: totVatomNumber
            projectedDosNmMatAtomV{layerIdx, orbitPickAtomVIdx} = ...
                cat(2, projectedDosNmMatAtomV{layerIdx, orbitPickAtomVIdx}, ...
                    nmdata{layerIdx, atomIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)));
        end
        projectedDosNmAverageAtomV{layerIdx, orbitPickAtomVIdx}.energy = nmdata{layerIdx, 1}.energy;
        projectedDosNmAverageAtomV{layerIdx, orbitPickAtomVIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)) = ...
            mean(projectedDosNmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, 1: totVatomNumber), 2);
    end
    % Calculate DOS for Se atom
    for orbitPickAtomSeIdx = 1: length(orbitalPickAtomSe)
        % Calculate the FM total DOS average data
        for atomIdx = totVatomNumber + 1: 3*totVatomNumber
            projectedDosUpFmMatAtomSe{layerIdx, orbitPickAtomSeIdx} = ...
                cat(2, projectedDosUpFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}, ...
                    fmdata{layerIdx, atomIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_up'));
            projectedDosDownFmMatAtomSe{layerIdx, orbitPickAtomSeIdx} = ...
                cat(2, projectedDosDownFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}, ...
                    fmdata{layerIdx, atomIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_down'));
        end
        projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy = fmdata{layerIdx, 1}.energy;
        projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_up') = ...
            mean(projectedDosUpFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, 1: 2*totVatomNumber), 2);
        projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy = fmdata{layerIdx, 1}.energy;
        projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_down') = ...
            mean(projectedDosDownFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, 1: 2*totVatomNumber), 2);
        % Calulate the AFM total DOS average data
        for atomIdx = totVatomNumber + 1: 3*totVatomNumber
            projectedDosUpAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx} = ...
                cat(2, projectedDosUpAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}, ...
                    afmdata{layerIdx, atomIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_up'));
            projectedDosDownAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx} = ...
                cat(2, projectedDosDownAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}, ...
                    afmdata{layerIdx, atomIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_down'));
        end
        projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy = afmdata{layerIdx, 1}.energy;
        projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_up') = ...
            mean(projectedDosUpAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, 1: 2*totVatomNumber), 2);
        projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy = afmdata{layerIdx, 1}.energy;
        projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_down') = ...
            mean(projectedDosDownAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, 1: 2*totVatomNumber), 2);
        % Calulate the NM total DOS average data
        for atomIdx = totVatomNumber + 1: 3*totVatomNumber
            projectedDosNmMatAtomSe{layerIdx, orbitPickAtomSeIdx} = ...
                cat(2, projectedDosNmMatAtomSe{layerIdx, orbitPickAtomSeIdx}, ...
                    nmdata{layerIdx, atomIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)));
        end
        projectedDosNmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy = nmdata{layerIdx, 1}.energy;
        projectedDosNmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)) = ...
            mean(projectedDosNmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, 1: 2*totVatomNumber), 2);
    end
end

%% Data Visualization for processed DOS
% Parameter Setting: Label of atoms that to be visualized, atom index obeys
% the same rule in VESTA.
visualPickAtomV = [1 2];
visualPickAtomSe = [3 4 5 6];

% Visualization for each layer distance (FM v.s. AFM v.s NM)
for layerIdx = 1: length(layerDistanceDatabase)
    % Visualize projected DOS with selected orbits for atom V
    fig = figure();
    fig.WindowState = 'maximized';
    idx = 1;
    for orbitPickAtomVIdx = 1: length(orbitalPickAtomV)
        for pickAtomVIdx = 1: length(visualPickAtomV)
            subplot(length(orbitalPickAtomV), length(visualPickAtomV), idx);
            hold on;
            fmLineUp = plot(fmdata{layerIdx, 1}.energy, ...
                projectedDosUpFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, visualPickAtomV(pickAtomVIdx)), 'Color', 'blue');
            plot(fmdata{layerIdx, 1}.energy, ...
                -projectedDosDownFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, visualPickAtomV(pickAtomVIdx)), 'Color', 'blue');
            afmLineUp = plot(afmdata{layerIdx, 1}.energy, ...
                projectedDosUpAFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, visualPickAtomV(pickAtomVIdx)), 'Color', 'red');
            plot(afmdata{layerIdx, 1}.energy, ...
                -projectedDosDownAFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, visualPickAtomV(pickAtomVIdx)), 'Color', 'red');
            nmLine = plot(nmdata{layerIdx, 1}.energy, ...
                projectedDosNmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, visualPickAtomV(pickAtomVIdx))*2/3, 'Color', 'black');
            xlim([-8 8]);
            plot([0 0], [-8 8], '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            legend([fmLineUp afmLineUp nmLine], {'FM', 'AFM', 'NM'});
            title({[num2str(orbitalPickAtomV(orbitPickAtomVIdx)), ' DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: V, Seq.Num ',num2str(visualPickAtomV(pickAtomVIdx))]});
            idx = idx + 1;
        end
    end
    imgName = ['./images/', 'pdos', cellType, num2str(orbitalPickAtomV(orbitPickAtomVIdx))...
        ,'Dis', num2str(layerDistanceDatabase(layerIdx)), 'Atom', 'V', 'SeqNum',num2str(visualPickAtomV(pickAtomVIdx)) ,'.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');
    % Visualize projected DOS with selected orbits for atom Se
    fig = figure();
    fig.WindowState = 'maximized';
    idx = 1;
    for orbitPickAtomSeIdx = 1: length(orbitalPickAtomSe)
        for pickAtomSeIdx = 1: length(visualPickAtomSe)
            subplot(length(orbitalPickAtomSe), length(visualPickAtomSe), idx);
            hold on;
            fmLineUp = plot(fmdata{layerIdx, 1}.energy, ...
                projectedDosUpFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx) - totVatomNumber), 'Color', 'blue');
            plot(fmdata{layerIdx, 1}.energy, ...
                -projectedDosDownFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx) - totVatomNumber), 'Color', 'blue');
            afmLineUp = plot(afmdata{layerIdx, 1}.energy, ...
                projectedDosUpAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx) - totVatomNumber), 'Color', 'red');
            plot(afmdata{layerIdx, 1}.energy, ...
                -projectedDosDownAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx) - totVatomNumber), 'Color', 'red');
            nmLine = plot(nmdata{layerIdx, 1}.energy, ...
                projectedDosNmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx) - totVatomNumber)*2/3, 'Color', 'black');
            xlim([-8 8]);
            plot([0 0], [-5 5], '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            legend([fmLineUp afmLineUp nmLine], {'FM', 'AFM', 'NM'});
            title({[num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ' DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: Se, Seq.Num ',num2str(visualPickAtomSe(pickAtomSeIdx))]});
            idx = idx + 1;
        end
    end
    imgName = ['./images/', 'pdos', cellType, num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ...
        'Dis', num2str(layerDistanceDatabase(layerIdx)), 'Atom', 'Se', 'SeqNum',num2str(visualPickAtomSe(pickAtomSeIdx)), '.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');
end

%% Data visualization for each layer distance (projected DOS for each atom v.s. average DOS)
for layerIdx = 1: length(layerDistanceDatabase)
    % Visualize projected DOS with selected orbits for atom V
    fig = figure();
    fig.WindowState = 'maximized';
    idx = 1;
    % Visualize FM data
    for orbitPickAtomVIdx = 1: length(orbitalPickAtomV)
        for pickAtomVIdx = 1: length(visualPickAtomV)
            subplot(length(orbitalPickAtomV), length(visualPickAtomV), idx);
            hold on;
            fmLineUpAtom = plot(fmdata{layerIdx, 1}.energy, ...
                projectedDosUpFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, visualPickAtomV(pickAtomVIdx)), 'Color', 'blue');
            plot(fmdata{layerIdx, 1}.energy, ...
                -projectedDosDownFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, visualPickAtomV(pickAtomVIdx)), 'Color', 'blue');
            fmLineUpAverage = plot(projectedDosFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.energy, ...
                projectedDosFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_up'), 'Color', 'red');
            plot(projectedDosFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.energy, ...
                -projectedDosFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_down'), 'Color', 'red');
            xlim([-8 8]);
            plot([0 0], [-8 8], '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            legend([fmLineUpAtom fmLineUpAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomV(orbitPickAtomVIdx)), '  FM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: V, Seq.Num ',num2str(visualPickAtomV(pickAtomVIdx))]});
            idx = idx + 1;
        end
    end
    imgName = ['./images/', 'magFm', cellType, num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ...
        'Dis', num2str(layerDistanceDatabase(layerIdx)), 'Atom', 'V', 'SeqNum',num2str(visualPickAtomV(pickAtomVIdx)), '.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');

    % Visualize AFM data
    fig = figure();
    fig.WindowState = 'maximized';
    idx = 1;
    for orbitPickAtomVIdx = 1: length(orbitalPickAtomV)
        for pickAtomVIdx = 1: length(visualPickAtomV)
            subplot(length(orbitalPickAtomV), length(visualPickAtomV), idx);
            hold on;
            afmLineUpAtom = plot(afmdata{layerIdx, 1}.energy, ...
                projectedDosUpAFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, visualPickAtomV(pickAtomVIdx)), 'Color', 'blue');
            plot(afmdata{layerIdx, 1}.energy, ...
                -projectedDosDownAFmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, visualPickAtomV(pickAtomVIdx)), 'Color', 'blue');
            afmLineUpAverage = plot(projectedDosAFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.energy, ...
                projectedDosAFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_up'), 'Color', 'red');
            plot(projectedDosAFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.energy, ...
                -projectedDosAFmAverageAtomV{layerIdx, orbitPickAtomVIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)+'_down'), 'Color', 'red');
            xlim([-8 8]);
            plot([0 0], [-8 8], '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            legend([afmLineUpAtom afmLineUpAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomV(orbitPickAtomVIdx)), ' AFM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: V, Seq.Num ',num2str(visualPickAtomV(pickAtomVIdx))]});
            idx = idx + 1;
        end
    end
    imgName = ['./images/', 'magAFm', cellType, num2str(orbitalPickAtomV(orbitPickAtomVIdx)), ...
        'Dis', num2str(layerDistanceDatabase(layerIdx)), 'Atom', 'V', 'SeqNum',num2str(visualPickAtomV(pickAtomVIdx)), '.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');

    % Visualize NM data
    fig = figure();
    fig.WindowState = 'maximized';
    idx = 1;
    for orbitPickAtomVIdx = 1: length(orbitalPickAtomV)
        for pickAtomVIdx = 1: length(visualPickAtomV)
            subplot(length(orbitalPickAtomV), length(visualPickAtomV), idx);
            hold on;
            nmLineAtom = plot(nmdata{layerIdx, 1}.energy, ...
                projectedDosNmMatAtomV{layerIdx, orbitPickAtomVIdx}(:, visualPickAtomV(pickAtomVIdx)), 'Color', 'blue');
            nmLineAverage = plot(projectedDosNmAverageAtomV{layerIdx, orbitPickAtomVIdx}.energy, ...
                projectedDosNmAverageAtomV{layerIdx, orbitPickAtomVIdx}.(orbitalPickAtomV(orbitPickAtomVIdx)), 'Color', 'red');
            xlim([-8 8]);
            plot([0 0], [-8 8], '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            legend([nmLineAtom nmLineAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomV(orbitPickAtomVIdx)), ' NM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: V, Seq.Num ',num2str(visualPickAtomV(pickAtomVIdx))]});
            idx = idx + 1;
        end
    end
    imgName = ['./images/', 'magNm', cellType, num2str(orbitalPickAtomV(orbitPickAtomVIdx)), ...
        'Dis', num2str(layerDistanceDatabase(layerIdx)), 'Atom', 'V', 'SeqNum',num2str(visualPickAtomV(pickAtomVIdx)), '.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');

    % Visualize projected DOS with selected orbits for atom Se
    fig = figure();
    fig.WindowState = 'maximized';
    idx = 1;
    % Visualize FM data
    for orbitPickAtomSeIdx = 1: length(orbitalPickAtomSe)
        for pickAtomSeIdx = 1: length(visualPickAtomSe)
            subplot(length(orbitalPickAtomSe), length(visualPickAtomSe), idx);
            hold on;
            fmLineUpAtom = plot(fmdata{layerIdx, 1}.energy, ...
                projectedDosUpFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx)), 'Color', 'blue');
            plot(fmdata{layerIdx, 1}.energy, ...
                -projectedDosDownFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx)), 'Color', 'blue');
            fmLineUpAverage = plot(projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy, ...
                projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_up'), 'Color', 'red');
            plot(projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy, ...
                -projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_down'), 'Color', 'red');
            xlim([-8 8]);
            plot([0 0], [-8 8], '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            legend([fmLineUpAtom fmLineUpAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), '  FM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: Se, Seq.Num ',num2str(visualPickAtomSe(pickAtomSeIdx))]});
            idx = idx + 1;
        end
    end
    imgName = ['./images/', 'magFm', cellType, num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ...
        'Dis', num2str(layerDistanceDatabase(layerIdx)), 'Atom', 'Se', 'SeqNum',num2str(visualPickAtomSe(pickAtomSeIdx)), '.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');

    % Visualize AFM data
    fig = figure();
    fig.WindowState = 'maximized';
    idx = 1;
    for orbitPickAtomSeIdx = 1: length(orbitalPickAtomSe)
        for pickAtomSeIdx = 1: length(visualPickAtomSe)
            subplot(length(orbitalPickAtomSe), length(visualPickAtomSe), idx);
            hold on;
            afmLineUpAtom = plot(afmdata{layerIdx, 1}.energy, ...
                projectedDosUpAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx)), 'Color', 'blue');
            plot(afmdata{layerIdx, 1}.energy, ...
                -projectedDosDownAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx)), 'Color', 'blue');
            afmLineUpAverage = plot(projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy, ...
                projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_up'), 'Color', 'red');
            plot(projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy, ...
                -projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_down'), 'Color', 'red');
            xlim([-8 8]);
            plot([0 0], [-8 8], '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            legend([afmLineUpAtom afmLineUpAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ' AFM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: Se, Seq.Num ',num2str(visualPickAtomSe(pickAtomSeIdx))]});
            idx = idx + 1;
        end
    end
    imgName = ['./images/', 'magAFm', cellType, num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ...
        'Dis', num2str(layerDistanceDatabase(layerIdx)), 'Atom', 'Se', 'SeqNum',num2str(visualPickAtomSe(pickAtomSeIdx)), '.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');

    % Visualize NM data
    fig = figure();
    fig.WindowState = 'maximized';
    idx = 1;
    for orbitPickAtomSeIdx = 1: length(orbitalPickAtomSe)
        for pickAtomSeIdx = 1: length(visualPickAtomSe)
            subplot(length(orbitalPickAtomSe), length(visualPickAtomSe), idx);
            hold on;
            nmLineAtom = plot(nmdata{layerIdx, 1}.energy, ...
                projectedDosNmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx)), 'Color', 'blue');
            nmLineAverage = plot(projectedDosNmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy, ...
                projectedDosNmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)), 'Color', 'red');
            xlim([-8 8]);
            plot([0 0], [-8 8], '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            legend([nmLineAtom nmLineAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ' NM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: Se, Seq.Num ',num2str(visualPickAtomSe(pickAtomSeIdx))]});
            idx = idx + 1;
        end
    end
    imgName = ['./images/', 'magNm', cellType, num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ...
        'Dis', num2str(layerDistanceDatabase(layerIdx)), 'Atom', 'Se', 'SeqNum',num2str(visualPickAtomSe(pickAtomSeIdx)), '.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');
end

%% GPU Benchmark
clear variables;
N=5000;
X = single(diag(rand(N,1)));
U = single(orth(rand(N,N)));
B = (U' * X * U).*300;
tic;
A = gpuArray(B);
[V,D] = eig(B);
toc;
%% Fermi Surface Plotting
clear variables;
% Control Flag Setting
% app.dirname = 'D:\OneDrive - tongji.edu.cn\Research\VSe2\Tphase\monolayer\nmAllBZ';
% app.dirname = 'D:\OneDrive - tongji.edu.cn\Research\VSe2\Tphase\bulk\nmOptVol';
app.dirname = 'D:\OneDrive - tongji.edu.cn\Research\VSe2\Tphase\bilayer\Rotation\21.7867892982618\change_D_Ueff1\nm\2.12';
app.tmp.spinPolarizedFlag = 0;
app.EIGENVAL_3D.bandIndexSelection = 1:392;
app.tmp.meshDensity = 200;
app.EIGENVAL_3D.kzPlaneSelection = 0;
app.tmp.fermiLevelOffset = 0.0;

% Import file
filename = strcat(app.dirname,'/EIGENVAL_3D');
if ~exist(filename,'file')
    [name, path] = uigetfile('*.*');
    filename = strcat(path,name);
    if ~exist(filename,'file')
        errordlg('Static calculation full-BZ zone EIGENVAL file not found','File Error');
        clear variables;
        return
    end
end

% Read POSCAR Data
filename = strcat(app.dirname,'/POSCAR');
if ~exist(filename,'file')
    errordlg('POSCAR file not found','File Error');
    clear variables;
    return
end

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
    app.POSCAR.stringMat(i, 1: length(tmp.numIdx)) = POSCARtemp(i, tmp.numIdx);
end
app.POSCAR.stringMat(ismissing(app.POSCAR.stringMat)) = "";
app.POSCAR.element.name = app.POSCAR.stringMat(6, :);
app.POSCAR.element.name( app.POSCAR.element.name=='' )= [];
app.POSCAR.element.num = str2double (app.POSCAR.stringMat(7, :) );
app.POSCAR.element.num(isnan(app.POSCAR.element.num)) = [];
app.POSCAR.element.totalAtomNumber = sum(app.POSCAR.element.num);
tmp.totAtomNumber = sum(app.POSCAR.element.num);
clear opts i POSCARtemp

% Read DOSCAR
filename = strcat(app.dirname,'/DOSCAR');
fileID = fopen(filename,'r');
if fileID == -1
    errordlg('NOSCAR file not found','File Error');
    clear variables;
    return
end

switch app.tmp.spinPolarizedFlag
    case 0
        opts = delimitedTextImportOptions("NumVariables", 4);
        try
            opts.DataLines = [6, 6+tmp.nedos];
        catch
            opts.DataLines = [6, 307];
        end
        opts.Delimiter = " ";
        opts.VariableNames = ["energy", "DOS", "intDOS", "fermiLevel"];
        opts.VariableTypes = ["double", "double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
        app.DOSCAR.dos = readtable(filename, opts);
        clear opts
    case 1
        opts = delimitedTextImportOptions("NumVariables", 5);
        try
            opts.DataLines = [6, 6+tmp.nedos];
        catch
            opts.DataLines = [6, 307];
        end
        opts.Delimiter = " ";
        opts.VariableNames = ["energy", "DOS_Up", "DOS_Down", "intDOS_Up", "intDOS_Down"];
        opts.VariableTypes = ["double", "double", "double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
        app.DOSCAR.dos = readtable(filename, opts);
        clear opts
end
app.DOSCAR.fermiLevel = app.DOSCAR.dos(1, 4);
app.DOSCAR.fermiLevel = table2array(app.DOSCAR.fermiLevel);
app.DOSCAR.dos(1, :) = [];

% Read EIGENVAL Data
filename = strcat(app.dirname,'/EIGENVAL_3D');
opts = delimitedTextImportOptions("NumVariables", 25);
opts.DataLines = [8, Inf];
opts.Delimiter = " ";
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5"];
opts.SelectedVariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];
opts = setvaropts(opts, [1, 2], "TrimNonNumeric", true);
opts = setvaropts(opts, [1, 2], "ThousandsSeparator", ",");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
app.EIGENVAL_3D.originalData = readtable(filename, opts);
clear opts

% Data Process
% Delete the whitespace line, get the number of bands & kpoints
% independently.
app.EIGENVAL_3D.processedData = table2cell(app.EIGENVAL_3D.originalData);
numIdx = cellfun(@(x) ~isnan(str2double(x)), app.EIGENVAL_3D.processedData);
app.EIGENVAL_3D.processedData(numIdx) = cellfun(@(x) {str2double(x)}, app.EIGENVAL_3D.processedData(numIdx));
numIdx = cellfun(@(x) isnan((x)), app.EIGENVAL_3D.processedData(:, 1));
app.EIGENVAL_3D.processedData(numIdx, :) = [];
numIdx = cellfun(@(x) mod((x),1), app.EIGENVAL_3D.processedData(:, 1));
i = 2;
while ~numIdx(i)
    i = i + 1;
end
app.EIGENVAL_3D.totBands = i - 2;
app.EIGENVAL_3D.totKpoints = length(app.EIGENVAL_3D.processedData)/(app.EIGENVAL_3D.totBands + 1);
% Initialize the data structure, read data formally
readBaseLine = 1;
app.EIGENVAL_3D.kpointsWeight = zeros (app.EIGENVAL_3D.totKpoints, 1);
app.EIGENVAL_3D.kpointsMesh = zeros (3, app.EIGENVAL_3D.totKpoints);
app.EIGENVAL_3D.bandEigenvalue = cell(app.EIGENVAL_3D.totBands, app.EIGENVAL_3D.totKpoints);
app.EIGENVAL_3D.bandOccupation = cell(app.EIGENVAL_3D.totBands, app.EIGENVAL_3D.totKpoints);

for orderKpoint = 1: app.EIGENVAL_3D.totKpoints
    app.EIGENVAL_3D.kpointsWeight(orderKpoint) = app.EIGENVAL_3D.processedData{readBaseLine, 4};
    app.EIGENVAL_3D.kpointsMesh(1, orderKpoint) = app.EIGENVAL_3D.processedData{readBaseLine, 1};
    app.EIGENVAL_3D.kpointsMesh(2, orderKpoint) = app.EIGENVAL_3D.processedData{readBaseLine, 2};
    app.EIGENVAL_3D.kpointsMesh(3, orderKpoint) = app.EIGENVAL_3D.processedData{readBaseLine, 3};
    readBaseLine = readBaseLine + 1;
    for orderBands = 1: app.EIGENVAL_3D.totBands
        app.EIGENVAL_3D.bandEigenvalue{orderBands, orderKpoint} = app.EIGENVAL_3D.processedData{readBaseLine + orderBands - 1, 2};
        app.EIGENVAL_3D.bandOccupation{orderBands, orderKpoint} = app.EIGENVAL_3D.processedData{readBaseLine + orderBands - 1, 3};
    end
    readBaseLine = readBaseLine + app.EIGENVAL_3D.totBands;
end
app.EIGENVAL_3D.kpointsMesh(4, :) = 1: length(app.EIGENVAL_3D.kpointsMesh);
% Extend the k-mesh to k_x in [-1, 1], k_y in [-1, 1].
% Mark the order of the EigenValue in the kpoints mat

% k_y < 0, mirror along x-axies
numIdx = app.EIGENVAL_3D.kpointsMesh(2, :) < 0;
app.EIGENVAL_3D.kpointsMesh = [app.EIGENVAL_3D.kpointsMesh(:, numIdx) + [0 1 0 0]', app.EIGENVAL_3D.kpointsMesh];
% k_y > 0, mirror along x-axies
numIdx = app.EIGENVAL_3D.kpointsMesh(2, :) > 0;
app.EIGENVAL_3D.kpointsMesh = [app.EIGENVAL_3D.kpointsMesh(:, numIdx) + [0 -1 0 0]', app.EIGENVAL_3D.kpointsMesh];
% K_X < 0, mirror along y-axies
numIdx = app.EIGENVAL_3D.kpointsMesh(1, :) < 0;
app.EIGENVAL_3D.kpointsMesh = [app.EIGENVAL_3D.kpointsMesh(:, numIdx) + [1 0 0 0]', app.EIGENVAL_3D.kpointsMesh];
% k_x > 0, mirror along y-axis
numIdx = app.EIGENVAL_3D.kpointsMesh(1, :) > 0;
app.EIGENVAL_3D.kpointsMesh = [app.EIGENVAL_3D.kpointsMesh(:, numIdx) + [-1 0 0 0]', app.EIGENVAL_3D.kpointsMesh];
% Sort the expanded result, index will be used in bandplot
app.EIGENVAL_3D.kpointsMesh = sortrows(app.EIGENVAL_3D.kpointsMesh');
app.EIGENVAL_3D.kpointsMesh = app.EIGENVAL_3D.kpointsMesh';
% K-mesh process (direct based to reciprocal based)
SCAL = str2double(app.POSCAR.stringMat(2,1));
a1 = str2double(app.POSCAR.stringMat(3, 1: 3))*SCAL;
a2 = str2double(app.POSCAR.stringMat(4, 1: 3))*SCAL;
a3 = str2double(app.POSCAR.stringMat(5, 1: 3))*SCAL;
V = dot(a1, cross(a2, a3));

b1 = 2*pi*cross(a2,a3)/V;
b2 = 2*pi*cross(a3,a1)/V;
b3 = 2*pi*cross(a1,a2)/V;
app.EIGENVAL_3D.transMat = [b1;b2;b3];
app.EIGENVAL_3D.kpointsRecipMesh = zeros(4, length(app.EIGENVAL_3D.kpointsMesh));
for i = 1: length(app.EIGENVAL_3D.kpointsRecipMesh)
    app.EIGENVAL_3D.kpointsRecipMesh(:, i) = [app.EIGENVAL_3D.kpointsMesh(1: 3, i)' * app.EIGENVAL_3D.transMat, app.EIGENVAL_3D.kpointsMesh(4, i)];
end
% Delete the duplicated elements
numIdx = false([1, length(app.EIGENVAL_3D.kpointsRecipMesh)]);
for i = 1: length(app.EIGENVAL_3D.kpointsRecipMesh) - 1
    if app.EIGENVAL_3D.kpointsRecipMesh(2, i) == app.EIGENVAL_3D.kpointsRecipMesh(2, i+1)
        numIdx(i) = 1;
    end
end
app.EIGENVAL_3D.kpointsRecipMesh(:, numIdx) = [];
% Set the restriction of the 1st Brillouin zone (Hex lattice only!!)
a = norm(a1);
xCriteria1 = abs(app.EIGENVAL_3D.kpointsRecipMesh(1, :)) < 2*pi/(3*a) ;
yCriteria1 = abs(app.EIGENVAL_3D.kpointsRecipMesh(2, :)) < 2*pi/(sqrt(3)*a);
xCriteria2 = ((app.EIGENVAL_3D.kpointsRecipMesh(1, :)) < -2*pi/(3*a)&((app.EIGENVAL_3D.kpointsRecipMesh(1, :)) >  -4*pi/(3*a)));
yCriteria2 = abs(app.EIGENVAL_3D.kpointsRecipMesh(2, :)) < sqrt(3)*app.EIGENVAL_3D.kpointsRecipMesh(1, :) + 4*pi/(sqrt(3)*a);
xCriteria3 = ((app.EIGENVAL_3D.kpointsRecipMesh(1, :)) > 2*pi/(3*a)&((app.EIGENVAL_3D.kpointsRecipMesh(1, :)) <  4*pi/(3*a)));
yCriteria3 = abs(app.EIGENVAL_3D.kpointsRecipMesh(2, :)) < -sqrt(3)*app.EIGENVAL_3D.kpointsRecipMesh(1, :) + 4*pi/(sqrt(3)*a);
numIdx = xCriteria1.*yCriteria1 + xCriteria2.*yCriteria2 + xCriteria3.*yCriteria3;
app.EIGENVAL_3D.brillouinZoneIndex = logical(numIdx);
% Set the trival z-value to zero for usage of selection panel,
% Process Fermi level.
numIdx = abs(app.EIGENVAL_3D.kpointsMesh(3, :)) < 1e-15;
app.EIGENVAL_3D.kpointsMesh(3, numIdx) = 0;
app.EIGENVAL_3D.kpointZvalue = unique(app.EIGENVAL_3D.kpointsMesh(3, :));
app.EIGENVAL_3D.bandEigenvalue = cellfun(@(x) ((x) - app.DOSCAR.fermiLevel + app.tmp.fermiLevelOffset), app.EIGENVAL_3D.bandEigenvalue);
% Get KzPlane Data
for numIdx = 1: length(app.EIGENVAL_3D.kpointZvalue)
    app.KzplaneDropDown.Items{numIdx} = num2str(app.EIGENVAL_3D.kpointZvalue(numIdx));
end
app.KzplaneDropDown.Value = app.KzplaneDropDown.Items(1);

%Plot data
% Data preparation
meshDensity = app.tmp.meshDensity;
kzPlaneSelection = app.EIGENVAL_3D.kzPlaneSelection;
kzPlaneSelection = kzPlaneSelection*app.EIGENVAL_3D.transMat;
kzPlaneSelection = kzPlaneSelection(3);
kpointsNumIdx = abs(app.EIGENVAL_3D.kpointsRecipMesh (3, :) - kzPlaneSelection) < 1e-6;
app.EIGENVAL_3D.kpointsRecipMeshSelected = app.EIGENVAL_3D.kpointsRecipMesh (:, kpointsNumIdx);
app.EIGENVAL_3D.bandeigenvalueSelected = cell(length(app.EIGENVAL_3D.bandIndexSelection), 1);
app.EIGENVAL_3D.interpolatedEigenvalue = cell(length(app.EIGENVAL_3D.bandIndexSelection), 1);
app.EIGENVAL_3D.atFermiLevelIndex = false(app.tmp.meshDensity);
% Plot 1st Brillouin Zone
figure();
brillouinZone = polyshape([2*pi/(3*a), 4*pi/(3*a), 2*pi/(3*a), -2*pi/(3*a), -4*pi/(3*a), -2*pi/(3*a)],...
    [2*pi/(sqrt(3)*a), 0, -2*pi/(sqrt(3)*a),  -2*pi/(sqrt(3)*a), 0, 2*pi/(sqrt(3)*a)]);
plot(brillouinZone, 'FaceAlpha', 0, 'EdgeColor', 'black', 'LineWidth', 4);
xlim([-1.5 1.5]);
ylim([-1.25 1.25]);
grid on;
title(['Fermi Surface with Fermi Level Offset = ', num2str(app.tmp.fermiLevelOffset)]);
% Determine the plot k-mesh boundary
kxMin = min(app.EIGENVAL_3D.kpointsRecipMeshSelected(1, :));
kxMax = max(app.EIGENVAL_3D.kpointsRecipMeshSelected(1, :));
kyMin = min(app.EIGENVAL_3D.kpointsRecipMeshSelected(2, :));
kyMax = max(app.EIGENVAL_3D.kpointsRecipMeshSelected(2, :));
% Particular data process for parallelization
for idx = 1: length(app.EIGENVAL_3D.bandIndexSelection)
    bandIndexSelection = app.EIGENVAL_3D.bandIndexSelection(idx);
    % Select the bandEigenValue from the index recorded in the
    % variable app.Eigenval_3D.kpointsRecipMesh(4, :)
    app.EIGENVAL_3D.bandeigenvalueSelected{idx} = zeros(1, length(app.EIGENVAL_3D.kpointsRecipMeshSelected));
    for numIdx = 1: length(app.EIGENVAL_3D.kpointsRecipMeshSelected)
        app.EIGENVAL_3D.bandeigenvalueSelected{idx}(numIdx) = app.EIGENVAL_3D.bandEigenvalue(bandIndexSelection, app.EIGENVAL_3D.kpointsRecipMeshSelected(4, numIdx));
    end
    % Projection to the 1st Brillouin Zone
    app.EIGENVAL_3D.bandeigenvalueSelected{idx}(~app.EIGENVAL_3D.brillouinZoneIndex) = NaN;
    
    % kpoints mesh creatation using VASP Kpoints data
    [app.EIGENVAL_3D.xmesh, app.EIGENVAL_3D.ymesh] = meshgrid(linspace(kxMin, kxMax, meshDensity), linspace(kyMin, kyMax, meshDensity));
    scatterModel = scatteredInterpolant(app.EIGENVAL_3D.kpointsRecipMeshSelected(1, :)', app.EIGENVAL_3D.kpointsRecipMeshSelected(2, :)', app.EIGENVAL_3D.bandeigenvalueSelected{idx}');
    scatterModel.Method = 'linear';
    app.EIGENVAL_3D.interpolatedEigenvalue{idx} = scatterModel(app.EIGENVAL_3D.xmesh, app.EIGENVAL_3D.ymesh);
    % Pick data at Fermi Level
%     app.EIGENVAL_3D.atFermiLevelIndex = app.EIGENVAL_3D.atFermiLevelIndex | (abs(app.EIGENVAL_3D.interpolatedEigenvalue{idx}) < 1e-2);    
    f = @(x, y) scatterModel(x, y);
    hold on;
    app.bandStructure_3D = fcontour(f, [kxMin, kxMax, kyMin, kyMax], '-r', 'LineWidth', 2, 'LevelList', [0.00 0.00], 'MeshDensity', app.tmp.meshDensity);
    hold off;
end

axis equal;
% Exports plot data
app.figureData.bandStructure_3D.bandIndexSelection = bandIndexSelection;
app.figureData.bandStructure_3D.kzPlaneSelection = kzPlaneSelection;
app.figureData.bandStructure_3D.kpointsRecipMeshSelected = app.EIGENVAL_3D.kpointsRecipMeshSelected;
app.figureData.bandStructure_3D.bandeigenvalueSelected = app.EIGENVAL_3D.bandeigenvalueSelected;
app.figureData.bandStructure_3D.xmesh = app.EIGENVAL_3D.xmesh;
app.figureData.bandStructure_3D.ymesh = app.EIGENVAL_3D.ymesh;
app.figureData.bandStructure_3D.interpolatedEigenvalue = app.EIGENVAL_3D.interpolatedEigenvalue;
app.figureData.bandStructure_3D.figure = app.bandStructure_3D;
% mesh(app.EIGENVAL_3D.xmesh, app.EIGENVAL_3D.ymesh, app.EIGENVAL_3D.interpolatedEigenvalue{6})

% Z-range settings
app.bandStructure_3D.Parent.XLabel.String = "K_x";
app.bandStructure_3D.Parent.YLabel.String = "K_y";
app.bandStructure_3D.Parent.ZLabel.String = "Energy (eV)";
app.ZminSpinner.Enable = "on";
app.ZmaxSpinner.Enable = "on";
app.tmp.bandStructure_3DZlimit = app.bandStructure_3D.Parent.ZLim;
app.ZminSpinner.Value = app.tmp.bandStructure_3DZlimit(1);
app.ZmaxSpinner.Value = app.tmp.bandStructure_3DZlimit(2);
% Enable Hubbard Model Fitting
app.OpenHubbardModelFitButton.Enable = "on";

%% Read POSCAR
clear variables;
% Parameter Setting
centerAtomIndex = 1;
% Read POSCAR
app.dirname = 'D:\tmp\test';
filename = strcat(app.dirname,'/POSCAR.vasp');
if ~exist(filename,'file')
    errordlg('POSCAR file not found','File Error');
    clear variables;
    return
end

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
app.POSCARtemp = readtable(filename, opts);
app.POSCARtemp = table2array(app.POSCARtemp);
app.POSCARtemp(ismissing(app.POSCARtemp)) = "";
for i = 1:length(app.POSCARtemp)
    tmp.numIdx = find( strlength(app.POSCARtemp(i, :)) );
    app.POSCAR.stringMat(i, 1: length(tmp.numIdx)) = app.POSCARtemp(i, tmp.numIdx);
end
app.POSCAR.stringMat(ismissing(app.POSCAR.stringMat)) = "";
app.POSCAR.element.name = app.POSCAR.stringMat(6, :);
app.POSCAR.element.name( app.POSCAR.element.name=='' )= [];
app.POSCAR.element.num = str2double (app.POSCAR.stringMat(7, :) );
app.POSCAR.element.num(isnan(app.POSCAR.element.num)) = [];
app.POSCAR.element.totalAtomNumber = sum(app.POSCAR.element.num);
tmp.totAtomNumber = sum(app.POSCAR.element.num);
clear opts i app.POSCARtemp
% Rotation
app.POSCAR.Position = str2double(app.POSCAR.stringMat(9: end, :));
% app.POSCAR.Position(4, :) = str2double(app.POSCAR.stringMat(12, :)) - str2double(app.POSCAR.stringMat(3, :)) - str2double(app.POSCAR.stringMat(4, :));
% app.POSCAR.Position(6, :) = str2double(app.POSCAR.stringMat(14, :)) - str2double(app.POSCAR.stringMat(3, :)) - str2double(app.POSCAR.stringMat(4, :));
translationMatrix = - app.POSCAR.Position(centerAtomIndex, :);
vec = str2double(app.POSCAR.stringMat(11, :)) ...
     - (str2double(app.POSCAR.stringMat(12, :)) - str2double(app.POSCAR.stringMat(3, :)) - str2double(app.POSCAR.stringMat(4, :)));
% phi = atan(vec(2)/vec(1));
% theta = atan(vec(3)/sqrt(vec(1).^2 + vec(2).^2)) + pi/2;
[phi, theta, r] = cart2sph(vec(1), vec(2), vec(3));
% Rotate along y-axis Definition
deg = theta;
rotMatY = [cos(deg) 0 sin(deg); 0 1 0; -sin(deg) 0 cos(deg)];
% Rotate along z-axis Definition
deg = -phi;
rotMatZ = [cos(deg) -sin(deg) 0; sin(deg) cos(deg) 0; 0 0 1];
% Rotation
app.POSCAR.rotated = zeros(length(app.POSCAR.Position), 3);
for numIdx = 1: length(app.POSCAR.Position)
    atomPos = app.POSCAR.Position(numIdx, :)' + translationMatrix';
    app.POSCAR.rotated(numIdx, :) = (rotMatY * rotMatZ * atomPos)';
end

%% Boat Test
clear variables;
a = 50;
theta = pi/6;
verticalTotLen = a/tan(theta);
totLen = a/sin(theta);
delta = 1;
Len = totLen;
i = 1;
while Len > a + delta
    Len = totLen - delta;
    verticalLen = sqrt(Len.^2 - a.^2);
    deltaX = verticalTotLen - verticalLen;
    verticalTotLen = verticalLen;
    result(i) = deltaX;
    i = i + 1;
    totLen = Len;
end
plot(1:length(result), result);
%% Temp Workspace
fig = figure();
hold on;
scatter(thetaMixUnique(sortedMixTypeUnique == 1), lengthMixUnique(sortedMixTypeUnique == 1), 5, 'filled');
scatter(thetaMixUnique(sortedMixTypeUnique == 2), lengthMixUnique(sortedMixTypeUnique == 2), 5, 'filled');
scatter(thetaMixUnique(sortedMixTypeUnique == 3), lengthMixUnique(sortedMixTypeUnique == 3), 5, 'filled');
scatter(thetaMixUnique(sortedMixTypeUnique == 4), lengthMixUnique(sortedMixTypeUnique == 4), 5, 'filled');
hold off;
legend({'case 1', 'case 2', 'case 3', 'case 4'});
%% Temp Workspace2
load('D:\OneDrive - tongji.edu.cn\Research\VSe2\Tphase\bulk\wannierFit\wannierFitBandNum11\wannier90-dos.dat');
fermiLevel = 4.55282876;
hold on;
plt = plot(wannier90_dos(:, 1) - fermiLevel, wannier90_dos(:, 2), 'r', 'LineWidth', 1.5, 'DisplayName', 'Wannier Fit');
hold off;
exportgraphics(fig, [figDir, '\', fileName, '.eps'], 'ContentType', 'vector');
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

