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
workMode = 'afm';
atomNumber = length(mag(:, 1))/3;
switch workMode
    case 'fm'
        magmom_V  = mean(mag(1: atomNumber, 5))
        magmom_Se = mean(mag(atomNumber + 1: end, 5))
        var_V     = var(mag(1: atomNumber, 5))
    case 'afm'
        magmom_V  = mean(abs(mag(1: atomNumber, 5)))
        magmom_Se = mean(abs(mag(atomNumber + 1: end, 5)))
        var_V     = var(abs(mag(1: atomNumber, 5)))
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
