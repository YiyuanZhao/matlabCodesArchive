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
%%
a.a1 = 1;
a.a2 = 2;
a.a3 = 3;
out = geta(a);
%%
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

function out = geta(a)
out = a.a1;
end
