clear variables;
%% Input parameters
% Lattice Construction
a1 = [3.3291   -0.0000    0.0000];
a2 = [-1.6645    2.8831   -0.0000];
a3 = [0.0000   -0.0000   23.1180];

% Hopping Matrix & Hopping Parameter Construction
hoppingMatrix{1} = [0, 0, 0];
hoppingMatrix{2} = [-1, 0, 0; -1, -1, 0; 0, 1, 0; 0, -1, 0; 1, 1, 0; 1, 0, 0];
hoppingMatrix{3} = [-2, -1, 0; -1, 1, 0; -1, -2, 0; 1, 2, 0; 1, -1, 0; 2, 1, 0];
hoppingMatrix{4} = [-2, 0, 0; -2, -2, 0; 0, 2, 0; 0, -2, 0; 2, 2, 0; 2, 0, 0];
hoppingParameter{1} = -2.430432;
hoppingParameter{2} = -0.040101;
hoppingParameter{3} =  0.097522;
hoppingParameter{4} = -0.070591;
%% Lattice data Processing
transMatA = [a1;a2;a3];
V = dot(a1, cross(a2, a3));
b1 = 2*pi*cross(a2,a3)/V;
b2 = 2*pi*cross(a3,a1)/V;
b3 = 2*pi*cross(a1,a2)/V;
transMatB = [b1;b2;b3];

%% Generates K-Path
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
%% OR import K-mesh
%  load('savedFigureData.mat');
%  selectionIndex = ~isnan(exportData.bandStructure_3D.bandeigenvalueSelected);
%  kPointsMesh = exportData.hubbardModel.kpointsRecipMeshSelected(1:3, selectionIndex);

%% Calculate H_kin
Hkin = zeros(1, length(kPointsMesh));
for loopIndex = 1:length(kPointsMesh)
    Hkin(loopIndex) = calculateKineticHamiltonian(transMatA, kPointsMesh(1:3, loopIndex)', hoppingParameter, hoppingMatrix);
end
 scatter3(kPointsMesh(1, :), kPointsMesh(2, :), Hkin, 1);
%  scatter(kPointsMesh(4, :), Hkin, 1);

%% External Function
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