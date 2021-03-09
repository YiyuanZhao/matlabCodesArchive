clear variables
workingMode = 2;                % 1: K-mesh at xy plane;  2: K-path at high-symmetry path
fermiLevel = -2.4474445;        % Should be defined in K-path Working Mode
fittingNeighbourOrder = 15;     % Fitting Order
fittingTypeDatabase = {'WannierFit', 'DFT'};     % options: 'DFT' / 'WannierFit'
calcDiffMethod = 2;
repeatTimes = 1;
optWeightMode = 0;
maxWeight = 2;
unityWidth = 1;

fittingNeighbourOrderMin = 5;
fittingNeighbourOrderStep = 1;
fittingNeighbourOrderMax = 5;

calcDiffMethodMin = 2;
calcDiffMethodMax = 2;

% numIdx = 1;
% for repeatTimes = 1: 5
%     for fittingNeighbourOrder = 5: 1: 5
%         initialHoppingParameter = generateInitialHoppingParameter(fittingNeighbourOrder, -0.1, 0.1);
%         for fittingTypeIndex = 2: length(fittingType)
%             for calcDiffMethod = 2: 2
%                 disp('fittingType   fittingNeighbourOrder   calcDiffMethod');
%                 disp({fittingType{fittingTypeIndex}, fittingNeighbourOrder, calcDiffMethod});
%                 [x{fittingTypeIndex, numIdx}, fval(fittingTypeIndex, numIdx), exitflag(fittingTypeIndex, numIdx), output{fittingTypeIndex, numIdx}] = ...
%                     HubbardFittingLiveFunction(workingMode, fermiLevel, fittingNeighbourOrder,...
%                     fittingType{fittingTypeIndex}, calcDiffMethod, initialHoppingParameter);
%                 numIdx = numIdx + 1;
%             end
%         end
%     end
%     save("result.mat");
% end

% x = cell(repeatTimes, 1);
% output = cell(repeatTimes, 1);
% fval = ones(repeatTimes, 1);
% exitflag = nan(repeatTimes, 1);
% for repeatTimesNumIdx = 1: repeatTimes
%     for fittingNeighbourOrder = fittingNeighbourOrderMin: fittingNeighbourOrderStep: fittingNeighbourOrderMax
%         initialHoppingParameter = generateInitialHoppingParameter(fittingNeighbourOrder, -0.1, 0.1);
%             for calcDiffMethod = calcDiffMethodMin: calcDiffMethodMax
%                 numIdx = (fittingNeighbourOrder - fittingNeighbourOrderMin) * (calcDiffMethodMax - calcDiffMethodMin + 1) ...
%                     + calcDiffMethod - 1
%                 disp('Repeat_Times   fittingNeighbourOrder   calcDiffMethod');
%                 disp({repeatTimesNumIdx, fittingNeighbourOrder, calcDiffMethod});
%                 [x{repeatTimesNumIdx, numIdx}, fval(repeatTimesNumIdx, numIdx), exitflag(repeatTimesNumIdx, numIdx), output{repeatTimesNumIdx, numIdx}] = ...
%                     HubbardFittingLiveFunction(workingMode, fermiLevel, fittingNeighbourOrder,...
%                     fittingType{2}, calcDiffMethod, initialHoppingParameter);
%             end
%     end
%     save("result.mat");
% end

initialHoppingParameter = generateInitialHoppingParameter(fittingNeighbourOrder, -0.1, 0.1);
fittingType = fittingTypeDatabase{2};
%% Function construction
% function [x, fval, exitflag, output] = HubbardFittingLiveFunction(workingMode, fermiLevel,...
%     fittingNeighbourOrder, fittingType, calcDiffMethod, initialHoppingParameter)
%% parameters Setting
% Start clock
initClock = clock;
initialHoppingParameterDataBase = [-2.430432, -0.040101, 0.097522, -0.070591, 0.0176250, -0.027414, 0.0046590, 0.001030, -0.010559, 0.000786, ...
    0.006749, -0.009858, 0.002722, -0.001447, 0.005350, -0.004718, 0.000434, 0.000423, 0.002563, 0.001689];
% if nargin == 5
%     initialHoppingParameter = initialHoppingParameterDataBase(1: fittingNeighbourOrder);
% end

%% Initialize fitting function (DO NOT EDIT THIS SECTION)
% Lattice Construction
a1 = [3.3291   -0.0000    0.0000];
a2 = [-1.6645    2.8831   -0.0000];
a3 = [0.0000   -0.0000   23.1180];
% Hopping Matrix & Hopping Parameter Construction
hoppingMatrix{1} = [0, 0, 0];
hoppingMatrix{2} = [-1, 0, 0; -1, -1, 0; 0, 1, 0; 0, -1, 0; 1, 1, 0; 1, 0, 0];
hoppingMatrix{3} = [-2, -1, 0; -1, 1, 0; -1, -2, 0; 1, 2, 0; 1, -1, 0; 2, 1, 0];
hoppingMatrix{4} = [-2, 0, 0; -2, -2, 0; 0, 2, 0; 0, -2, 0; 2, 2, 0; 2, 0, 0];
hoppingMatrix{5} = [-3,-1,0;-3,-2,0;-2,1,0;-2,-3,0;-1,2,0;-1,-3,0;1,3,0;1,-2,0;2,3,0;2,-1,0;3,2,0;3,1,0];
hoppingMatrix{6} = [-3,0,0;-3,-3,0;0,3,0;0,-3,0;3,3,0;3,0,0];
hoppingMatrix{7} = [-4,-2,0;-2,2,0;-2,-4,0;2,4,0;2,-2,0;4,2,0];
hoppingMatrix{8} = [-4,-1,0;-4,-3,0;-3,1,0;-3,-4,0;-1,3,0;-1,-4,0;1,4,0;1,-3,0;3,4,0;3,-1,0;4,3,0;4,1,0];
hoppingMatrix{9} = [-4,0,0;-4,-4,0;0,4,0;0,-4,0;4,4,0;4,0,0];
hoppingMatrix{10}= [-5,-2,0;-5,-3,0;-3,2,0;-3,-5,0;-2,3,0;-2,-5,0;2,5,0;2,-3,0;3,5,0;3,-2,0;5,3,0;5,2,0];
hoppingMatrix{11}= [-5,-1,0;-5,-4,0;-4,1,0;-4,-5,0;-1,4,0;-1,-5,0;1,5,0;1,-4,0;4,5,0;4,-1,0;5,4,0;5,1,0];
hoppingMatrix{12}=[-5,0,0;-5,-5,0;0,5,0;0,-5,0;5,5,0;5,0,0];
hoppingMatrix{13}=[-6,-3,0;-3,3,0;-3,-6,0;3,6,0;3,-3,0;6,3,0];
hoppingMatrix{14}=[-6,-2,0;-6,-4,0;-4,2,0;-4,-6,0;-2,4,0;-2,-6,0;2,6,0;2,-4,0;4,6,0;4,-2,0;6,4,0;6,2,0];
hoppingMatrix{15}=[-6,-1,0;-6,-5,0;-5,1,0;-5,-6,0;-1,5,0;-1,-6,0;1,6,0;1,-5,0;5,6,0;5,-1,0;6,5,0;6,1,0];
hoppingMatrix{16}=[-6,0,0;-6,-6,0;0,6,0;0,-6,0;6,6,0;6,0,0];
hoppingMatrix{17}=[-7,-3,0;-7,-4,0;-4,3,0;-4,-7,0;-3,4,0;-3,-7,0;3,7,0;3,-4,0;4,7,0;4,-3,0;7,4,0;7,3,0];
hoppingMatrix{18}=[-7,-2,0;-7,-5,0;-5,2,0;-5,-7,0;-2,5,0;-2,-7,0;2,7,0;2,-5,0;5,7,0;5,-2,0;7,5,0;7,2,0];
hoppingMatrix{19}=[-7,-1,0;-7,-6,0;-6,1,0;-6,-7,0;-1,6,0;-1,-7,0;1,7,0;1,-6,0;6,7,0;6,-1,0;7,6,0;7,1,0];
hoppingMatrix{20}=[-8,-4,0;-4,4,0;-4,-8,0;4,8,0;4,-4,0;8,4,0];
hoppingParameter{1} = -2.430432;
hoppingParameter{2} = -0.0401010000000000;
hoppingParameter{3} = 0.0975220000000000;
hoppingParameter{4} = -0.0705910000000000;
hoppingParameter{5} = 0.0176250000000000;
hoppingParameter{6} = -0.0274140000000000;
hoppingParameter{7} = 0.00465900000000000;
hoppingParameter{8} = 0.00103000000000000;
hoppingParameter{9} = -0.0105590000000000;
hoppingParameter{10}= 0.000786000000000000;
hoppingParameter{11}= 0.00674900000000000;
hoppingParameter{12}= -0.00985800000000000;
hoppingParameter{13}= 0.00272200000000000;
hoppingParameter{14}= -0.00144700000000000;
hoppingParameter{15}= 0.00535000000000000;
hoppingParameter{16}= -0.00471800000000000;
hoppingParameter{17}= 0.000434000000000000;
hoppingParameter{18}= 0.000423000000000000;
hoppingParameter{19}= 0.00256300000000000;
hoppingParameter{20}= 0.00168900000000000;

% Lattice data Processing
transMatA = [a1;a2;a3];
V = dot(a1, cross(a2, a3));
b1 = 2*pi*cross(a2,a3)/V;
b2 = 2*pi*cross(a3,a1)/V;
b3 = 2*pi*cross(a1,a2)/V;
transMatB = [b1;b2;b3];

if workingMode == 2
% Generates K-Path (Path Working Mode)
    kpoints{1} = [0 0 0]* transMatB;
    kpoints{2} = [0 0.5 0]* transMatB;
    kpoints{3} = [-0.3333 0.6667 0.0]* transMatB;
    kpoints{4} = [0 0 0]* transMatB;
    Ntot = 50;
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
    switch fittingType
        case 'WannierFit'
            load('wannier90_band.dat');
            % Interpolante the wannier fitting results for calculating difference.
            interpGrid = griddedInterpolant(wannier90_band(:, 1), wannier90_band(:, 2), 'linear');
            targetHkin = interpGrid(kPointsMesh(4, :)) - fermiLevel;
        case 'DFT'
            load('savedFigureData.mat');
            targetHkin = exportData.bandStructure.energy{12};            
    end
    
else
% OR import K-mesh (Plane Working Mode)
     load('savedFigureData.mat');
     selectionIndex = ~isnan(exportData.bandStructure_3D.bandeigenvalueSelected);
     kPointsMesh = exportData.hubbardModel.kpointsRecipMeshSelected(1:3, selectionIndex);
     targetHkin = exportData.hubbardModel.bandeigenvalueSelected(selectionIndex);
end

switch workingMode
    case 1
        f = @(t)calculateTargetFunction(workingMode, transMatA, kPointsMesh, ...
            t, hoppingMatrix, exportData.fermiLevel, targetHkin, calcDiffMethod,...
            optWeightMode, maxWeight, unityWidth);
    case 2
        f = @(t)calculateTargetFunction(workingMode, transMatA, kPointsMesh, ...
            t, hoppingMatrix, fermiLevel, targetHkin, calcDiffMethod,...
            optWeightMode, maxWeight, unityWidth);
end

% Create Problem, Linear Constrains, and Solution (Without constrain)
fminuncOptions = optimoptions('fminunc', "Algorithm","quasi-Newton", "MaxIterations", 1000,...
    "Display","iter",'FunctionTolerance',1e-9,...
    "MaxFunctionEvaluations", 4000, "UseParallel",true, "FiniteDifferenceType", 'central');
% fminsearchOptions = optimset('Display', 'iter', 'FunValCheck', 'on', 'PlotFcns', @optimplotfval);
[x,fval,exitflag,output,~,~] = fminunc(f, initialHoppingParameter, fminuncOptions);
% [x,fval,exitflag,output] = fminsearch(f, initialHoppingParameter, fminsearchOptions);

% Create Problem, Linear Constrains, and Solution (With constrain)

% load("constrain.mat");
% constrainNumber = 1;
% Aeq = [];
% for constrainIndex = 1: constrainNumber
%     Aeq = [Aeq; constrain{constrainIndex}(1: length(initialHoppingParameter))];
% end
% % Aeq = [1     6     6     6    12     6     6    12     6    12; ...
% %      1    -2    -2     6    -4    -2     6    -4     6    -4; ...
% %      1    -3     6    -3    -6     6     6    -6    -3    -6];
% % Aeq = [1     6     6     6    12     6     6    12     6    12];
% 
% switch constrainNumber
%     case 1
%         beq = real(targetHkin(1)) + fermiLevel;
%     case 3
%         beq = [real(targetHkin(1)), real(targetHkin(51)), real(targetHkin(101))] + fermiLevel;
% end
% A = [];
% b = [];
% lb = [];
% ub = [];
% nonlcon = [];
% fminconOptions = optimoptions("fmincon","Algorithm","interior-point","Display","iter",...
%     "PlotFcn","optimplotfval","UseParallel",true, "FiniteDifferenceType", 'central');
% [x,fval,exitflag,output] = fmincon(f, initialHoppingParameter, A, b, Aeq, beq, lb, ub, nonlcon, fminconOptions);

% Visualize the calulate result, uncomment if you need
updatedHoppingParameter = num2cell(x);
calcHkin = zeros(1, length(kPointsMesh));
initHkin = zeros(1, length(kPointsMesh));
parfor loopIndex = 1:length(kPointsMesh)
    calcHkin(loopIndex) = calculateKineticHamiltonian(transMatA, kPointsMesh(1:3, loopIndex)', updatedHoppingParameter, hoppingMatrix);
    initHkin(loopIndex) = calculateKineticHamiltonian(transMatA, kPointsMesh(1:3, loopIndex)', hoppingParameter, hoppingMatrix);
end

switch workingMode
    case 1
        switch optWeightMode
            case 0
                diff = calculateDifference(1, kPointsMesh(1, :)', fermiLevel ,calcHkin, targetHkin);
            case 1
                diff = calculateDifference(1, kPointsMesh(1, :)', fermiLevel ,calcHkin, targetHkin, maxWeight, unityWidth);
        end
        figure();
        scatter3(kPointsMesh(1, :), kPointsMesh(2, :), real(calcHkin) - exportData.fermiLevel, 1);
        figure();
        scatter3(kPointsMesh(1, :), kPointsMesh(2, :), real(diff), 1);
    case 2
        switch optWeightMode
            case 0
                diff = calculateDifference(1, kPointsMesh(4, :)', fermiLevel ,calcHkin, targetHkin);
            case 1
                diff = calculateDifference(1, kPointsMesh(4, :)', fermiLevel ,calcHkin, targetHkin, maxWeight, unityWidth);
        end
        figure();
        plot(kPointsMesh(4, :), real(calcHkin) - fermiLevel, 'b');
        hold on;
        plot(kPointsMesh(4, :), real(targetHkin), 'r');
        plot(kPointsMesh(4, :), real(initHkin) - fermiLevel, 'k--');
        legend(["Calculate", "Target", "Initial"]);
        hold off;
        figure();
        plot(kPointsMesh(4, :), diff);
end
FinClock = clock;
timeDuration = etime(FinClock, initClock);
disp(timeDuration);
% Export the updated hopping parameter
fileId = fopen(join(['.\summary\', join(string(FinClock(1: 5)),""), '.dat'], ""), "w");
for numIdx = 1: length(updatedHoppingParameter)
    fprintf(fileId, "%12.6f", updatedHoppingParameter{numIdx});
end
fclose(fileId);
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
    for NumIdx = 1: hoppingLength
        a1part{NumIdx} = hoppingMatrixProcessed{hoppingOrder}(NumIdx, 1).*transMat(1, :);
        a2part{NumIdx} = hoppingMatrixProcessed{hoppingOrder}(NumIdx, 2).*transMat(2, :);
        a3part{NumIdx} = hoppingMatrixProcessed{hoppingOrder}(NumIdx, 3).*transMat(3, :);
        dotpart(NumIdx) = dot(kpointsRecip, a1part{NumIdx} + a2part{NumIdx} + a3part{NumIdx});
        sum = sum + exp(1i*dotpart(NumIdx));
    end
    Hkin = Hkin + hoppingParameter{hoppingOrder}*sum;
end
end

% This calculateDifference function can be customized depending on your
% fitting error method.
function diff = calculateDifference(workingMode, kpointsRecip, fermiLevel ,calcHkin, targetHkin, maxWeight, unityWidth)
difference = zeros(1, length(kpointsRecip));

if nargin == 7
    weight = calculateWeightProportialToEnergy(maxWeight, unityWidth, targetHkin, 0);
    for LoopIndex = 1: length(kpointsRecip)
        difference(LoopIndex) = calcHkin(LoopIndex) - fermiLevel - targetHkin(LoopIndex);
    end
    difference = difference.*weight;
else
    for LoopIndex = 1: length(kpointsRecip)
        difference(LoopIndex) = calcHkin(LoopIndex) - fermiLevel - targetHkin(LoopIndex);
    end
end

switch workingMode
    case 1
        diff = real(difference);
    case 2
        [~, maxIdx] = max(abs(difference));
        diff = real(difference(maxIdx));
    case 3
        diff = mean(real(difference));
    case 4
        diff = var(real(difference));
end

end

function fittingDifference = calculateTargetFunction(workingMode, transMatA, kPointsMesh, hoppingParameter,...
    hoppingMatrix, fermiLevel, targetHkin, calcDiffMethod, optWeightMode, maxWeight, unityWidth)

hoppingParameterCell = num2cell(hoppingParameter);
Hkin = zeros(1, length(kPointsMesh));
for LoopIndex = 1:length(kPointsMesh)
    Hkin(LoopIndex) = calculateKineticHamiltonian(transMatA, kPointsMesh(1:3, LoopIndex)', hoppingParameterCell, hoppingMatrix);
end
switch workingMode
    case 1
        switch optWeightMode
            case 0
                Diff = calculateDifference(calcDiffMethod, kPointsMesh(1, :)', fermiLevel, Hkin, targetHkin);
            case 1
                Diff = calculateDifference(calcDiffMethod, kPointsMesh(1, :)', fermiLevel, Hkin, targetHkin, maxWeight, unityWidth);
        end
        fittingDifference = abs(Diff);
    case 2
        switch optWeightMode
            case 0
                Diff = calculateDifference(calcDiffMethod, kPointsMesh(4, :)', fermiLevel, Hkin, targetHkin);
            case 1
                Diff = calculateDifference(calcDiffMethod, kPointsMesh(4, :)', fermiLevel, Hkin, targetHkin, maxWeight, unityWidth);
        end
        fittingDifference = abs(Diff);
end
end

function initParam = generateInitialHoppingParameter(order, paraMin, paraMax)
    initParam = zeros(1, order);
    initParam(1) = -2.40;
    randMatrix = paraMin + (paraMax - paraMin) * rand(1, order - 1);
    initParam(2:end) = randMatrix;
end

function weight = calculateWeightProportialToEnergy(maxWeight, unityWidth, targetHkin, FermiLevel)
    % Gauss distribution, peaks at Fermilevel with input variable maxWeight,
    % and decays to unity at unityWidth.
    % Distribution Density: $$ y = \frac{1}{\sqrt{2\pi}\sigma} exp{[-\frac{(E_{tar} - E_f)^2}{k \cdot 2\sigma^2}]}$$
    % where k is a renormalization parameter which satisfies $$ k = \frac{\pi L^2 N^2}{ln N}$$, others are standard gauss
    % distribution probability density function(pdf).
    sigma = 1/(maxWeight * sqrt(2*pi()));
    k = pi()*unityWidth^2*maxWeight^2/log(maxWeight);
    weight = zeros(1, length(targetHkin));
    for numIdx = 1: length(targetHkin)
        weight(numIdx) = 1/(sqrt(2*pi())*sigma)*exp(-(targetHkin(numIdx)-FermiLevel)^2/(k*2*sigma^2));
        if weight(numIdx) < 1
            weight(numIdx) = 1;
        end
    end
    end