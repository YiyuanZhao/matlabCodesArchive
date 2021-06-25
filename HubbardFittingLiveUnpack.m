clear variables
workingMode = 1;                % 1: K-mesh at xy plane;  2: K-path at high-symmetry path
fermiLevel = -2.4474445;        % Should be defined in K-path Working Mode
fittingNeighbourOrder = 30;     % Fitting Order
fittingTypeDatabase = {'WannierFit', 'DFT'};     % options: 'DFT' / 'WannierFit'
fittingType = fittingTypeDatabase{2};
calcDiffMethod = 5;
repeatTimes = 1;
continueOptFlag = 0;
optWeightOptions.optWeightMode = 0;
optWeightOptions.maxWeight = 2;
optWeightOptions.unityWidth = 1;
% rng(5);

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
for repNumIdx = 1: repeatTimes
initialHoppingParameter = generateInitialHoppingParameter(fittingNeighbourOrder, -0.1, 0.1);
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
if continueOptFlag == 1
    load('updatedHoppingParameter.mat');
    initialHoppingParameter = cell2mat(updatedHoppingParameter);
end
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
hoppingMatrix{21}=[-7,0,0;-7,-7,0;0,7,0;0,-7,0;7,7,0;7,0,0];
hoppingMatrix{22}=[-8,-3,0;-8,-5,0;-5,3,0;-5,-8,0;-3,5,0;-3,-8,0;3,8,0;3,-5,0;5,8,0;5,-3,0;8,5,0;8,3,0];
hoppingMatrix{23}=[-8,-2,0;-8,-6,0;-6,2,0;-6,-8,0;-2,6,0;-2,-8,0;2,8,0;2,-6,0;6,8,0;6,-2,0;8,6,0;8,2,0];
hoppingMatrix{24}=[-8,-1,0;-8,-7,0;-7,1,0;-7,-8,0;-1,7,0;-1,-8,0;1,8,0;1,-7,0;7,8,0;7,-1,0;8,7,0;8,1,0];
hoppingMatrix{25}=[-9,-4,0;-9,-5,0;-5,4,0;-5,-9,0;-4,5,0;-4,-9,0;4,9,0;4,-5,0;5,9,0;5,-4,0;9,5,0;9,4,0];
hoppingMatrix{26}=[-9,-3,0;-9,-6,0;-6,3,0;-6,-9,0;-3,6,0;-3,-9,0;3,9,0;3,-6,0;6,9,0;6,-3,0;9,6,0;9,3,0];
hoppingMatrix{27}=[-8,0,0;-8,-8,0;0,8,0;0,-8,0;8,8,0;8,0,0];
hoppingMatrix{28}=[-9,-2,0;-9,-7,0;-7,2,0;-7,-9,0;-2,7,0;-2,-9,0;2,9,0;2,-7,0;7,9,0;7,-2,0;9,7,0;9,2,0];
hoppingMatrix{29}=[-9,-1,0;-9,-8,0;-8,1,0;-8,-9,0;-1,8,0;-1,-9,0;1,9,0;1,-8,0;8,9,0;8,-1,0;9,8,0;9,1,0];
hoppingMatrix{30}=[-10,-5,0;-5,5,0;-5,-10,0;5,10,0;5,-5,0;10,5,0];
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
            optWeightOptions);
    case 2
        f = @(t)calculateTargetFunction(workingMode, transMatA, kPointsMesh, ...
            t, hoppingMatrix, fermiLevel, targetHkin, calcDiffMethod,...
            optWeightOptions);
end

%% Create Problem, Linear Constrains, and Solution (Without constrain)
% fminuncOptions = optimoptions('fminunc', "Algorithm","quasi-Newton", "MaxIterations", 1000,...
%     "Display","iter-detailed",'FunctionTolerance',1e-9,...
%     "MaxFunctionEvaluations", 40000, "UseParallel",true, ...
%     "FiniteDifferenceType", 'central', 'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-7);
% % fminsearchOptions = optimset('Display', 'iter', 'FunValCheck', 'on', 'PlotFcns', @optimplotfval);
% [x,fval,exitflag,output,~,~] = fminunc(f, initialHoppingParameter, fminuncOptions);
% % [x,fval,exitflag,output] = fminsearch(f, initialHoppingParameter, fminsearchOptions);

%% Or least-square-question fitting
% lsqnonlinOptions = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt', 'Display', 'iter', 'UseParallel', true, "MaxIterations", 1000,...
%     "Display", "iter", "MaxFunctionEvaluations", 40000);
% [x,resnorm,residual,exitflag,output] = lsqnonlin(f, initialHoppingParameter,[], [], lsqnonlinOptions);
%% Or PSO algorithm
particleswarmOptions = optimoptions('particleswarm', 'MaxIterations', 1000, 'UseParallel', true, ...
    'Display', 'iter', 'OutputFcn', @pswplotranges, 'SwarmSize', 1000);
lb = [-5, -ones(1, fittingNeighbourOrder - 1)];
ub = [ 5,  ones(1, fittingNeighbourOrder - 1)];
[x,~,~,~] = particleswarm(f, fittingNeighbourOrder, lb, ub, particleswarmOptions);
updatedHoppingParameter = x;
fminuncOptions = optimoptions('fminunc', "Algorithm","quasi-Newton", "MaxIterations", 1000,...
    "Display","iter-detailed",'FunctionTolerance',1e-9,...
    "MaxFunctionEvaluations", 40000, "UseParallel",true, ...
    "FiniteDifferenceType", 'central', 'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-7);
[x,fval,exitflag,output,~,~] = fminunc(f, updatedHoppingParameter, fminuncOptions);
%% Create Problem, Linear Constrains, and Solution (With constrain)

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

%% Visualize the calulate result, uncomment if you need
updatedHoppingParameter = num2cell(x);
calcHkin = zeros(1, length(kPointsMesh));
initHkin = zeros(1, length(kPointsMesh));
kpointsMesh1 = kPointsMesh(1:3, :);
kpointsMesh2 = kPointsMesh(1:3, :);
parfor loopIndex = 1:length(kPointsMesh)
    calcHkin(loopIndex) = calculateKineticHamiltonian(transMatA, kpointsMesh1(:, loopIndex)', updatedHoppingParameter, hoppingMatrix);
    initHkin(loopIndex) = calculateKineticHamiltonian(transMatA, kpointsMesh2(:, loopIndex)', hoppingParameter, hoppingMatrix);
end
switch workingMode
    case 1
        switch optWeightOptions.optWeightMode
            case 0
                diff = calculateDifference(1, kPointsMesh(1, :)', fermiLevel ,calcHkin, targetHkin);
            case 1
                diff = calculateDifference(1, kPointsMesh(1, :)', fermiLevel ,calcHkin, targetHkin, optWeightOptions);
        end
        figure();
        scatter3(kPointsMesh(1, :), kPointsMesh(2, :), real(calcHkin) - exportData.fermiLevel, 1);
        figure();
        scatter3(kPointsMesh(1, :), kPointsMesh(2, :), real(diff), 1);
    case 2
        switch optWeightOptions.optWeightMode
            case 0
                diff = calculateDifference(1, kPointsMesh(4, :)', fermiLevel ,calcHkin, targetHkin);
            case 1
                diff = calculateDifference(1, kPointsMesh(4, :)', fermiLevel ,calcHkin, targetHkin, optWeightOptions);
        end
        figure();
        plot(kPointsMesh(4, :), real(calcHkin) - fermiLevel, 'b');
        hold on;
        plot(kPointsMesh(4, :), real(targetHkin), 'r');
        legend(["Calculate", "Target"]);
%         plot(kPointsMesh(4, :), real(initHkin) - fermiLevel, 'k--');
%         legend(["Calculate", "Target", "Initial"]);
        hold off;
%         figure();
%         plot(kPointsMesh(4, :), diff);
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
save('updatedHoppingParameter.mat', 'updatedHoppingParameter');
end
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
function diff = calculateDifference(workingMode, kpointsRecip, fermiLevel ,calcHkin, targetHkin, optWeightOptions)
difference = zeros(1, length(kpointsRecip));

if nargin == 6
    weight = calculateWeightProportialToEnergy(optWeightOptions.maxWeight, optWeightOptions.unityWidth, targetHkin, 0);
    for LoopIndex = 1: length(kpointsRecip)
        difference(LoopIndex) = calcHkin(LoopIndex) - fermiLevel - targetHkin(LoopIndex);
    end
    difference = difference.*weight;
else
    parfor LoopIndex = 1: length(kpointsRecip)
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
    case 5
        diff = sum(difference.^2)./length(difference);
end

end

function fittingDifference = calculateTargetFunction(workingMode, transMatA, kPointsMesh, hoppingParameter,...
    hoppingMatrix, fermiLevel, targetHkin, calcDiffMethod, optWeightOptions)

hoppingParameterCell = num2cell(hoppingParameter);
Hkin = zeros(1, length(kPointsMesh));

switch workingMode
    case 1
        kPointsMesh = kPointsMesh(1:3, :);
        parfor LoopIndex = 1:length(kPointsMesh)
            Hkin(LoopIndex) = calculateKineticHamiltonian(transMatA, kPointsMesh(:, LoopIndex)', hoppingParameterCell, hoppingMatrix);
        end
        switch optWeightOptions.optWeightMode
            case 0
                Diff = calculateDifference(calcDiffMethod, kPointsMesh(1, :)', fermiLevel, Hkin, targetHkin);
            case 1
                Diff = calculateDifference(calcDiffMethod, kPointsMesh(1, :)', fermiLevel, Hkin, targetHkin, optWeightOptions);
        end
        fittingDifference = abs(Diff);
    case 2
        kpointsPath = kPointsMesh(4, :);
        kPointsMesh = kPointsMesh(1:3, :);
        parfor LoopIndex = 1:length(kPointsMesh)
            Hkin(LoopIndex) = calculateKineticHamiltonian(transMatA, kPointsMesh(:, LoopIndex)', hoppingParameterCell, hoppingMatrix);
        end
        switch optWeightOptions.optWeightMode
            case 0
                Diff = calculateDifference(calcDiffMethod, kpointsPath', fermiLevel, Hkin, targetHkin);
            case 1
                Diff = calculateDifference(calcDiffMethod, kpointsPath', fermiLevel, Hkin, targetHkin, optWeightOptions);
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

function stop = pswplotranges(optimValues,state)
stop = false; % This function does not stop the solver
switch state
    case 'init'
        nplot = size(optimValues.swarm,2); % Number of dimensions
        if nplot <= 6
            for i = 1:nplot % Set up axes for plot
                subplot(nplot,1,i);
                tag = sprintf('psoplotrange_var_%g',i); % Set a tag for the subplot
                semilogy(optimValues.iteration,0,'-k','Tag',tag); % Log-scaled plot
                ylabel(num2str(i))
                xlabel('Iteration','interp','none'); % Iteration number at the bottom
                subplot(nplot,1,1) % Title at the top
            end
        else
            rowNum = ceil(sqrt(nplot));
            colNum = ceil(nplot/rowNum);
            for i = 1:nplot % Set up axes for plot
                subplot(rowNum, colNum, i);
                tag = sprintf('psoplotrange_var_%g',i); % Set a tag for the subplot
                semilogy(optimValues.iteration,0,'-k','Tag',tag); % Log-scaled plot
                ylabel(num2str(i))
                xlabel('Iteration','interp','none'); % Iteration number at the bottom
                subplot(rowNum,colNum,1) % Title at the top
            end
        end
        
        title('Log range of particles by component')
        setappdata(gcf,'t0',tic); % Set up a timer to plot only when needed
    case 'iter'
        nplot = size(optimValues.swarm,2); % Number of dimensions
        if nplot <= 6
            for i = 1:nplot
                subplot(nplot,1,i);
                % Calculate the range of the particles at dimension i
                irange = max(optimValues.swarm(:,i)) - min(optimValues.swarm(:,i));
                tag = sprintf('psoplotrange_var_%g',i);
                plotHandle = findobj(get(gca,'Children'),'Tag',tag); % Get the subplot
                xdata = plotHandle.XData; % Get the X data from the plot
                newX = [xdata optimValues.iteration]; % Add the new iteration
                plotHandle.XData = newX; % Put the X data into the plot
                ydata = plotHandle.YData; % Get the Y data from the plot
                newY = [ydata irange]; % Add the new value
                plotHandle.YData = newY; % Put the Y data into the plot
            end
        else
            for i = 1:nplot
                rowNum = ceil(sqrt(nplot));
                colNum = ceil(nplot/rowNum);
                subplot(rowNum, colNum, i);
                % Calculate the range of the particles at dimension i
                irange = max(optimValues.swarm(:,i)) - min(optimValues.swarm(:,i));
                tag = sprintf('psoplotrange_var_%g',i);
                plotHandle = findobj(get(gca,'Children'),'Tag',tag); % Get the subplot
                xdata = plotHandle.XData; % Get the X data from the plot
                newX = [xdata optimValues.iteration]; % Add the new iteration
                plotHandle.XData = newX; % Put the X data into the plot
                ydata = plotHandle.YData; % Get the Y data from the plot
                newY = [ydata irange]; % Add the new value
                plotHandle.YData = newY; % Put the Y data into the plot
            end
        end
        if toc(getappdata(gcf,'t0')) > 1/30 % If 1/30 s has passed
            drawnow % Show the plot
            setappdata(gcf,'t0',tic); % Reset the timer
        end
    case 'done'
        % No cleanup necessary
end
end