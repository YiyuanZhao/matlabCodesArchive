% clear variables;
a = 3.441;  % a-axies of the lattice
b = a;      % b-axies of the lattice
gamma = 120;% angle of <a, b>
sizeLattice = 51;  % Scale of the system (should be odd)

%% Initialize of the variables
centerOrder = (sizeLattice + 1) ./ 2;
lattice.x = zeros(sizeLattice);
lattice.y = zeros(sizeLattice);
lattice.distance = zeros(sizeLattice);
% Difference in x and y in the primitive cell
deltaA = [a, 0];
deltaB = [b*cosd(gamma), b*sind(gamma)];
% Set the zero of the axies
lattice.x(sizeLattice, 1) = 0;
lattice.y(sizeLattice, 1) = 0;
% Initialize the position of the center
lattice.x(centerOrder, centerOrder) = ((centerOrder - 1)*deltaA(1) + (sizeLattice - centerOrder)*deltaB(1));
lattice.y(centerOrder, centerOrder) = ((centerOrder - 1)*deltaA(2) + (sizeLattice - centerOrder)*deltaB(2));
%% Calculate the distance from the center
for row = sizeLattice: -1: 1
    for column = 1: sizeLattice
        lattice.x(row, column) = ((column - 1)*deltaA(1) + (sizeLattice - row)*deltaB(1));
        lattice.y(row, column) = ((column - 1)*deltaA(2) + (sizeLattice - row)*deltaB(2));
        lattice.hoppingA(row, column) = column - centerOrder;
        lattice.hoppingB(row, column) = centerOrder - row;
        lattice.distance(row, column) = sqrt((lattice.x(row, column) - lattice.x(centerOrder, centerOrder)).^2 + (lattice.y(row, column) - lattice.y(centerOrder, centerOrder)).^2);
    end
end
%% Get the N^(th) neighbour
output.distanceTemp = lattice.distance;
output.neighbour = cell(1, floor(numel(lattice.x)/6));
i = 1;
while min(min(output.distanceTemp)) ~= Inf
%     minValue = min(min(output.distanceTemp));
%     output.neighbour{i} = (output.distanceTemp - minValue) < 1e-6;
%     output.neighbourDistance{1, i} = output.distanceTemp(output.neighbour{i});
    [output.neighbour{i}, output.neighbourDistance, output.distanceTemp] = judgeNeighbour(lattice, output.neighbour{i}, output.distanceTemp, i, 12, 1e-1);
    output.neighbourDistance{2, i} = length(output.neighbourDistance{1, i});
    output.hopping{i}(:, 1) = lattice.hoppingA(output.neighbour{i});
    output.hopping{i}(:, 2) = lattice.hoppingB(output.neighbour{i});
    output.distanceTemp(output.neighbour{i}) = Inf;
    i = i+1;
end

%% Plot the neighour
color = lines(2^10);
fig = figure;
limit = [lattice.x(centerOrder, centerOrder) - 50, lattice.x(centerOrder, centerOrder) + 50, lattice.y(centerOrder, centerOrder) - 50, lattice.y(centerOrder, centerOrder) + 50];
scat = scatter(lattice.x(:), lattice.y(:), [], ones(length(lattice.x).^2, 3), 'filled');
axis(limit);
fig.PaperPositionMode = 'auto';
for i = 1: 22 %length(output.neighbour)
    scat.CData(output.neighbour{i}, :) = repmat(color(i, :), length(scat.CData(output.neighbour{i}, 1)), 1);
    %saveas(fig, ['C:\Users\SCES\Desktop\temp\img\', num2str(i),'.png']);
end

%% Temporary Codes
opts = delimitedTextImportOptions("NumVariables", 61);
opts.DataLines = [1, Inf];
opts.Delimiter = " ";
opts.VariableNames = ["written", "on", "Dec2020", "at", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var39", "Var40", "Var41", "Var42", "Var43", "Var44", "Var45", "Var46", "Var47", "Var48", "Var49", "Var50", "Var51", "Var52", "Var53", "Var54", "Var55", "Var56", "Var57", "Var58", "Var59", "Var60", "Var61"];
opts.SelectedVariableNames = ["written", "on", "Dec2020", "at", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2], "TrimNonNumeric", true);
opts = setvaropts(opts, [1, 2], "ThousandsSeparator", ",");
opts = setvaropts(opts, [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
grap = readtable("D:\OneDrive - tongji.edu.cn\vscode_workspace\DESKTOP-DQVLUVG\VScode_WorkSpace\model\data\grap.dat", opts);
tmp.grap = table2array(grap);
clear opts
% Preprocess of wanner_hr.dat, gaining degeneracy and hopping matrix
tmp.hoppingNumber = tmp.grap(3, 1);
tmp.residue = mod(tmp.hoppingNumber, 15);
tmp.loopMaxIndex = floor(tmp.hoppingNumber./15);
tmp.degeneracy = [reshape(tmp.grap(4: 3 + tmp.loopMaxIndex, :)', 1, []), tmp.grap(4 + tmp.loopMaxIndex, 1: tmp.residue)];
tmp.hopping = tmp.grap(tmp.loopMaxIndex + 5: end, 1: 7);
fileID = fopen('degeneracy.dat', 'w');
fprintf(fileID, "%5i", tmp.degeneracy);
fclose(fileID);
fileID = fopen('hoppingParameterProcessed.dat', 'w');
for numIdx = 1: length(tmp.hopping(:, 1))
    fprintf(fileID, "%5i %5i %5i %5i %5i %12.6f %12.6f\n", tmp.hopping(numIdx, :));
end
fclose(fileID);
%% External Functions
function [neighbour, neighbourDistance, distanceTemp] = judgeNeighbour(lattice, ~, distanceTemp, i, degeneracy, resolution)
minValue = min(min(distanceTemp));
neighbour = (distanceTemp - minValue) < 1e-6;
neighbourDistance{1, i} = distanceTemp(neighbour);
% Create N = degeneracy center-symmetric area in hex lattice
N = degeneracy;
degreeStep = 360 / N;
% Create boundary
xCenterOrder = (length(lattice.x(1, :)) + 1)/2;
yCenterOrder = (length(lattice.x(1, :)) + 1)/2;
xCenter = lattice.x(xCenterOrder, xCenterOrder);
yCenter = lattice.y(yCenterOrder, yCenterOrder);
xBound1 = [xCenter, xCenter];
xBound2 = [xCenter, max(max(lattice.x))];

yBoundMax1 = max(max(lattice.y));
yBoundMax2 = (xBound2(1) - xBound2(2)) / cosd(90 - degreeStep) + xBound2(2);
yBound = max([yBoundMax1, yBoundMax2]);
% Judge if the point is in/on the boundary
% [in, on] = inpolygon(lattice.x(neighbour), lattice.y(neighbour), [xBound1, xBound2], [yCenter, yBound, yBound, yCenter]);
griddata = griddedInterpolant(xBound2, [yCenter, yBound]);
xq = lattice.x(neighbour);
yq = griddata(xq);
onVertical = abs(lattice.y(neighbour) - yCenter) < resolution;
onBound2 = abs(lattice.y(neighbour) - yq) < resolution;
in = lattice.y(neighbour) - yq >= resolution;
end