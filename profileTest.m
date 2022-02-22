%% Fermi Surface Plotting
clear variables;
% Control Flag Setting
% app.dirname = 'D:\OneDrive - tongji.edu.cn\Research\VSe2\Tphase\monolayer\nmAllBZ';
% app.dirname = 'D:\OneDrive - tongji.edu.cn\Research\VSe2\Tphase\bulk\nmOptVol';
% app.dirname = 'D:\OneDrive - tongji.edu.cn\Research\VSe2\Tphase\bilayer\Rotation\21.7867892982618\change_D_Ueff1\nm\2.12';
% app.dirname = 'D:\OneDrive - tongji.edu.cn\Research\VSe2\Tphase\bilayer\Rotation\21.7867892982618\change_D_Ueff1\nm\3.12AddEdos';
app.dirname = 'D:\OneDrive - tongji.edu.cn\Research\VSe2\Tphase\bilayer\nm\allBZ';
app.tmp.spinPolarizedFlag = 0;
app.EIGENVAL_3D.bandIndexSelection = 1:96;%100:200;%1:392;
app.tmp.meshDensity = 200;
app.EIGENVAL_3D.kzPlaneSelection = 0;
app.tmp.fermiLevelOffset = 0.0;
tmp.nedos = 2000;
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

figure;
plot(app.DOSCAR.dos.energy - app.DOSCAR.fermiLevel, app.DOSCAR.dos.DOS, '-x');
xlim([-8 8]);
grid on;
set(gca, 'XMinorGrid','on');
%% Data Process
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
app.EIGENVAL_3D.kpointsRecipMesh = sortrows(app.EIGENVAL_3D.kpointsRecipMesh')';
for i = 1: length(app.EIGENVAL_3D.kpointsRecipMesh) - 1
    if (abs(app.EIGENVAL_3D.kpointsRecipMesh(1, i) - app.EIGENVAL_3D.kpointsRecipMesh(1, i+1)) < 1e-8) && ...
       (abs(app.EIGENVAL_3D.kpointsRecipMesh(2, i) - app.EIGENVAL_3D.kpointsRecipMesh(2, i+1)) < 1e-8) && ...
       (abs(app.EIGENVAL_3D.kpointsRecipMesh(3, i) - app.EIGENVAL_3D.kpointsRecipMesh(3, i+1)) < 1e-8)
        numIdx(i) = 1;
    end
end
app.EIGENVAL_3D.kpointsRecipMesh(:, numIdx) = [];
% Set the restriction of the 1st Brillouin zone (Hex lattice only!!)
a = norm(a1);
pickTolerance = 1e-2;
polyEdgeX = [2*pi/(3*a) + pickTolerance, 4*pi/(3*a) + pickTolerance, 2*pi/(3*a) + pickTolerance,...
    -2*pi/(3*a) - pickTolerance, -4*pi/(3*a) - pickTolerance, -2*pi/(3*a) - pickTolerance, 2*pi/(3*a) + pickTolerance];
polyEdgeY = [2*pi/(sqrt(3)*a) + pickTolerance, 0, -2*pi/(sqrt(3)*a) - pickTolerance, -2*pi/(sqrt(3)*a) - pickTolerance,...
    0, 2*pi/(sqrt(3)*a) + pickTolerance, 2*pi/(sqrt(3)*a) + pickTolerance];
% xCriteria1 = abs(app.EIGENVAL_3D.kpointsRecipMesh(1, :)) < 2*pi/(3*a) + pickPrec;
% yCriteria1 = abs(app.EIGENVAL_3D.kpointsRecipMesh(2, :)) < 2*pi/(sqrt(3)*a) + pickPrec;
% xCriteria2 = ((app.EIGENVAL_3D.kpointsRecipMesh(1, :)) < -2*pi/(3*a) & ((app.EIGENVAL_3D.kpointsRecipMesh(1, :)) >  -4*pi/(3*a)) );
% yCriteria2 = (absapp.EIGENVAL_3D.kpointsRecipMesh(2, :)) < sqrt(3)*app.EIGENVAL_3D.kpointsRecipMesh(1, :) + 4*pi/(sqrt(3)*a);
% xCriteria3 = ((app.EIGENVAL_3D.kpointsRecipMesh(1, :) > 2*pi/(3*a))&((app.EIGENVAL_3D.kpointsRecipMesh(1, :)) < 4*pi/(3*a)));
% yCriteria3 = abs(app.EIGENVAL_3D.kpointsRecipMesh(2, :)) < -sqrt(3)*app.EIGENVAL_3D.kpointsRecipMesh(1, :) + 4*pi/(sqrt(3)*a);
% numIdx = xCriteria1.*yCriteria1 + xCriteria2.*yCriteria2 + xCriteria3.*yCriteria3;
% app.EIGENVAL_3D.brillouinZoneIndex = logical(numIdx);
app.EIGENVAL_3D.brillouinZoneIndex = inpolygon(app.EIGENVAL_3D.kpointsRecipMesh(1, :), app.EIGENVAL_3D.kpointsRecipMesh(2, :), polyEdgeX, polyEdgeY);
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
fig = figure();
brillouinZone = polyshape([2*pi/(3*a), 4*pi/(3*a), 2*pi/(3*a), -2*pi/(3*a), -4*pi/(3*a), -2*pi/(3*a)],...
    [2*pi/(sqrt(3)*a), 0, -2*pi/(sqrt(3)*a),  -2*pi/(sqrt(3)*a), 0, 2*pi/(sqrt(3)*a)]);
plot(brillouinZone, 'FaceAlpha', 0, 'EdgeColor', 'black', 'LineWidth', 4);
figAxis = fig.Children;
% xlim([-1.5 1.5]);
% ylim([-1.25 1.25]);
xlim([-0.75 0.75]);
ylim([-0.5 0.5]);
grid on;
title(['Fermi Surface with Fermi Level Offset = ', num2str(app.tmp.fermiLevelOffset)]);
% Determine the plot k-mesh boundary
kxMin = min(app.EIGENVAL_3D.kpointsRecipMeshSelected(1, :));
kxMax = max(app.EIGENVAL_3D.kpointsRecipMeshSelected(1, :));
kyMin = min(app.EIGENVAL_3D.kpointsRecipMeshSelected(2, :));
kyMax = max(app.EIGENVAL_3D.kpointsRecipMeshSelected(2, :));
% kpoints mesh creatation using VASP Kpoints data
[app.EIGENVAL_3D.xmesh, app.EIGENVAL_3D.ymesh] = meshgrid(linspace(kxMin, kxMax, meshDensity), linspace(kyMin, kyMax, meshDensity));
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
end
% Particular data process for parallelization
kpointsRecipMeshSelectedX = app.EIGENVAL_3D.kpointsRecipMeshSelected(1, :);
kpointsRecipMeshSelectedY = app.EIGENVAL_3D.kpointsRecipMeshSelected(2, :);
bandeigenvalueSelected = app.EIGENVAL_3D.bandeigenvalueSelected;
xmesh = app.EIGENVAL_3D.xmesh;
ymesh = app.EIGENVAL_3D.ymesh;
MeshDensity = app.tmp.meshDensity;
interpolatedEigenvalue = cell(1, length(app.EIGENVAL_3D.bandIndexSelection));
f = cell(1, length(app.EIGENVAL_3D.bandIndexSelection));

parfor idx = 1: length(app.EIGENVAL_3D.bandIndexSelection)
    scatterModel = scatteredInterpolant(kpointsRecipMeshSelectedX', kpointsRecipMeshSelectedY', bandeigenvalueSelected{idx}');
    scatterModel.Method = 'linear';
    interpolatedEigenvalue{idx} = scatterModel(xmesh, ymesh);
    % Pick data at Fermi Level
%     app.EIGENVAL_3D.atFermiLevelIndex = app.EIGENVAL_3D.atFermiLevelIndex | (abs(app.EIGENVAL_3D.interpolatedEigenvalue{idx}) < 1e-2);    
    f{idx} = @(x, y) scatterModel(x, y);
end
hold on;
for idx = 1: length(app.EIGENVAL_3D.bandIndexSelection)
    bandStructure_3D = fcontour(figAxis, f{idx}, [kxMin, kxMax, kyMin, kyMax], '-r', 'LineWidth', 2, 'LevelList', [0.00 0.00], 'MeshDensity', MeshDensity);
end
hold off;
clear kpointsRecipMeshSelectedX kpointsRecipMeshSelectedY bandeigenvalueSelected xmesh ymesh MeshDensity;
axis equal;
% Exports plot data
% app.figureData.bandStructure_3D.bandIndexSelection = bandIndexSelection;
% app.figureData.bandStructure_3D.kzPlaneSelection = kzPlaneSelection;
% app.figureData.bandStructure_3D.kpointsRecipMeshSelected = app.EIGENVAL_3D.kpointsRecipMeshSelected;
% app.figureData.bandStructure_3D.bandeigenvalueSelected = app.EIGENVAL_3D.bandeigenvalueSelected;
% app.figureData.bandStructure_3D.xmesh = app.EIGENVAL_3D.xmesh;
% app.figureData.bandStructure_3D.ymesh = app.EIGENVAL_3D.ymesh;
% app.figureData.bandStructure_3D.interpolatedEigenvalue = app.EIGENVAL_3D.interpolatedEigenvalue;
% app.figureData.bandStructure_3D.figure = app.bandStructure_3D;
% mesh(app.EIGENVAL_3D.xmesh, app.EIGENVAL_3D.ymesh, app.EIGENVAL_3D.interpolatedEigenvalue{6})

% Z-range settings
app.bandStructure_3D.Parent.XLabel.String = "K_x";
app.bandStructure_3D.Parent.YLabel.String = "K_y";
app.bandStructure_3D.Parent.ZLabel.String = "Energy (eV)";
% app.ZminSpinner.Enable = "on";
% app.ZmaxSpinner.Enable = "on";
% app.tmp.bandStructure_3DZlimit = app.bandStructure_3D.Parent.ZLim;
% app.ZminSpinner.Value = app.tmp.bandStructure_3DZlimit(1);
% app.ZmaxSpinner.Value = app.tmp.bandStructure_3DZlimit(2);
% Enable Hubbard Model Fitting
% app.OpenHubbardModelFitButton.Enable = "on";
