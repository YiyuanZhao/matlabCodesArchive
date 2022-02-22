%% Read Data from MySQL database
clear variables;
% Parameter Setting: material, cellType, layerDistance
conn = database('local_nm_dos','','');
material = '''VSe2''';
% cellType = '''Primary''';
cellType = '''Super''';
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
        query = ['select * from nm_DOSCAR_isym0 ' ...
            'where  material = ', material, ' and cellType = ', cellType, ' and '...
            'magConfiguration = ',magConfiguration, ' and atomIdx = ', num2str(atomNumber), ...
            ' and round(layerDistance, 2) = ', num2str(layerDistance) ,' ;'];
        nmdata{layerIdx, atomNumber} = fetch(conn, query);
        dispStr = ['Fetched ', num2str(atomNumber), '/', num2str(totAtomNumber)];
        disp(dispStr);
    end
end
close(conn);

%% Calculate DOS accumulation
% Get number of V atoms
if strcmp(cellType,'''Primary''')
    totAtomNumber = 6;
    totVatomNumber = 2;
else
    totAtomNumber = 42;
    totVatomNumber = 14;
end
% Visualize occupation
occupiedDOS = cell(2, totAtomNumber);
for layerIdx = 1: totLayerDistanceNumber
    layerDistance = layerDistanceDatabase(layerIdx);
    closestFermiEnergy = min(abs(nmdata{layerIdx, 1}.energy));
    energyIdx = find(abs(nmdata{layerIdx, 1}.energy) == closestFermiEnergy);
        for atomNumber = 1: totAtomNumber
            occupiedDOS{layerIdx, atomNumber}.s_int = sum(nmdata{layerIdx, atomNumber}.s(1: energyIdx));
            occupiedDOS{layerIdx, atomNumber}.p_y_int = sum(nmdata{layerIdx, atomNumber}.p_y(1: energyIdx));
            occupiedDOS{layerIdx, atomNumber}.p_z_int = sum(nmdata{layerIdx, atomNumber}.p_z(1: energyIdx));
            occupiedDOS{layerIdx, atomNumber}.p_x_int = sum(nmdata{layerIdx, atomNumber}.p_x(1: energyIdx));
            occupiedDOS{layerIdx, atomNumber}.d_xy_int = sum(nmdata{layerIdx, atomNumber}.d_xy(1: energyIdx));
            occupiedDOS{layerIdx, atomNumber}.d_yz_int = sum(nmdata{layerIdx, atomNumber}.d_yz(1: energyIdx));
            occupiedDOS{layerIdx, atomNumber}.d_x2_r2_int = sum(nmdata{layerIdx, atomNumber}.d_x2_r2(1: energyIdx));
            occupiedDOS{layerIdx, atomNumber}.d_xz_int = sum(nmdata{layerIdx, atomNumber}.d_xz(1: energyIdx));
            occupiedDOS{layerIdx, atomNumber}.d_x2_y2_int = sum(nmdata{layerIdx, atomNumber}.d_x2_y2(1: energyIdx));
            tab = struct2table(occupiedDOS{layerIdx, atomNumber});
            occupiedDOS{layerIdx, atomNumber}.sum_int = sum(tab{1, :});
            occupiedDOS{layerIdx, atomNumber}.sum_int_tot = sum(sum(table2array(nmdata{layerIdx, atomNumber}(:, 7:15))));
        end
end
for layerIdx = 1: length(nmdata(:, 1))
    % Visualize V atom Data
    figure();
    mat = cellfun(@(x) x.s_int, occupiedDOS);
    subplot(3, 4, 1);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('s-int');
    grid on;
    xlim([1 totVatomNumber]);
    mat = cellfun(@(x) x.p_y_int, occupiedDOS);
    subplot(3, 4, 2);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('p-y-int');
    grid on;
    xlim([1 totVatomNumber]);
    mat = cellfun(@(x) x.p_z_int, occupiedDOS);
    subplot(3, 4, 3);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('p-z-int');
    grid on;
    xlim([1 totVatomNumber]);
    mat = cellfun(@(x) x.p_x_int, occupiedDOS);
    subplot(3, 4, 4);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('p-x-int');
    grid on;
    xlim([1 totVatomNumber]);
    mat = cellfun(@(x) x.d_xy_int, occupiedDOS);
    subplot(3, 4, 5);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('d-xy-int');
    grid on;
    xlim([1 totVatomNumber]);
    mat = cellfun(@(x) x.d_yz_int, occupiedDOS);
    subplot(3, 4, 6);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('d-yz-int');
    grid on;
    xlim([1 totVatomNumber]);
    mat = cellfun(@(x) x.d_x2_r2_int, occupiedDOS);
    subplot(3, 4, 7);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('d-x2-r2-int');
    grid on;
    xlim([1 totVatomNumber]);
    mat = cellfun(@(x) x.d_xz_int, occupiedDOS);
    subplot(3, 4, 8);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('d-xz-int');
    grid on;
    xlim([1 totVatomNumber]);
    mat = cellfun(@(x) x.d_x2_y2_int, occupiedDOS);
    subplot(3, 4, 9);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('d-x2-y2-int');
    grid on;
    xlim([1 totVatomNumber]);
    mat = cellfun(@(x) x.sum_int, occupiedDOS);
    subplot(3, 4, 10);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('sum-int');
    grid on;
    xlim([1 totVatomNumber]);
    mat = cellfun(@(x) x.sum_int_tot, occupiedDOS);
    subplot(3, 4, 11);
    plot(1:totVatomNumber, mat(layerIdx, 1: totVatomNumber));
    title('sum-int-tot');
    grid on;
    xlim([1 totVatomNumber]);
    % Visualize Se atom Data
    figure();
    mat = cellfun(@(x) x.s_int, occupiedDOS);
    subplot(3, 4, 1);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('s-int');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
    mat = cellfun(@(x) x.p_y_int, occupiedDOS);
    subplot(3, 4, 2);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('p-y-int');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
    mat = cellfun(@(x) x.p_z_int, occupiedDOS);
    subplot(3, 4, 3);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('p-z-int');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
    mat = cellfun(@(x) x.p_x_int, occupiedDOS);
    subplot(3, 4, 4);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('p-x-int');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
    mat = cellfun(@(x) x.d_xy_int, occupiedDOS);
    subplot(3, 4, 5);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('d-xy-int');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
    mat = cellfun(@(x) x.d_yz_int, occupiedDOS);
    subplot(3, 4, 6);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('d-yz-int');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
    mat = cellfun(@(x) x.d_x2_r2_int, occupiedDOS);
    subplot(3, 4, 7);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('d-x2-r2-int');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
    mat = cellfun(@(x) x.d_xz_int, occupiedDOS);
    subplot(3, 4, 8);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('d-xz-int');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
    mat = cellfun(@(x) x.d_x2_y2_int, occupiedDOS);
    subplot(3, 4, 9);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('d-x2-y2-int');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
    mat = cellfun(@(x) x.sum_int, occupiedDOS);
    subplot(3, 4, 10);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('sum-int');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
    mat = cellfun(@(x) x.sum_int_tot, occupiedDOS);
    subplot(3, 4, 11);
    plot(totVatomNumber + 1: atomNumber, mat(layerIdx, totVatomNumber + 1 : atomNumber));
    title('sum-int-tot');
    grid on;
    xlim([totVatomNumber + 1, atomNumber]);
end
%% Process total dos
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
    nmLine = plot(totalDosNmAverage{layerIdx, 1}.energy, totalDosNmAverage{layerIdx, 1}.tot, 'Color', 'black');
    xlim([-6 6]);
    plot([0 0], [-150 150], '--', 'Color', [0.25 0.25 0.25]);
    hold off;
    grid on;
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
visualPickAtomV = [1 2 4];
visualPickAtomSe = [15 16 17 23];
pdosVatomVerticalLineLim = [-8 8];
pdosVatomSeerticalLineLim = [-2 2];
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
            plot([0 0], pdosVatomVerticalLineLim, '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            grid on;
            legend([fmLineUp afmLineUp nmLine], {'FM', 'AFM', 'NM'});
            title({[num2str(orbitalPickAtomV(orbitPickAtomVIdx)), ' DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: V, Seq.Num ',num2str(visualPickAtomV(pickAtomVIdx))]});
            set(gca,'XMinorGrid','on')
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
            plot([0 0], pdosVatomSeerticalLineLim, '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            grid on;
            legend([fmLineUp afmLineUp nmLine], {'FM', 'AFM', 'NM'});
            title({[num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ' DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: Se'];...
                ['Seq.Num ',num2str(visualPickAtomSe(pickAtomSeIdx))]});
            set(gca,'XMinorGrid','on')
            idx = idx + 1;
        end
    end
    imgName = ['./images/', 'pdos', cellType, num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ...
        'Dis', num2str(layerDistanceDatabase(layerIdx)), 'Atom', 'Se', 'SeqNum',num2str(visualPickAtomSe(pickAtomSeIdx)), '.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');
end
close all;
% %% Data visualization for each layer distance (projected DOS for each atom v.s. average DOS)
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
            plot([0 0], pdosVatomVerticalLineLim, '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            grid on;
            legend([fmLineUpAtom fmLineUpAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomV(orbitPickAtomVIdx)), '  FM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: V, Seq.Num ',num2str(visualPickAtomV(pickAtomVIdx))]});
            set(gca,'XMinorGrid','on')
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
            plot([0 0], pdosVatomVerticalLineLim, '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            grid on;
            legend([afmLineUpAtom afmLineUpAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomV(orbitPickAtomVIdx)), ' AFM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: V, Seq.Num ',num2str(visualPickAtomV(pickAtomVIdx))]});
            set(gca,'XMinorGrid','on')
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
            plot([0 0], pdosVatomVerticalLineLim, '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            grid on;
            legend([nmLineAtom nmLineAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomV(orbitPickAtomVIdx)), ' NM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: V, Seq.Num ',num2str(visualPickAtomV(pickAtomVIdx))]});
            set(gca,'XMinorGrid','on')
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
                projectedDosUpFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx)- totVatomNumber), 'Color', 'blue');
            plot(fmdata{layerIdx, 1}.energy, ...
                -projectedDosDownFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx)- totVatomNumber), 'Color', 'blue');
            fmLineUpAverage = plot(projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy, ...
                projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_up'), 'Color', 'red');
            plot(projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy, ...
                -projectedDosFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_down'), 'Color', 'red');
            xlim([-8 8]);
            plot([0 0], pdosVatomSeerticalLineLim, '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            grid on;
            legend([fmLineUpAtom fmLineUpAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), '  FM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: Se'];...
                ['Seq.Num ',num2str(visualPickAtomSe(pickAtomSeIdx))]});
            set(gca,'XMinorGrid','on')
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
                projectedDosUpAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx)- totVatomNumber), 'Color', 'blue');
            plot(afmdata{layerIdx, 1}.energy, ...
                -projectedDosDownAFmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx)- totVatomNumber), 'Color', 'blue');
            afmLineUpAverage = plot(projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy, ...
                projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_up'), 'Color', 'red');
            plot(projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy, ...
                -projectedDosAFmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)+'_down'), 'Color', 'red');
            xlim([-8 8]);
            plot([0 0], pdosVatomSeerticalLineLim, '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            grid on;
            legend([afmLineUpAtom afmLineUpAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ' AFM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: Se'];...
                ['Seq.Num ',num2str(visualPickAtomSe(pickAtomSeIdx))]});
            set(gca,'XMinorGrid','on')
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
                projectedDosNmMatAtomSe{layerIdx, orbitPickAtomSeIdx}(:, visualPickAtomSe(pickAtomSeIdx)- totVatomNumber), 'Color', 'blue');
            nmLineAverage = plot(projectedDosNmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.energy, ...
                projectedDosNmAverageAtomSe{layerIdx, orbitPickAtomSeIdx}.(orbitalPickAtomSe(orbitPickAtomSeIdx)), 'Color', 'red');
            xlim([-8 8]);
            plot([0 0], pdosVatomSeerticalLineLim, '--', 'Color', [0.25 0.25 0.25]);
            hold off;
            grid on;
            legend([nmLineAtom nmLineAverage], {'Atom', 'Average'});
            title({[num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ' NM DOS at D = ', num2str(layerDistanceDatabase(layerIdx))];...
                ['Cell Type: ', cellType, ', Atom: Se'];...
                ['Seq.Num ',num2str(visualPickAtomSe(pickAtomSeIdx))]});
            set(gca,'XMinorGrid','on')
            idx = idx + 1;
        end
    end
    imgName = ['./images/', 'magNm', cellType, num2str(orbitalPickAtomSe(orbitPickAtomSeIdx)), ...
        'Dis', num2str(layerDistanceDatabase(layerIdx)), 'Atom', 'Se', 'SeqNum',num2str(visualPickAtomSe(pickAtomSeIdx)), '.emf'];
    exportgraphics(fig, imgName, 'ContentType', 'vector');
    close all;
end