%% Determine the supercell length
clear variables;
nMax = 40;
mMax = 40;
a = 1;
cutMaxA = 10* a;

theta1 = zeros(nMax, mMax);
theta2 = zeros(nMax, mMax);
theta3 = zeros(nMax, mMax);
theta4 = zeros(nMax, mMax);
length1 = zeros(nMax, mMax);
length2 = zeros(nMax, mMax);
length3 = zeros(nMax, mMax);
length4 = zeros(nMax, mMax);

for numIdxN = 1: nMax
    for numIdxM = 1: mMax
        theta1(numIdxN, numIdxM) = 2*atand((2*numIdxN - 1)/(sqrt(3)* (2*numIdxM - 1)));
        theta2(numIdxN, numIdxM) = 2*atand(numIdxN/(sqrt(3)* numIdxM));
        theta3(numIdxN, numIdxM) = 2*atand((sqrt(3)* (2*numIdxM - 1))/(2*numIdxN - 1));
        theta4(numIdxN, numIdxM) = 2*atand((sqrt(3)* numIdxM)/numIdxN);
        length1(numIdxN, numIdxM) = a/2 * sqrt((2*numIdxN - 1)^2 + 3*(2*numIdxM - 1)^2);
        length2(numIdxN, numIdxM) = a* sqrt(numIdxN^2 + 3* numIdxM^2);
    end
end
length3 = length1;
length4 = length2;
% [sortedTheta1, sortedIndex1] = sort(theta1(:));
% [sortedn1, sortedm1] = ind2sub(size(theta1), sortedIndex1);
% [sortedTheta2, sortedIndex2] = sort(theta2(:));
% [sortedn2, sortedm2] = ind2sub(size(theta2), sortedIndex2);
% 
% figure();
% hold on;
% scatter(sortedTheta1, length1(sortedIndex1), 5, 'b');
% scatter(sortedTheta2, length2(sortedIndex2), 5, 'r');
% hold off;

thetaMix = [theta1, theta2, theta3, theta4];
lengthMix = [length1, length2, length3, length4];
[sortedThetaMix, sortedIndexMix] = sort(thetaMix(:));
[sortednMix, sortedmMix] = ind2sub(size(thetaMix), sortedIndexMix);
while ~prod(sortedmMix <= 20)
    sortedmMix(sortedmMix > 20) = sortedmMix(sortedmMix > 20) - 20;
end
sortedMixType = nan(size(thetaMix));
sortedMixType = 1*(sortedIndexMix(:) <= (numel(thetaMix)/4)) + 2*(sortedIndexMix(:) > (numel(thetaMix)/4) & sortedIndexMix(:) <= (numel(thetaMix)/2)) ...
    + 3* (sortedIndexMix(:) > (numel(thetaMix)/2) & sortedIndexMix(:) <= (numel(thetaMix)*3/4)) + 4*(sortedIndexMix(:) > (numel(thetaMix)*3/4));
% scatter(sortedThetaMix, lengthMix(sortedIndexMix), 5);

[thetaMixUnique, sortedIndexMixUnique] = uniquetol(thetaMix);
lengthMixUnique = nan(size(thetaMixUnique));
sortedmMixUnique = nan(size(thetaMixUnique));
sortednMixUnique = nan(size(thetaMixUnique));
sortedMixTypeUnique = nan(size(thetaMixUnique));
for numIdx = 1: length(thetaMixUnique)
    searchIndex = abs(thetaMixUnique(numIdx) -  sortedThetaMix) < 1e-8;
    lengthMixUnique(numIdx) = min(lengthMix(sortedIndexMix(searchIndex)));
    tmpIndex = searchIndex & abs(lengthMix(sortedIndexMix) - lengthMixUnique(numIdx))< 1e-8;
    sortedmMixUnique(numIdx) = sortedmMix(tmpIndex);
    sortednMixUnique(numIdx) = sortednMix(tmpIndex);
    sortedMixTypeUnique(numIdx) = sortedMixType(tmpIndex);
end

% figure;
% scatter(thetaMixUnique, lengthMixUnique, 5);
% plot(thetaMixUnique, lengthMixUnique);

cutIndex = lengthMixUnique < cutMaxA;
thetaMixUniqueCut = thetaMixUnique(cutIndex);
lengthMixUniqueCut = lengthMixUnique(cutIndex);
sortedmMixUniqueCut = sortedmMixUnique(cutIndex);
sortednMixUniqueCut = sortednMixUnique(cutIndex);
sortedMixTypeUniqueCut = sortedMixTypeUnique(cutIndex);

figure;
scatter(thetaMixUniqueCut, lengthMixUniqueCut, 5);
xlim([0, 60]);

% Clear variables
clear length1 length2 length3 length4 lengthMix theta1 theta2 theta3 theta4 thetaMix
clear cutMaxA mMax nMax numIdx numIdxM numIdxN tmpIndex cutIndex
clear lengthMixUnique searchIndex sortedmMix sortedMixType sortedMixTypeUnique sortedIndexMix sortedmMix sortedmMixUnique
clear sortedIndexMixUnique sortednMix sortednMixUnique sortedThetaMix thetaMixUnique
%% Generates the lattice
a = 1;  % a-axies of the lattice
b = a;      % b-axies of the lattice
gamma = 120;% angle of <a, b>
sizeLattice = 21;  % Scale of the system (should be odd)

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

%% Rotation
% Construct layers
d = 0.05;
D = 0.5;
% rotationNumIdx = 20;
for rotationNumIdx = 19: 19
rotationDegree = thetaMixUniqueCut(rotationNumIdx)/2;
supercellLattice = lengthMixUniqueCut(rotationNumIdx);
supercellBoundx = [0, supercellLattice, supercellLattice/2, -supercellLattice/2, 0];
supercellBoundy = [0, 0, sqrt(3)/2 * supercellLattice, sqrt(3)/2 * supercellLattice, 0];
offsetSe1 = [deltaA', deltaB']*[2/3; 1/3];
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
layer{2}.Se1.z = layer{1}.Se1.z + D;
layer{2}.Se2.z = layer{1}.Se2.z + D;
layer{2}.V.z = layer{1}.V.z + D;

% Rotation Process
for i = 1: 2
    % layer0: counterclock; layer1: clock
    switch i
        case 1
            rotationMatrix = ...
                [cosd(rotationDegree), -sind(rotationDegree), 0; sind(rotationDegree), cosd(rotationDegree), 0; 0 0 1];
        case 2
            rotationMatrix = ...
                [cosd(-rotationDegree), -sind(-rotationDegree), 0; sind(-rotationDegree), cosd(-rotationDegree), 0; 0 0 1];
    end
    switch sortedMixTypeUniqueCut(rotationNumIdx)
        case {1, 2}
            correctMatrix = [cosd(-30), -sind(-30), 0; sind(-30), cosd(-30), 0; 0 0 1];
        case {3, 4}
            correctMatrix = [1 0 0; 0 1 0; 0 0 1];
    end
    for numIdx = 1: numel(layer{i}.V.x)
        tmpMat = [layer{i}.Se1.x(numIdx), layer{i}.Se1.y(numIdx), layer{i}.Se1.z(numIdx)]';
        tmpMat = correctMatrix * rotationMatrix * tmpMat;
        layer{i}.Se1.x(numIdx) = tmpMat(1);
        layer{i}.Se1.y(numIdx) = tmpMat(2);
        layer{i}.Se1.z(numIdx) = tmpMat(3);
        
        tmpMat = [layer{i}.Se2.x(numIdx), layer{i}.Se2.y(numIdx), layer{i}.Se2.z(numIdx)]';
        tmpMat = correctMatrix * rotationMatrix * tmpMat;
        layer{i}.Se2.x(numIdx) = tmpMat(1);
        layer{i}.Se2.y(numIdx) = tmpMat(2);
        layer{i}.Se2.z(numIdx) = tmpMat(3);
        
        tmpMat = [layer{i}.V.x(numIdx), layer{i}.V.y(numIdx), layer{i}.V.z(numIdx)]';
        tmpMat = correctMatrix * rotationMatrix * tmpMat;
        layer{i}.V.x(numIdx) = tmpMat(1);
        layer{i}.V.y(numIdx) = tmpMat(2);
        layer{i}.V.z(numIdx) = tmpMat(3);
    end
    % Cut the supercell
    layer{i}.V.inSupercell = inpolygon(layer{i}.V.x, layer{i}.V.y, supercellBoundx, supercellBoundy);
end

% 
% fig = figure;
% hold on;
% % scatter3(layer{1}.Se1.x(:), layer{1}.Se1.y(:), layer{1}.Se1.z(:), 'green');
% scatter3(layer{1}.V.x(:), layer{1}.V.y(:), layer{1}.V.z(:), 'red');
% % scatter3(layer{1}.Se2.x(:), layer{1}.Se2.y(:), layer{1}.Se2.z(:), 'green');
% % scatter3(layer{2}.Se1.x(:), layer{2}.Se1.y(:), layer{2}.Se1.z(:), 'green');
% scatter3(layer{2}.V.x(:), layer{2}.V.y(:), layer{2}.V.z(:), 'blue');
% % scatter3(layer{2}.Se2.x(:), layer{2}.Se2.y(:), layer{2}.Se2.z(:), 'green');
% hold off
% maxLim = max(abs([fig.Children.XLim, fig.Children.YLim]));
% fig.Children.XLim = [-maxLim, maxLim];
% fig.Children.YLim = [-maxLim, maxLim];

fig = figure;
hold on;
% scatter3(layer{1}.Se1.x(:), layer{1}.Se1.y(:), layer{1}.Se1.z(:), 'green');
scatter3(layer{1}.V.x(layer{i}.V.inSupercell), layer{1}.V.y(layer{i}.V.inSupercell), layer{1}.V.z(layer{i}.V.inSupercell), 'red');
% scatter3(layer{1}.Se2.x(:), layer{1}.Se2.y(:), layer{1}.Se2.z(:), 'green');
% scatter3(layer{2}.Se1.x(:), layer{2}.Se1.y(:), layer{2}.Se1.z(:), 'green');
scatter3(layer{2}.V.x(layer{i}.V.inSupercell), layer{2}.V.y(layer{i}.V.inSupercell), layer{2}.V.z(layer{i}.V.inSupercell), 'blue');
% scatter3(layer{2}.Se2.x(:), layer{2}.Se2.y(:), layer{2}.Se2.z(:), 'green');
hold off
maxLim = max(abs([fig.Children.XLim, fig.Children.YLim]));
fig.Children.XLim = [-maxLim, maxLim];
fig.Children.YLim = [-maxLim, maxLim];

end
