clear variables;
a = 3.32987;
scal = 2;
xOrigin = 1;
yOrigin = 1;
zOrigin = 0;

a1 = [a 0 0];
a2 = [-a/2 sqrt(3)*a/2 0];
groupPosition = [xOrigin - a/4; yOrigin - sqrt(3)*a/4; zOrigin + 0];
% Build the symmetry
SymmetryOperation = cell(1, 24);
SymmetryOperation{1} = [ 1  0  0;  0  1  0;  0  0  1];
SymmetryOperation{2} = [ 0 -1  0;  1 -1  0;  0  0  1];
SymmetryOperation{3} = [-1  1  0; -1  0  0;  0  0  1];
SymmetryOperation{4} = [-1  0  0;  0 -1  0;  0  0  1];
SymmetryOperation{5} = [ 1 -1  0;  1  0  0;  0  0  1];
SymmetryOperation{6} = [ 0  1  0; -1  1  0;  0  0  1];
SymmetryOperation{7} = [ 0 -1  0; -1  0  0;  0  0  1];
SymmetryOperation{8} = [-1  1  0;  0  1  0;  0  0  1];
SymmetryOperation{9} = [ 1  0  0;  1 -1  0;  0  0  1];
SymmetryOperation{10}= [ 0  1  0;  1  0  0;  0  0  1];
SymmetryOperation{11}= [ 1 -1  0;  0 -1  0;  0  0  1];
SymmetryOperation{12}= [-1  0  0; -1  1  0;  0  0  1];
SymmetryOperation{13}= [-1  0  0;  0 -1  0;  0  0 -1];
SymmetryOperation{14}= [ 0  1  0; -1  1  0;  0  0 -1];
SymmetryOperation{15}= [ 1 -1  0;  1  0  0;  0  0 -1];
SymmetryOperation{16}= [ 1  0  0;  0  1  0;  0  0 -1];
SymmetryOperation{17}= [-1  1  0; -1  0  0;  0  0 -1];
SymmetryOperation{18}= [ 0 -1  0;  1 -1  0;  0  0 -1];
SymmetryOperation{19}= [ 0  1  0;  1  0  0;  0  0 -1];
SymmetryOperation{20}= [ 1 -1  0;  0 -1  0;  0  0 -1];
SymmetryOperation{21}= [-1  0  0; -1  1  0;  0  0 -1];
SymmetryOperation{22}= [ 0 -1  0; -1  0  0;  0  0 -1];
SymmetryOperation{23}= [-1  1  0;  0  1  0;  0  0 -1];
SymmetryOperation{24}= [ 1  0  0;  1 -1  0;  0  0 -1];
fig = figure();
hold on;
% Calculate the position
equiv = cell(1, length(SymmetryOperation));
for numIdx = 1: length(SymmetryOperation)
    equiv{numIdx} = SymmetryOperation{numIdx} * groupPosition;
    scatter(equiv{numIdx}(1) + a/4, equiv{numIdx}(2) + sqrt(3)*a/4);
end
latticePos = cell(1, 4);
latticePos{1} = [0 0 0];
latticePos{2} = [a 0 0];
latticePos{3} = [a/2 a*sqrt(3)/2 0];
latticePos{4} = [-a/2 a*sqrt(3)/2 0];

xLineLength = 2* scal * a;
yLineLength = 2* scal * sqrt(3)*a /2;
for numIdx = -scal: scal
    lineMin = -scal * a1 + numIdx * a2;
    plot([lineMin(1), lineMin(1) + xLineLength], [lineMin(2), lineMin(2)], 'r');
end
for numIdx = -scal: scal
    lineMin = numIdx * a1 - scal * a2;
    plot([lineMin(1), lineMin(1) - a/2 * scal * 2], [lineMin(2), lineMin(2) + yLineLength], 'r');
end
hold off;
DisplayMin = [-scal * a1 + scal * a2];
displayRange = [DisplayMin(1), -DisplayMin(1), DisplayMin(1), -DisplayMin(1)];
axis(displayRange);
axis square;