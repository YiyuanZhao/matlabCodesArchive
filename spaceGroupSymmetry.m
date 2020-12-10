clear variables;
a = 3.32987;
xOrigin = 0;
yOrigin = 0;
zOrigin = 0;

groupPosition = [a/4 + xOrigin; sqrt(3)*a/4 + yOrigin; 0 + zOrigin];
% Build the symmetry
SymmetryOperation = cell(1, 24);
SymmetryOperation{1} = [ 1  0 0;  0  1 0; 0 0 1];
SymmetryOperation{2} = [ 0 -1 0;  1 -1 0; 0 0 1];
SymmetryOperation{3} = [-1  1 0; -1  0 0; 0 0 1];
SymmetryOperation{4} = [-1  0 0;  0 -1 0; 0 0 1];
SymmetryOperation{5} = [1  -1 0;  1  0 0; 0 0 1];
SymmetryOperation{6} = [0   1 0; ]
% Calculate the position
for numIdx = 1: 1
    equiv{numIdx} = SymmetryOperation{numIdx} * groupPosition;
end