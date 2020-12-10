queryHopping = hoppingMatrixOrderLoop;
sourceHoppingMatrix = hoppingMatrixOrigin;
sourceHoppingParameter = hoppingParameterOrderLoop;
[mat, para] = queryHoppingParameter(queryHopping, sourceHoppingMatrix, sourceHoppingParameter);

function [sortedHoppingMatrix, sortedHoppingParameter] = queryHoppingParameter(queryHopping, sourceHoppingMatrix, sourceHoppingParameter)
% sortedHoppingMatrix = cell(1, length(queryHopping));
% sortedHoppingParameter = cell(1, length(queryHopping));
numIdxLoopTerminationFlag = 0;
for numIdxLoop1 = 2: length(queryHopping)
    for numIdxLoop2 = 1: length(sourceHoppingMatrix)
        for numIdxLoop3 = 1: length(sourceHoppingMatrix{numIdxLoop2}(:, 1))
            if numIdxLoopTerminationFlag == 1
                numIdxLoopTerminationFlag = 0;
                break;
            end
            if prod(queryHopping{numIdxLoop1}(1, 1: 2) == sourceHoppingMatrix{numIdxLoop2}(numIdxLoop3, :))
                sortedHoppingMatrix{numIdxLoop2} = sourceHoppingMatrix{numIdxLoop2};
%                 sortedHoppingParameter{numIdxLoop2} = sourceHoppingParameter{numIdxLoop2};
                numIdxLoopTerminationFlag = 1;
            end
        end
    end
end
end
