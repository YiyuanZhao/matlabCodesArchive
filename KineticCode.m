function Hkin = KineticCode(transMat, kpointsRecip, hoppingParameter, hoppingMatrixProcessed)
Hkin = complex(0, 0);
for hoppingOrder = 1: length(hoppingParameter)
    sum = complex(0, 0);
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