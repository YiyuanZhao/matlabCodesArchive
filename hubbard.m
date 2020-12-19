clear variables;
load('savedFigureData.mat');
selectedNumIdx = ~isnan(exportData.hubbardModel.bandeigenvalueSelected);
notNanRecipMesh = exportData.hubbardModel.kpointsRecipMeshSelected(1: 3, selectedNumIdx);
Hkin = zeros(1, length(notNanRecipMesh));
for numIdx = 1: length(notNanRecipMesh)
    Hkin(numIdx) = calculateKineticHamiltonian(exportData.hubbardModel.latticeTransMat, notNanRecipMesh(1:3, numIdx), ...
        exportData.hubbardModel.hoppingParameter, exportData.hubbardModel.hoppingMatrixProcessed);
end
exportData.hubbardModel.Hkin = nan(1, length(exportData.hubbardModel.bandeigenvalueSelected));
if prod(imag(Hkin) < 1e-10)
    exportData.hubbardModel.Hkin(selectedNumIdx) = real(Hkin);
end
scatter3(notNanRecipMesh(1, :), notNanRecipMesh(2, :), exportData.hubbardModel.Hkin(selectedNumIdx), 6);

function Hkin = calculateKineticHamiltonian(transMat, kpointsRecip, hoppingParameter, hoppingMatrixProcessed)
            Hkin = 0;
            for hoppingOrder = 1: length(hoppingParameter)
                sum = 0;
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