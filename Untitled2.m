orbit = 3;
c = getWeight(2, 7, PROCAR.totkpoints,orbit,PROCAR);

function weight = getWeight (atomNumber, nbands, totkpoints, orbitIndex, PROCAR)
% orbitIndex:
%     2:  s
%     3:  py
%     4:  pz
%     5:  px
%     6:  dxy
%     7:  dyz
%     8:  dz^2
%     9:  dxz
%     10: dx2_y2
%     11: tot
    weight = zeros (totkpoints, 1);
    for numIdx = 1: totkpoints
        weight(numIdx) = PROCAR.processed{nbands, numIdx}.weight{atomNumber, orbitIndex};        
    end  
end