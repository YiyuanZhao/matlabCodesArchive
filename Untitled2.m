for i = 2: 11
    hoppingMatrixTemp{i} = hoppingMatrixOrigin{i};
    len = length(hoppingMatrixTemp{i}(:, 1));
    hoppingOrder = i - 1;
    hoppingMatrixTemp{i} = [hoppingOrder * ones(len, 1), hoppingMatrixTemp{i}, zeros(len, 1), ones(len, 2)];
    if i == 2
        out = hoppingMatrixTemp{i};
    else
        out = [out ;hoppingMatrixTemp{i}];
    end
%     writecell({[hoppingOrder, len]}, 'hoppingPara.txt', 'WriteMode', 'append', 'Delimiter', 'tab');
%     writemcell(hoppingMatrixTemp(i), 'hoppingPara.txt', 'WriteMode', 'append', 'Delimiter', 'tab');
end