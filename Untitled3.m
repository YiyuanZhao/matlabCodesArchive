%% ���ı��ļ��е�������
% ���ڴ������ı��ļ��е������ݵĽű�:
%
%    filename: C:\temp\EIGENVAL2
%
% �� MATLAB �� 2020-09-14 11:39:32 �Զ�����

%% ���õ���ѡ��
opts = delimitedTextImportOptions("NumVariables", 25);

% ָ����Χ�ͷָ���
opts.DataLines = [8, Inf];
opts.Delimiter = " ";

% ָ�������ƺ�����
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5"];
opts.SelectedVariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];
opts = setvaropts(opts, [1, 2], "TrimNonNumeric", true);
opts = setvaropts(opts, [1, 2], "ThousandsSeparator", ",");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% ��������
EIGENVAL2 = readtable("C:\temp\EIGENVAL3", opts);

%% ת��Ϊ�������
EIGENVAL2 = table2cell(EIGENVAL2);
numIdx = cellfun(@(x) ~isnan(str2double(x)), EIGENVAL2);
EIGENVAL2(numIdx) = cellfun(@(x) {str2double(x)}, EIGENVAL2(numIdx));
numIdx = cellfun(@(x) isnan((x)), EIGENVAL2(:, 1));
EIGENVAL2(numIdx, :) = [];
numIdx = cellfun(@(x) mod((x),1), EIGENVAL2(:, 1));
i = 2;
while ~numIdx(i)
    i = i + 1;
end
 totBands = i - 2;
 totKpoints = length(EIGENVAL2)/(totBands + 1);
%% �����ʱ����
clear opts numIdx i

%% Test Parts
readBaseLine = 1;
kpointsWeight = zeros (totKpoints, 1);
kpointsMesh = zeros (3, totKpoints);
bandEigenvalue = cell(totBands, totKpoints);
bandOccupation = cell(totBands, totKpoints);

for orderKpoint = 1: totKpoints
    kpointsWeight(orderKpoint) = EIGENVAL2{readBaseLine, 4};
    kpointsMesh(1, orderKpoint) = EIGENVAL2{readBaseLine, 1};
    kpointsMesh(2, orderKpoint) = EIGENVAL2{readBaseLine, 2};
    kpointsMesh(3, orderKpoint) = EIGENVAL2{readBaseLine, 3};
    readBaseLine = readBaseLine + 1;
    for orderBands = 1: totBands
        bandEigenvalue{orderBands, orderKpoint} = EIGENVAL2{readBaseLine + orderBands - 1, 2};
        bandOccupation{orderBands, orderKpoint} = EIGENVAL2{readBaseLine + orderBands - 1, 3};
    end
    readBaseLine = readBaseLine + totBands;
end
numIdx = kpointsMesh(3, :) < 1e-15;
kpointsMesh(3, numIdx) = 0;
kpointZvalue = unique(kpointsMesh(3, :));