%% Normal Single/Multiple Lines Plot
clear variables;
close all;
% Variables Setup
data = importdata('plotdata.csv', ',');
% Parameter Setting
PlotType = 'plot';  % 'plot', 'area';
% xlab = "Lattice Parameter";
ylab = "$\rho_\mathrm{atom} - \bar{\rho}$ $(e^{-})$";
dataDir = 'D:\OneDrive - tongji.edu.cn\Research\thesis\figureData';
figDir = 'D:\OneDrive - tongji.edu.cn\Research\thesis\TongjiThesis-master\fig';
fileName = 'CDWNew';
XMinorGrid = 0;
YMinorGrid = 0;
printMarker = 0;
treatAsStr = 0;
Interpreter = 'latex';
TickLabelInterpreter = 'latex';
XTickLabelRotation = 0;
legendLocation = 'northwest'; % https://ww2.mathworks.cn/help/matlab/ref/legend.html?s_tid=doc_ta#bt6ef_q-1-lcn
% xlimit = [0 11];
ylimit = [-0.1 0.1];

% Plot figure
fig = figure;
ax = axes;
if isfield(data, 'textdata') && length(data.textdata(:, 1)) > 1
    x = 1: length(data.textdata(2: end, 1));
    y = data.data;
elseif isfield(data, 'textdata')
    if treatAsStr == 0
        x = data.data(:, 1);
    else
        x = 1:length(data.data(:, 1));
        XTickLabel = data.data(:, 1);
    end
    y = data.data(:, 2: end);
else
    if treatAsStr == 0
        x = data(:, 1);
    else
         x = 1: length(data.data(:, 1));
         XTickLabel = data(:, 1);
    end
    y = data(:, 2: end);
end
hold on;
switch PlotType
    case 'plot'
        for numIdx = 1: length(y(1, :))
            if printMarker == 1
                pt = plot(x, y(:, numIdx), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6);
            else
                pt = plot(x, y(:, numIdx), 'LineWidth', 1);
            end
        end
    case 'area'
        for numIdx = 1: length(y(1, :))
            pt = area(x, y(:, numIdx), 'LineWidth', 1);
            pt.FaceAlpha = 0.5;
        end
    otherwise
        error("Variables 'PlotType' should be 'plot' or 'area'.");
end
hold off;
% Enable Grid
grid on;
% XMinorGrid Setting
if XMinorGrid == 1
    ax.XMinorGrid = 'on';
end
% YMinorGrid Setting
if YMinorGrid == 1
    ax.YMinorGrid = 'on';
end
% Label Setting
if exist('xlab', 'var')
    xlabel(xlab, 'FontSize', 14);
else
    xlabel(data.textdata(1, 1), 'FontSize', 14);
end
ylabel(ylab, 'FontSize', 14);
ax.XLabel.Interpreter = Interpreter;
ax.YLabel.Interpreter = Interpreter;
% Legend Setting
legend(data.textdata(1, 2: end), 'FontSize', 14, 'Location', legendLocation);
ax.Legend.Interpreter = Interpreter;
% Limits Setting (If exist)
if exist('xlimit', 'var')
    xlim(xlimit)
end
if exist('ylimit', 'var')
    ylim(ylimit)
end
% XTickLabel Setting (If exist)
if isfield(data, 'textdata') && length(data.textdata(:, 1)) > 1
    ax.XTickLabel = data.textdata(2: end, 1);
elseif treatAsStr == 1
    ax.XTick = x;
    ax.XTickLabel = XTickLabel;
end
% Add Box of Figure
ax.Box = 'on';
ax.LineWidth = 1;

% Tick Setting
ax.TickLabelInterpreter = TickLabelInterpreter;
ax.XTickLabelRotation = XTickLabelRotation;

% Save Figure
if exist('fileName', 'var')
    saveas(fig, [dataDir, '\', fileName, '.fig']);
else
    saveas(fig, [dataDir, '\untitled.fig']);
end
% Export Figure (.eps format)
if exist('fileName', 'var')
    exportgraphics(fig, [figDir, '\', fileName, '.eps'], 'ContentType', 'vector');
else
    exportgraphics(fig, [figDir, '\untitled.eps'], 'ContentType', 'vector');
end


%% Normal Stack/Lines Plot (Column 1, 2 are stacking, column 3 and after are lines)
clear variables;
close all;
% Variables Setup
data = importdata('plotdata.csv', ',');
rightColorMap = [...
    0.6350 0.0780 0.1840; ...
    0.3010 0.7450 0.9330; ...
    0.4660 0.6740 0.1880; ...
    0.4940 0.1840 0.5560; ...
    0.9290 0.6940 0.1250; ...
    0.8500 0.3250 0.0980; ...
    0 0.4470 0.7410; ...
];
% Parameter Setting
% xlab = "Lattice Parameter";
ylabLeft = "E$_\mathrm{tot}$ (eV)";
ylabRight = "$\Delta$ E (meV)";
dataDir = 'D:\OneDrive - tongji.edu.cn\Research\thesis\figureData';
figDir = 'D:\OneDrive - tongji.edu.cn\Research\thesis\TongjiThesis-master\fig';
fileName = 'phaseTransitionOrder';
XMinorGrid = 0;
YMinorGrid = 0;
printMarker = 1;
printLabel = 0;
treatAsStr = 1;
Interpreter = 'latex';
TickLabelInterpreter = 'latex';
XTickLabelRotation = 0;
legendLocation = 'northwest'; % https://ww2.mathworks.cn/help/matlab/ref/legend.html?s_tid=doc_ta#bt6ef_q-1-lcn
% xlimit = [0 11];
ylimitLeft = [-5.95 -5.65];
ylimitRight = [-1.2 1.2];

% Plot figure
fig = figure;
ax = axes;
if isfield(data, 'textdata') && length(data.textdata(:, 1)) > 1
    x = 1: length(data.textdata(2: end, 1));
    y = data.data;
elseif isfield(data, 'textdata')
    if treatAsStr == 0
        x = data.data(:, 1);
    else
        x = 1:length(data.data(:, 1));
        XTickLabel = data.data(:, 1);
    end
    y = data.data(:, 2: end);
else
    if treatAsStr == 0
        x = data(:, 1);
    else
         x = 1: length(data.data(:, 1));
         XTickLabel = data(:, 1);
    end
    y = data(:, 2: end);
end
hold on;
% Bar Plotting
yyaxis left;
colororder('default');
if printLabel == 1
    pt = bar(x, y(:, 1:2));
    xtips1 = pt(1).XEndPoints;
    ytips1 = pt(1).YEndPoints;
    labels1 = string(pt(1).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','top');
    xtips2 = pt(2).XEndPoints;
    ytips2 = pt(2).YEndPoints;
    labels2 = string(pt(2).YData);
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','top');
else
    pt = bar(x, y(:, 1:2));
end
pt(1).FaceAlpha = 0.8;
pt(2).FaceAlpha = 0.8;
% Line Plotting
yyaxis right;
colororder(rightColorMap);
if length(y(1, :)) > 2
    for numIdx = 3: length(y(1, :))
        if printMarker == 1
            pt = plot(x, y(:, numIdx), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6);
        else
            pt = plot(x, y(:, numIdx), 'LineWidth', 1);
        end
        ax.YColor = ax.ColorOrder(ax.ColorOrderIndex - 1, :);
        ax.YLabel.Color = ax.ColorOrder(ax.ColorOrderIndex - 1, :);
    end
end
hold off;
% Enable Grid
grid on;
% XMinorGrid Setting
if XMinorGrid == 1
    ax.XMinorGrid = 'on';
end
% YMinorGrid Setting
if YMinorGrid == 1
    ax.YMinorGrid = 'on';
end
% Label Setting
if exist('xlab', 'var')
    xlabel(xlab, 'FontSize', 14);
else
    xlabel(data.textdata(1, 1), 'FontSize', 14);
end
ax.XLabel.Interpreter = Interpreter;
yyaxis left
ylabel(ylabLeft, 'FontSize', 14);
ax.YLabel.Interpreter = Interpreter;
yyaxis right
ylabel(ylabRight, 'FontSize', 14);
ax.YLabel.Interpreter = Interpreter;
% Legend Setting
legend(data.textdata(1, 2: end), 'FontSize', 14, 'Location', legendLocation);
ax.Legend.Interpreter = Interpreter;
% Limits Setting (If exist)
if exist('xlimit', 'var')
    xlim(xlimit)
end
if exist('ylimitLeft', 'var')
    yyaxis left
    ylim(ylimitLeft)
end
if exist('ylimitRight', 'var')
    yyaxis right
    ylim(ylimitRight)
end
% XTickLabel Setting (If exist)
if isfield(data, 'textdata') && length(data.textdata(:, 1)) > 1
    ax.XTickLabel = data.textdata(2: end, 1);
elseif treatAsStr == 1
    ax.XTick = x;
    ax.XTickLabel = XTickLabel;
end
% Add Box of Figure
ax.Box = 'on';
ax.LineWidth = 1;

% Tick Setting
ax.TickLabelInterpreter = TickLabelInterpreter;
ax.XTickLabelRotation = XTickLabelRotation;

% Save Figure
if exist('fileName', 'var')
    saveas(fig, [dataDir, '\', fileName, '.fig']);
else
    saveas(fig, [dataDir, '\untitled.fig']);
end

% Export Figure (.eps format)
if exist('fileName', 'var')
    exportgraphics(fig, [figDir, '\', fileName, '.eps'], 'ContentType', 'vector');
else
    exportgraphics(fig, [figDir, '\untitled.eps'], 'ContentType', 'vector');
end
%% Edit Figure from VASP Plotter
clear variables;
fig = gcf;
ax = gca;
fig.Position = [218         522        1077         420];
ax.Units = 'normalized';
ax.InnerPosition = [0.1300    0.1100    0.7750    0.8150];
ax.OuterPosition = [0     0     1     1];

% xlab = "Lattice Parameter";
ylab = "Energy (eV)";
xlab = "K-Path";
dataDir = 'D:\OneDrive - tongji.edu.cn\Research\thesis\figureData';
figDir = 'D:\OneDrive - tongji.edu.cn\Research\thesis\TongjiThesis-master\fig';
fileName = 'wannierBulkBandNum11';
% legendStr = {'$d_{xy}$', '$d_{yz}$', 'd_{xz}'};
XMinorGrid = 0;
YMinorGrid = 1;
Interpreter = 'latex';
legendLocation = 'northwest'; % https://ww2.mathworks.cn/help/matlab/ref/legend.html?s_tid=doc_ta#bt6ef_q-1-lcn
% xlimit = [0 60];
% ylimit = [-0.1 0.1];

ax.LineWidth = 1;

grid on;
% XMinorGrid Setting
if XMinorGrid == 1
    ax.XMinorGrid = 'on';
end
% YMinorGrid Setting
if YMinorGrid == 1
    ax.YMinorGrid = 'on';
end
% Label Setting
if exist('xlab', 'var')
    xlabel(xlab, 'FontSize', 14);
end
ylabel(ylab, 'FontSize', 14);
ax.XLabel.Interpreter = Interpreter;
ax.YLabel.Interpreter = Interpreter;
% Legend Setting
if exist('legendStr', 'var')
    legend(legendStr, 'FontSize', 10, 'Location', legendLocation);
    ax.Legend.Interpreter = Interpreter;
else
    legend('off');
end
% Limits Setting (If exist)
if exist('xlimit', 'var')
    xlim(xlimit)
end
if exist('ylimit', 'var')
    ylim(ylimit)
end
% XTickLabel Setting (If exist)
% if isfield(data, 'textdata') && length(data.textdata(:, 1)) > 1
%     ax.XTickLabel = data.textdata(2: end, 1);
% elseif treatAsStr == 1
%     ax.XTick = x;
%     ax.XTickLabel = XTickLabel;
% end
% Add Box of Figure
ax.Box = 'on';
ax.LineWidth = 1;

% Tick Setting
% ax.TickLabelInterpreter = TickLabelInterpreter;
% ax.XTickLabelRotation = XTickLabelRotation;
ax.FontSize = 22;

% For Fat Band
% ax.Legend.Interpreter = 'latex';
% legend([ax.Children(653), ax.Children(655), ax.Children(656)], {'$d_{xz}$', '$d_{yz}$', '$d_{xy}$'},'NumColumns',3, 'Position', [0.6996    0.8456    0.2182    0.0687]);
% legend([ax.Children(652), ax.Children(654)], {'$d_{x^2-y^2}$', '$d_{z^2-r^2}$'},'NumColumns',2, 'Position', [0.7283    0.8480    0.1887    0.0687]);
% legend([ax.Children(657), ax.Children(658), ax.Children(659)], {'$p_{x}$', '$p_{z}$', '$p_{y}$'},'NumColumns',3, 'Position', [0.6996    0.8456    0.2182    0.0687]);
% ax.Legend.Interpreter = 'latex';


% Save Figure
if exist('fileName', 'var')
    saveas(fig, [dataDir, '\', fileName, '.fig']);
else
    saveas(fig, [dataDir, '\untitled.fig']);
end
% Export Figure (.eps format)
if exist('fileName', 'var')
    exportgraphics(fig, [figDir, '\', fileName, '.png'], 'ContentType', 'image');
else
    exportgraphics(fig, [figDir, '\untitled.eps'], 'ContentType', 'vector');
end
%% Area Plot
clear variables;
fig = gcf;
ax = gca;
fig.Position = [218         522        1077         420];
ax.Units = 'normalized';
ax.InnerPosition = [0.1300    0.1100    0.7750    0.8150];
ax.OuterPosition = [0     0     1     1];
dotLineX = ax.Children(1).XData;
dotLineY = ax.Children(1).YData;
dosLineX = ax.Children(2).XData;
dosLineY = ax.Children(2).YData;
% xlab = "Lattice Parameter";
ylab = "DOS";
xlab = "Energy (eV)";
dataDir = 'D:\OneDrive - tongji.edu.cn\Research\thesis\figureData';
figDir = 'D:\OneDrive - tongji.edu.cn\Research\thesis\TongjiThesis-master\fig';
fileName = 'wannierBulkBandNum11dos';
legendStr = 'DFT';
XMinorGrid = 0;
YMinorGrid = 0;
printMarker = 0;
treatAsStr = 0;
Interpreter = 'latex';
TickLabelInterpreter = 'latex';
XTickLabelRotation = 0;
legendLocation = 'northwest'; % https://ww2.mathworks.cn/help/matlab/ref/legend.html?s_tid=doc_ta#bt6ef_q-1-lcn
xlimit = [-8 8];
ylimit = [0 13];

fig = figure;
fig.Position = [180   435   939   420];
ax = axes;
ax.XMinorGrid = 'on';
ax.XMinorTick = 'on';
hold on;
ar = area(dosLineX, dosLineY);
ar.FaceAlpha = 0.2;
ar.EdgeColor = ar.FaceColor;
ar.LineWidth = 1;
plt = plot(dotLineX, dotLineY, '--k');
plt.LineWidth = 1;
hold off;

grid on;
% XMinorGrid Setting
if XMinorGrid == 1
    ax.XMinorGrid = 'on';
end
% YMinorGrid Setting
if YMinorGrid == 1
    ax.YMinorGrid = 'on';
end
% Label Setting
if exist('xlab', 'var')
    xlabel(xlab, 'FontSize', 14);
end
ylabel(ylab, 'FontSize', 14);
ax.XLabel.Interpreter = Interpreter;
ax.YLabel.Interpreter = Interpreter;
% Legend Setting
if exist('legendStr', 'var')
    legend(legendStr, 'FontSize', 18, 'Location', legendLocation);
    ax.Legend.Interpreter = Interpreter;
    ax.Legend.Box = 'off';
else
    legend('off');
end
% Limits Setting (If exist)
if exist('xlimit', 'var')
    xlim(xlimit)
end
if exist('ylimit', 'var')
    ylim(ylimit)
end
% XTickLabel Setting (If exist)
% if isfield(data, 'textdata') && length(data.textdata(:, 1)) > 1
%     ax.XTickLabel = data.textdata(2: end, 1);
% elseif treatAsStr == 1
%     ax.XTick = x;
%     ax.XTickLabel = XTickLabel;
% end
% Add Box of Figure
ax.Box = 'on';
ax.LineWidth = 1;

% Tick Setting
% ax.TickLabelInterpreter = TickLabelInterpreter;
% ax.XTickLabelRotation = XTickLabelRotation;
ax.FontSize = 18;
% Save Figure
if exist('fileName', 'var')
    saveas(fig, [dataDir, '\', fileName, '.fig']);
else
    saveas(fig, [dataDir, '\untitled.fig']);
end
% Export Figure (.eps format)
if exist('fileName', 'var')
    exportgraphics(fig, [figDir, '\', fileName, '.eps'], 'ContentType', 'vector');
else
    exportgraphics(fig, [figDir, '\untitled.eps'], 'ContentType', 'vector');
end