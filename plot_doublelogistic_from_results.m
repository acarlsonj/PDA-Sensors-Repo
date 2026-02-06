function plot_doublelogistic_from_results(resultsCsv, dataDir)


close all; clc;
applyPublicationDefaults();

if nargin < 2 || isempty(dataDir)
    dataDir = pwd;
end
if exist(dataDir,'dir') ~= 7
    error('Data directory not found: %s', dataDir);
end


if nargin < 1 || isempty(resultsCsv)
    resultsCsv = 'ManualGlobalSearch_DoubleLogistic_FitResults.csv';
end

if exist(resultsCsv,'file') ~= 2
    alt = fullfile(dataDir, resultsCsv);
    if exist(alt,'file') ~= 2
        error('Results CSV not found: %s', resultsCsv);
    end
    resultsCsv = alt;
end

T = readtable(resultsCsv);

required = {'Conc_mg_mL','Time_s','a','b','c','e','f','x01','x02'};
for k = 1:numel(required)
    if ~ismember(required{k}, T.Properties.VariableNames)
        error('Results CSV is missing required column: %s', required{k});
    end
end


T_hotplate = [90:2.5:140, 150:10:250];
T_actual   = [50:2.5:100, 110:10:210];


concList = sort(unique(T.Conc_mg_mL(:)));
timeList = sort(unique(T.Time_s(:)));


concColors = lines(numel(concList));


expMarker  = 'o';
mSize      = 5;


xLimFixed = [50, 210];


padFracY = 0.03;


for ti = 1:numel(timeList)
    tsec = timeList(ti);
    showLegend = (ti == 1);

    figure('Name', sprintf('%ds: Experimental (dots) vs Model (dashed)', tsec));
    hold on;

   showLegend = (ti == 1);


if showLegend
    nC = numel(concList);
    legH = gobjects(2*nC, 1);
    legL = cell(2*nC, 1);

    li = 1;


    for ci = 1:nC
        conc = concList(ci);
        legH(li) = plot(nan, nan, '--', ...
            'Color', concColors(ci,:));
        legL{li} = sprintf('%d mg/mL Modeled', conc);
        li = li + 1;
    end


    for ci = 1:nC
        conc = concList(ci);
        legH(li) = plot(nan, nan, ...
            'LineStyle','none', ...
            'Marker', expMarker, ...
            'MarkerSize', mSize, ...
            'MarkerFaceColor', concColors(ci,:), ...
            'MarkerEdgeColor', concColors(ci,:));
        legL{li} = sprintf('%d mg/mL Experimental', conc);
        li = li + 1;
    end
end


    yAll = [];

    for ci = 1:numel(concList)
        conc = concList(ci);


        rowMask = (T.Conc_mg_mL == conc) & (T.Time_s == tsec);
        if ~any(rowMask)
            continue;
        end
        r = T(find(rowMask,1,'first'),:);


        p = [r.a, r.b, r.c, r.x01, r.e, r.f, r.x02];


        fname = fullfile(dataDir, sprintf('%dmg_%ds_Background_Adjusted.xlsx', conc, tsec));S = readLabFile(fname);

        rawTemps = S.raw(:,1);
        x = correctTempsFromTable(rawTemps, T_hotplate, T_actual);
        y = S.L(:);

        [x,y] = sanitizeXY(x,y);


        inWin = (x >= xLimFixed(1)) & (x <= xLimFixed(2));
        xw = x(inWin);
        yw = y(inWin);


        if numel(xw) < 2
            continue;
        end


        plot(xw, yw, ...
    'LineStyle','none', ...
    'Marker', expMarker, ...
    'MarkerSize', mSize, ...
    'MarkerFaceColor', concColors(ci,:), ...
    'MarkerEdgeColor', concColors(ci,:), ...
    'HandleVisibility','off');


        xFit = linspace(xLimFixed(1), xLimFixed(2), 800).';
        yFit = modelOnly(p, xFit);

        plot(xFit, yFit, '--', ...
    'Color', concColors(ci,:), ...
    'HandleVisibility','off');


        yAll = [yAll; yw(:); yFit(:)];
    end

    xlabel('Temperature ($^{\circ}$C)');
    ylabel('Avg $L^*$');
    grid minor;


    xlim(xLimFixed);


    if ~isempty(yAll)
        yMin = min(yAll);
        yMax = max(yAll);
        dy   = max(yMax - yMin, eps);
        ylim([yMin - padFracY*dy, yMax + padFracY*dy]);
    end


    tightenAxesWithPadding(gca, 0.02);


    if showLegend
    legend(legH, legL, 'Location','southeast');
    end

end

end


function yfit = modelOnly(p, x)
a = p(1); b = p(2); c = p(3); x01 = p(4);
e = p(5); f = p(6); x02 = p(7);

[S1, ~] = safeSigmoidAndDeriv(c.*(x - x01));
[S2, ~] = safeSigmoidAndDeriv(f.*(x - x02));

yfit = b + (a-b).*S1 + (e-a).*S2;
end

function [S, Sp] = safeSigmoidAndDeriv(z)
z = max(min(z, 60), -60);
S  = 1 ./ (1 + exp(-z));
Sp = S .* (1 - S);
end


function S = readLabFile(fname)
if exist(fname,'file') ~= 2
    error('Experimental data file not found: %s', fname);
end

if exist('readmatrix','file') == 2
    raw = readmatrix(fname);
else
    tmp = importdata(fname);
    if isstruct(tmp) && isfield(tmp,'data')
        raw = tmp.data;
    else
        raw = tmp;
    end
end

if isempty(raw) || size(raw,2) < 4
    error('File %s does not have >= 4 numeric columns.', fname);
end

S = struct();
S.raw = raw;
S.L   = raw(:,2);
S.a   = raw(:,3);
S.b   = raw(:,4);
end

function [x,y] = sanitizeXY(x,y)
x = x(:); y = y(:);
n = min(numel(x), numel(y));
x = x(1:n); y = y(1:n);

mask = isfinite(x) & isfinite(y);
x = x(mask); y = y(mask);

[x, idx] = sort(x);
y = y(idx);

[x, ia] = unique(x,'stable');
y = y(ia);

if numel(x) < 8
    error('Not enough valid points after cleaning (n=%d).', numel(x));
end
end

function temps_corr = correctTempsFromTable(temps_in, T_hotplate, T_actual)
temps_in   = temps_in(:);
T_hotplate = T_hotplate(:);
T_actual   = T_actual(:);

[T_hotplate, idx] = sort(T_hotplate);
T_actual = T_actual(idx);

temps_corr = interp1(T_hotplate, T_actual, temps_in, 'linear', 'extrap');
end


function applyPublicationDefaults()
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultFigurePaperPositionMode', 'auto');
set(groot, 'defaultFigureRenderer', 'painters');


pubFont = 'Helvetica';
baseFS  = 18;
lw      = 2;
ms      = 20;

set(groot, ...
    'defaultFigureColor', 'w', ...
    'defaultFigureUnits', 'inches', ...
    'defaultFigurePosition', [1 1 6 4], ...
    'defaultAxesFontName', pubFont, ...
    'defaultAxesFontSize', baseFS, ...
    'defaultAxesLineWidth', 1, ...
    'defaultAxesBox', 'on', ...
    'defaultAxesTickDir', 'out', ...
    'defaultAxesTickLength', [0.015 0.015], ...
    'defaultAxesXGrid', 'on', ...
    'defaultAxesYGrid', 'on', ...
    'defaultAxesLabelFontSizeMultiplier', 1.1, ...
    'defaultAxesTitleFontSizeMultiplier', 1.1, ...
    'defaultLineLineWidth', lw, ...
    'defaultLineMarkerSize', ms, ...
    'defaultLegendBox', 'on', ...
    'defaultLegendFontSize', baseFS);
end

function tightenAxesWithPadding(ax, pad)

if nargin < 2 || isempty(pad), pad = 0.02; end

u0 = ax.Units;
ax.Units = 'normalized';

ti = ax.TightInset;
left   = ti(1) + pad;
bottom = ti(2) + pad;
width  = 1 - (ti(1) + ti(3)) - 2*pad;
height = 1 - (ti(2) + ti(4)) - 2*pad;


width  = max(width, 0.1);
height = max(height, 0.1);

ax.Position = [left, bottom, width, height];
ax.Units = u0;
end
