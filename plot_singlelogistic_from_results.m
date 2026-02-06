function plot_singlelogistic_from_results(resultsCsv, dataDir)
if nargin < 1 || isempty(resultsCsv)
    resultsCsv = 'ManualGlobalSearch_SingleLogistic_FitResults.csv';
end
if nargin < 2 || isempty(dataDir)
    dataDir = pwd;
end

close all;
clc;

applyPublicationDefaults();

if exist(resultsCsv,'file') ~= 2
    alt = fullfile(dataDir, resultsCsv);
    if exist(alt,'file') == 2
        resultsCsv = alt;
    else
        error('Results CSV not found: %s', resultsCsv);
    end
end

T = readtable(resultsCsv);

required = {'Conc_mg_mL','Time_s','a','b','c','x0'};
for k = 1:numel(required)
    if ~ismember(required{k}, T.Properties.VariableNames)
        error('Results CSV is missing required column: %s', required{k});
    end
end

T_hotplate = [90:2.5:140, 150:10:250];
T_actual   = [50:2.5:100, 110:10:210];

concList = unique(T.Conc_mg_mL(:));
timeList = unique(T.Time_s(:));
concList = sort(concList(:)).';
timeList = sort(timeList(:)).';

C = lines(numel(concList));
ms = 18;

for ti = 1:numel(timeList)
    tsec = timeList(ti);

    figure('Name', sprintf('Single Logistic Fits - %d s', tsec));
    hold on;

    h = gobjects(0);
    labels = {};

    for ci = 1:numel(concList)
        conc = concList(ci);

        idx = (T.Conc_mg_mL == conc) & (T.Time_s == tsec);
        if nnz(idx) ~= 1
            error('Expected exactly one row in results for conc=%g, time=%g.', conc, tsec);
        end
        R = T(idx, :);
        p = [R.a, R.b, R.c, R.x0];

        fname = fullfile(dataDir, sprintf('%gmg_%gs_Background_Adjusted.xlsx', conc, tsec));
        if exist(fname,'file') ~= 2
            fname = fullfile(dataDir, sprintf('%dmg_%ds_Background_Adjusted.xlsx', round(conc), round(tsec)));
        end
        S = readLabFile(fname);

        x_raw = S.raw(:,1);
        x = correctTempsFromTable(x_raw, T_hotplate, T_actual);
        y = S.L;

        [x,y] = sanitizeXY(x,y);

        hx = plot(x, y, 'LineStyle','none', 'Marker','.', 'MarkerSize', ms);
        set(hx, 'Color', C(ci,:));
        h(end+1) = hx; %#ok<AGROW>
        labels{end+1} = sprintf('%g mg/mL Exp', conc); %#ok<AGROW>

        xg = linspace(50, 210, 800);
        yg = logisticModel(p, xg);

        hm = plot(xg, yg, '--', 'LineWidth', 2);
        set(hm, 'Color', C(ci,:));
        h(end+1) = hm; %#ok<AGROW>
        labels{end+1} = sprintf('%g mg/mL Model', conc); %#ok<AGROW>
    end

    xlim([50 210]);
    xlabel('Temperature (^{\circ}C)');
    ylabel('L^{*}');
    title(sprintf('Single Logistic Fits (%g s)', tsec));
    legend(h, labels, 'Location','best');
    grid on;
    grid minor;

    hold off;
end

end

function y = logisticModel(p, x)
a  = p(1);
b  = p(2);
c  = p(3);
x0 = p(4);
y = b + (a - b) ./ (1 + exp(-c .* (x - x0)));
end

function S = readLabFile(fname)
if exist(fname,'file') ~= 2
    error('Missing file: %s', fname);
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
x = x(:);
y = y(:);
n = min(numel(x), numel(y));
x = x(1:n);
y = y(1:n);

mask = isfinite(x) & isfinite(y);
x = x(mask);
y = y(mask);

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
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

set(0, ...
    'DefaultFigureColor', 'w', ...
    'DefaultAxesFontName', 'Helvetica', ...
    'DefaultAxesFontSize', 18, ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesBox', 'on', ...
    'DefaultAxesTickDir', 'out', ...
    'DefaultAxesXGrid', 'on', ...
    'DefaultAxesYGrid', 'on', ...
    'DefaultLineLineWidth', 2, ...
    'DefaultLegendBox', 'off');
end
