function response_table = UV_Response_Adjusted_Data_v2(dataDir)
if nargin < 1 || isempty(dataDir)
    dataDir = pwd;
end

close all;
clc;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:numel(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

pubFont   = 'Helvetica';
baseFS    = 18;
lw        = 2;
ms        = 20;

set(0, ...
    'DefaultFigureColor',                'w', ...
    'DefaultFigureUnits',                'inches', ...
    'DefaultFigurePosition',             [1 1 6 4], ...
    'DefaultAxesFontName',               pubFont, ...
    'DefaultAxesFontSize',               baseFS, ...
    'DefaultAxesLineWidth',              1, ...
    'DefaultAxesBox',                    'on', ...
    'DefaultAxesTickDir',                'out', ...
    'DefaultAxesTickLength',             [0.015 0.015], ...
    'DefaultAxesXGrid',                  'on', ...
    'DefaultAxesYGrid',                  'on', ...
    'DefaultAxesLabelFontSizeMultiplier', 1.1, ...
    'DefaultAxesTitleFontSizeMultiplier', 1.1, ...
    'DefaultLineLineWidth',              lw, ...
    'DefaultLineMarkerSize',             ms, ...
    'DefaultLegendBox',                  'off', ...
    'DefaultLegendFontSize',             baseFS);

concList = [5 10 12 15 20];
nC = numel(concList);

raw0 = readNumeric(fullfile(dataDir, sprintf('%dmg_UV_Background_Adjusted.xlsx', concList(1))));
if isempty(raw0) || size(raw0,2) < 4
    error('Invalid UV file for %d mg/mL.', concList(1));
end
timeVals = raw0(:,1);
nT = numel(timeVals);

Lmat = nan(nT, nC);
amat = nan(nT, nC);
bmat = nan(nT, nC);

for ci = 1:nC
    conc = concList(ci);
    fname = fullfile(dataDir, sprintf('%dmg_UV_Background_Adjusted.xlsx', conc));
    raw = readNumeric(fname);
    if isempty(raw) || size(raw,2) < 4
        error('File does not have >=4 numeric columns: %s', fname);
    end
    t = raw(:,1);
    if numel(t) ~= nT || any(abs(t - timeVals) > 0)
        error('Time column mismatch in %s', fname);
    end
    Lmat(:,ci) = raw(:,2);
    amat(:,ci) = raw(:,3);
    bmat(:,ci) = raw(:,4);
end

concColors = lines(nC);

plot2D(timeVals, Lmat, concList, concColors, 'Time (s)', 'L*', 'UV Response (L*)');
plot2D(timeVals, amat, concList, concColors, 'Time (s)', 'a*', 'UV Response (a*)');
plot2D(timeVals, bmat, concList, concColors, 'Time (s)', 'b*', 'UV Response (b*)');

plot3Dcurves(timeVals, Lmat, concList, concColors, 'Time (s)', 'Concentration (mg/mL)', 'L*',  'UV Response 3D (L*)');
plot3Dcurves(timeVals, amat, concList, concColors, 'Time (s)', 'Concentration (mg/mL)', 'a*',  'UV Response 3D (a*)');
plot3Dcurves(timeVals, bmat, concList, concColors, 'Time (s)', 'Concentration (mg/mL)', 'b*',  'UV Response 3D (b*)');

row1  = 1;
row32 = min(32, nT);

L0  = Lmat(row1,:);
L75 = Lmat(row32,:);
delta_L = L75 - L0;
abs_delta_L = abs(delta_L);

response_table = table( ...
    concList(:), delta_L(:), abs_delta_L(:), ...
    'VariableNames', {'Concentration_mg_ml', 'DeltaL_0to75', 'Abs_DeltaL_0to75'});

disp('Magnitude of L* response from 0 to 75 s (using rows 1 and 32 where available):');
disp(response_table);

figure('Name','L* Response Magnitude 0–75 s');
plot(concList, abs_delta_L, 'LineStyle', '-', 'Marker','.', 'Color','k');
xlabel('Concentration (mg/mL)');
ylabel('Magnitude of Response');
grid minor;
grid on;

end

function plot2D(x, Y, concList, concColors, xlab, ylab, ttl)
figure('Name', ttl);
hold on;
for ci = 1:numel(concList)
    plot(x, Y(:,ci), 'LineStyle','-', 'Marker','.', 'Color', concColors(ci,:));
end
xlabel(xlab);
ylabel(ylab);
title(ttl);
legend(arrayfun(@(c)sprintf('%g mg/mL',c), concList, 'UniformOutput', false), 'Location','best');
grid minor;
grid on;
hold off;
end

function plot3Dcurves(x, Y, concList, concColors, xlab, ylab, zlab, ttl)
figure('Name', ttl);
hold on;
for ci = 1:numel(concList)
    y = concList(ci) .* ones(size(x));
    plot3(x, y, Y(:,ci), 'LineStyle','-', 'Marker','.', 'Color', concColors(ci,:));
end
xlabel(xlab);
ylabel(ylab);
zlabel(zlab);
title(ttl);
grid on;
grid minor;
view(40, 22);
hold off;
end

function raw = readNumeric(fname)
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
raw = raw(~all(isnan(raw),2), :);
end
