close all;
clear;
clc;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

pubFont   = 'Helvetica';
baseFS    = 18;
titleFS   = 13;
lw        = 2;
ms        = 20;

set(0, ...
    'DefaultFigureColor',              'w', ...
    'DefaultFigureUnits',              'inches', ...
    'DefaultFigurePosition',           [1 1 6 4], ...
    'DefaultAxesFontName',             pubFont, ...
    'DefaultAxesFontSize',             baseFS, ...
    'DefaultAxesLineWidth',            1, ...
    'DefaultAxesBox',                  'on', ...
    'DefaultAxesTickDir',              'out', ...
    'DefaultAxesTickLength',           [0.015 0.015], ...
    'DefaultAxesXGrid',                'on', ...
    'DefaultAxesYGrid',                'on', ...
    'DefaultAxesLabelFontSizeMultiplier', 1.1, ...
    'DefaultAxesTitleFontSizeMultiplier', 1.1, ...
    'DefaultLineLineWidth',            lw, ...
    'DefaultLineMarkerSize',           ms, ...
    'DefaultLegendBox',                'off', ...
    'DefaultLegendFontSize',           baseFS);

concColors = lines(5);
concMarkers = {'.','.','.','.','.'};

fiveMg_import = importdata('5mg_UV_Background_Adjusted.xlsx');
fiveMg_data   = fiveMg_import.data;

fiveMg_L = fiveMg_data(:,2);
fiveMg_a = fiveMg_data(:,3);
fiveMg_b = fiveMg_data(:,4);

tenMg_import = importdata('10mg_UV_Background_Adjusted.xlsx');
tenMg_data   = tenMg_import.data;

tenMg_L = tenMg_data(:,2);
tenMg_a = tenMg_data(:,3);
tenMg_b = tenMg_data(:,4);

twelveMg_import = importdata('12mg_UV_Background_Adjusted.xlsx');
twelveMg_data   = twelveMg_import.data;

twelveMg_L = twelveMg_data(:,2);
twelveMg_a = twelveMg_data(:,3);
twelveMg_b = twelveMg_data(:,4);

fifteenMg_import = importdata('15mg_UV_Background_Adjusted.xlsx');
fifteenMg_data   = fifteenMg_import.data;

fifteenMg_L = fifteenMg_data(:,2);
fifteenMg_a = fifteenMg_data(:,3);
fifteenMg_b = fifteenMg_data(:,4);

twentyMg_import = importdata('20mg_UV_Background_Adjusted.xlsx');
twentyMg_data   = twentyMg_import.data;

twentyMg_L = twentyMg_data(:,2);
twentyMg_a = twentyMg_data(:,3);
twentyMg_b = twentyMg_data(:,4);

timeVals = fiveMg_data(:,1);

figure('Name','UV Response (L*)');
hold on;
plot(timeVals, fiveMg_L, 'LineStyle','-', 'Marker', concMarkers{1});
plot(timeVals, tenMg_L, 'LineStyle','-', 'Marker', concMarkers{2});
plot(timeVals, twelveMg_L, 'LineStyle','-', 'Marker', concMarkers{3});
plot(timeVals, fifteenMg_L, 'LineStyle','-', 'Marker', concMarkers{4});
plot(timeVals, twentyMg_L, 'LineStyle','-', 'Marker', concMarkers{5});
xlabel('Time (s)');
ylabel('L*');
title('UV Response (L*)');
legend({'5 mg/mL','10 mg/mL','12 mg/mL','15 mg/mL','20 mg/mL'}, 'Location','best');
grid minor;
grid on;
hold off;

figure('Name','UV Response (a*)');
hold on;
plot(timeVals, fiveMg_a, 'LineStyle','-', 'Marker', concMarkers{1});
plot(timeVals, tenMg_a, 'LineStyle','-', 'Marker', concMarkers{2});
plot(timeVals, twelveMg_a, 'LineStyle','-', 'Marker', concMarkers{3});
plot(timeVals, fifteenMg_a, 'LineStyle','-', 'Marker', concMarkers{4});
plot(timeVals, twentyMg_a, 'LineStyle','-', 'Marker', concMarkers{5});
xlabel('Time (s)');
ylabel('a*');
title('UV Response (a*)');
legend({'5 mg/mL','10 mg/mL','12 mg/mL','15 mg/mL','20 mg/mL'}, 'Location','best');
grid minor;
grid on;
hold off;

figure('Name','UV Response (b*)');
hold on;
plot(timeVals, fiveMg_b, 'LineStyle','-', 'Marker', concMarkers{1});
plot(timeVals, tenMg_b, 'LineStyle','-', 'Marker', concMarkers{2});
plot(timeVals, twelveMg_b, 'LineStyle','-', 'Marker', concMarkers{3});
plot(timeVals, fifteenMg_b, 'LineStyle','-', 'Marker', concMarkers{4});
plot(timeVals, twentyMg_b, 'LineStyle','-', 'Marker', concMarkers{5});
xlabel('Time (s)');
ylabel('b*');
title('UV Response (b*)');
legend({'5 mg/mL','10 mg/mL','12 mg/mL','15 mg/mL','20 mg/mL'}, 'Location','best');
grid minor;
grid on;
hold off;

row1  = 1;
row32 = 32;

L0  = [fiveMg_L(row1),  tenMg_L(row1),  twelveMg_L(row1),  fifteenMg_L(row1),  twentyMg_L(row1)];
L75 = [fiveMg_L(row32), tenMg_L(row32), twelveMg_L(row32), fifteenMg_L(row32), twentyMg_L(row32)];

delta_L = L75 - L0;
abs_delta_L = abs(delta_L);

concs_mg = [5, 10, 12, 15, 20];

response_table = table( ...
    concs_mg(:), delta_L(:), abs_delta_L(:), ...
    'VariableNames', {'Concentration_mg_ml', 'DeltaL_0to75', 'Abs_DeltaL_0to75'});

disp('Magnitude of L* response from 0 to 75 s (using rows 1 and 32):');
disp(response_table);

figure('Name','L* Response Magnitude 0–75 s');
plot(concs_mg, abs_delta_L, 'LineStyle', '-', 'Marker','.');
xlabel('Concentration (mg/mL)');
ylabel('Magnitude of Response');
grid minor;
grid on;
