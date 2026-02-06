function out = Heat_Response_LateX_Data(tauInputs)
if nargin < 1 || isempty(tauInputs)
    tauInputs = 'tau_inputs.mat';
end
if ischar(tauInputs) || isstring(tauInputs)
    f = char(tauInputs);
    if exist(f,'file') ~= 2
        error('Required input file not found: %s', f);
    end
    S = load(f, 'temps','L_byTime','dLdT_byTime','timeVals','timeLabels','timeShort','concVals','concLabels');
elseif isstruct(tauInputs)
    S = tauInputs;
else
    error('First argument must be a .mat filepath or a struct with required fields.');
end
close all; clc;
temps       = S.temps;
L_byTime    = S.L_byTime;
dLdT_byTime = S.dLdT_byTime;
timeVals    = S.timeVals(:);
timeLabels  = S.timeLabels;
timeShort   = S.timeShort;
concVals    = S.concVals(:).';
concLabels  = S.concLabels;


tauVec  = 0.05:0.01:0.25;
tauPick = 0.09;


primaryRule = 'maxSensitivity';


useAbsSlope = true;
smoothSpan  = 1;


minDeltaT = 2.5;
minMagDL  = 0;


pubFont   = 'Helvetica';
baseFS    = 11;
lw        = 1.8;
ms        = 4;

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
    'DefaultLineLineWidth',            lw, ...
    'DefaultLineMarkerSize',           ms, ...
    'DefaultLegendBox',                'on', ...
    'DefaultLegendFontSize',           baseFS);

concColors  = lines(numel(concVals));
concMarkers = {'o','s','^','d','v'};


assert(isvector(temps) && isnumeric(temps), 'temps must be a numeric vector');
assert(isequal(size(L_byTime),[numel(timeVals) numel(concVals)]), 'L_byTime size mismatch');
assert(isequal(size(dLdT_byTime),[numel(timeVals) numel(concVals)]), 'dLdT_byTime size mismatch');


nT   = numel(timeVals);
nC   = numel(concVals);
nTau = numel(tauVec);

MOR_all  = nan(nT,nC,nTau);
DR_all   = nan(nT,nC,nTau);
Sens_all = nan(nT,nC,nTau);

for kt = 1:nTau
    tau = tauVec(kt);

    [MOR_tau, DR_tau, Sens_tau] = computeMetricsForTau( ...
        temps, L_byTime, dLdT_byTime, tau, ...
        primaryRule, useAbsSlope, smoothSpan, minDeltaT, minMagDL);

    MOR_all(:,:,kt)  = MOR_tau;
    DR_all(:,:,kt)   = DR_tau;
    Sens_all(:,:,kt) = Sens_tau;
end


[~,kPick] = min(abs(tauVec - tauPick));
tauPickActual = tauVec(kPick);

MOR  = MOR_all(:,:,kPick);
DR   = DR_all(:,:,kPick);
Sens = Sens_all(:,:,kPick);


MOR_5mg  = MOR(:,1);  MOR_10mg = MOR(:,2);  MOR_12mg = MOR(:,3);  MOR_15mg = MOR(:,4);  MOR_20mg = MOR(:,5);
DR_5mg   = DR(:,1);   DR_10mg  = DR(:,2);   DR_12mg  = DR(:,3);   DR_15mg  = DR(:,4);   DR_20mg  = DR(:,5);
Sens_5mg = Sens(:,1); Sens_10mg= Sens(:,2); Sens_12mg= Sens(:,3); Sens_15mg= Sens(:,4); Sens_20mg= Sens(:,5);


figure('Name',sprintf('Magnitude of Response (Tau-based, tau=%.2f)', tauPickActual));
plot(timeVals, MOR_5mg,  '-', 'Color',concColors(1,:), 'Marker',concMarkers{1}, 'MarkerFaceColor',concColors(1,:)); hold on;
plot(timeVals, MOR_10mg, '-', 'Color',concColors(2,:), 'Marker',concMarkers{2}, 'MarkerFaceColor',concColors(2,:));
plot(timeVals, MOR_12mg, '-', 'Color',concColors(3,:), 'Marker',concMarkers{3}, 'MarkerFaceColor',concColors(3,:));
plot(timeVals, MOR_15mg, '-', 'Color',concColors(4,:), 'Marker',concMarkers{4}, 'MarkerFaceColor',concColors(4,:));
plot(timeVals, MOR_20mg, '-', 'Color',concColors(5,:), 'Marker',concMarkers{5}, 'MarkerFaceColor',concColors(5,:));
xlim([min(timeVals) max(timeVals)]);
xlabel('UV Time (s)');
ylabel('Magnitude of Response ($|\Delta L^*|$)');
legend(concLabels,'Location','northeast');
grid minor;


figure('Name',sprintf('Dynamic Range (Tau-based, tau=%.2f)', tauPickActual));
plot(timeVals, DR_5mg,  '-', 'Color',concColors(1,:), 'Marker',concMarkers{1}, 'MarkerFaceColor',concColors(1,:)); hold on;
plot(timeVals, DR_10mg, '-', 'Color',concColors(2,:), 'Marker',concMarkers{2}, 'MarkerFaceColor',concColors(2,:));
plot(timeVals, DR_12mg, '-', 'Color',concColors(3,:), 'Marker',concMarkers{3}, 'MarkerFaceColor',concColors(3,:));
plot(timeVals, DR_15mg, '-', 'Color',concColors(4,:), 'Marker',concMarkers{4}, 'MarkerFaceColor',concColors(4,:));
plot(timeVals, DR_20mg, '-', 'Color',concColors(5,:), 'Marker',concMarkers{5}, 'MarkerFaceColor',concColors(5,:));
xlim([min(timeVals) max(timeVals)]);
xlabel('UV Time (s)');
ylabel('Dynamic Range ($\Delta T$ ($^{\circ}$C))');
legend(concLabels,'Location','northeast');
grid minor;


figure('Name',sprintf('Sensitivity (Tau-based, tau=%.2f)', tauPickActual));
plot(timeVals, Sens_5mg,  '-', 'Color',concColors(1,:), 'Marker',concMarkers{1}, 'MarkerFaceColor',concColors(1,:)); hold on;
plot(timeVals, Sens_10mg, '-', 'Color',concColors(2,:), 'Marker',concMarkers{2}, 'MarkerFaceColor',concColors(2,:));
plot(timeVals, Sens_12mg, '-', 'Color',concColors(3,:), 'Marker',concMarkers{3}, 'MarkerFaceColor',concColors(3,:));
plot(timeVals, Sens_15mg, '-', 'Color',concColors(4,:), 'Marker',concMarkers{4}, 'MarkerFaceColor',concColors(4,:));
plot(timeVals, Sens_20mg, '-', 'Color',concColors(5,:), 'Marker',concMarkers{5}, 'MarkerFaceColor',concColors(5,:));
xlim([min(timeVals) max(timeVals)]);
xlabel('UV Time (s)');
ylabel('Sensitivity ($|\Delta L^*|/\Delta T$)');
legend(concLabels,'Location','northeast');
grid minor;


for ci = 1:nC
    figure('Name', sprintf('Tau Sweep: Sensitivity heatmap (%s)', concLabels{ci}));
    M = squeeze(Sens_all(:,ci,:)).';   % [tau x time]
    imagesc(timeVals, tauVec, M);
    set(gca,'YDir','normal');
    xlabel('UV Time (s)');
    ylabel('\tau threshold');
    title(sprintf('Sensitivity heatmap for %s (%s)', concLabels{ci}, primaryRule), 'FontSize', baseFS);
    colorbar;
end


ti15 = find(timeVals == 15, 1);
if ~isempty(ti15)
    fprintf('\nPeak sensitivity tau by concentration at 15 s:\n');
    for ci = 1:nC
        sCurve = squeeze(Sens_all(ti15,ci,:));
        [sMax, idx] = max(sCurve, [], 'omitnan');
        fprintf('  %s: max Sens=%.4g at tau=%.2f\n', concLabels{ci}, sMax, tauVec(idx));
    end
end


[MOR_tbl, DR_tbl, Sens_tbl, T_start_tbl, T_end_tbl, L_start_tbl, L_end_tbl] = computeTableFromPrimaryWindow( ...
    temps, L_byTime, dLdT_byTime, tauPickActual, primaryRule, minDeltaT, minMagDL);


PCDA = repmat(concVals(:), nT, 1);
UV   = reshape(repmat(timeVals(:), 1, nC), [], 1);

T12  = reshape(T_start_tbl, [], 1);
T23  = reshape(T_end_tbl,   [], 1);
L12  = reshape(L_start_tbl, [], 1);
L23  = reshape(L_end_tbl,   [], 1);
MORv = reshape(MOR_tbl,     [], 1);
DRv  = reshape(DR_tbl,      [], 1);
Sv   = reshape(Sens_tbl,    [], 1);

HeatMetrics = table(PCDA, UV, MORv, DRv, Sv, T12, T23, L12, L23, ...
    'VariableNames', {'PCDA_mg_mL','UV_s','MOR','DR','Sensitivity','T_start_C','T_end_C','L_start','L_end'});

outCsv = sprintf('HeatMetrics_primaryWindow_tau%.2f.csv', tauPickActual);
writetable(HeatMetrics, outCsv);
fprintf('\nWrote %s\n', outCsv);


fprintf('\nLaTeX rows (PCDA & UV & MOR & DR & Sens) for tau=%.2f:\n', tauPickActual);
for ci = 1:nC
    for ti = 1:nT
        mor = MOR_tbl(ti,ci);
        dr  = DR_tbl(ti,ci);
        ss  = Sens_tbl(ti,ci);

        if ~isfinite(mor), morStr = '--'; else, morStr = sprintf('%.3f', mor); end
        if ~isfinite(dr),  drStr  = '--'; else, drStr  = sprintf('%.2f',  dr);  end
        if ~isfinite(ss),  sStr   = '--'; else, sStr   = sprintf('%.3f', ss);  end

        fprintf('%g & %g & %s & %s & %s \\\\\n', concVals(ci), timeVals(ti), morStr, drStr, sStr);
    end
end


out = struct();
out.tauVec = tauVec;
out.tauPick = tauPick;
out.primaryRule = primaryRule;
out.tauPickActual = tauPickActual;
out.MOR_all = MOR_all;
out.DR_all = DR_all;
out.Sens_all = Sens_all;
if exist('HeatMetrics','var')
    out.HeatMetrics = HeatMetrics;
end
end

function [MOR_tau, DR_tau, Sens_tau] = computeMetricsForTau( ...
        temps, L_byTime, dLdT_byTime, tau, primaryRule, useAbsSlope, smoothSpan, minDeltaT, minMagDL)

    nT = size(L_byTime,1);
    nC = size(L_byTime,2);

    DeltaL = nan(nT,nC);
    DeltaT = nan(nT,nC);

    for ti = 1:nT
        for ci = 1:nC
            Lvec = L_byTime{ti,ci};
            dvec = dLdT_byTime{ti,ci};

            if numel(Lvec) ~= numel(temps) || numel(dvec) ~= (numel(temps)-1)
                continue;
            end

            winTbl = tauWindowsFromDerivative(temps, Lvec, dvec, tau, 'leftpoint');

            if isempty(winTbl)
                continue;
            end


            keep = isfinite(winTbl.DeltaT_C) & isfinite(winTbl.MagDeltaL) & ...
                   (winTbl.DeltaT_C >= minDeltaT) & (winTbl.MagDeltaL >= minMagDL);

            if ~any(keep)
                continue;
            end

            winTbl = winTbl(keep,:);


            switch lower(primaryRule)
                case 'maxsensitivity'
                    s = winTbl.MagDeltaL ./ winTbl.DeltaT_C;
                    [~,k] = max(s);
                case 'longestdeltat'
                    [~,k] = max(winTbl.DeltaT_C);
                otherwise
                    error('Unknown primaryRule: %s', primaryRule);
            end

            DeltaL(ti,ci) = winTbl.DeltaL(k);
            DeltaT(ti,ci) = winTbl.DeltaT_C(k);
        end
    end

    MOR_tau = abs(DeltaL);
    DR_tau  = DeltaT;

    Sens_tau = MOR_tau ./ DR_tau;
    Sens_tau(~isfinite(Sens_tau)) = NaN;
end

function [MOR_tbl, DR_tbl, Sens_tbl, T_start_tbl, T_end_tbl, L_start_tbl, L_end_tbl] = computeTableFromPrimaryWindow( ...
        temps, L_byTime, dLdT_byTime, tau, primaryRule, minDeltaT, minMagDL)


    nT = size(L_byTime,1);
    nC = size(L_byTime,2);

    MOR_tbl  = nan(nT,nC);
    DR_tbl   = nan(nT,nC);
    Sens_tbl = nan(nT,nC);

    T_start_tbl = nan(nT,nC);
    T_end_tbl   = nan(nT,nC);
    L_start_tbl = nan(nT,nC);
    L_end_tbl   = nan(nT,nC);

    for ti = 1:nT
        for ci = 1:nC
            Lvec = L_byTime{ti,ci};
            dvec = dLdT_byTime{ti,ci};

            if numel(Lvec) ~= numel(temps) || numel(dvec) ~= (numel(temps)-1)
                continue;
            end

            winTbl = tauWindowsFromDerivative(temps, Lvec, dvec, tau, 'leftpoint');
            if isempty(winTbl)
                continue;
            end


            keep = isfinite(winTbl.DeltaT_C) & isfinite(winTbl.MagDeltaL) & ...
                   (winTbl.DeltaT_C >= minDeltaT) & (winTbl.MagDeltaL >= minMagDL);

            if ~any(keep)
                continue;
            end
            winTbl = winTbl(keep,:);


            switch lower(primaryRule)
                case 'maxsensitivity'
                    s = winTbl.MagDeltaL ./ winTbl.DeltaT_C;
                    [~,k] = max(s);
                case 'longestdeltat'
                    [~,k] = max(winTbl.DeltaT_C);
                otherwise
                    error('Unknown primaryRule: %s', primaryRule);
            end

            T_start_tbl(ti,ci) = winTbl.T_start_C(k);
            T_end_tbl(ti,ci)   = winTbl.T_end_C(k);
            L_start_tbl(ti,ci) = winTbl.L_start(k);
            L_end_tbl(ti,ci)   = winTbl.L_end(k);

            DR_tbl(ti,ci)  = winTbl.DeltaT_C(k);
            MOR_tbl(ti,ci) = abs(winTbl.DeltaL(k));
            Sens_tbl(ti,ci)= MOR_tbl(ti,ci) ./ DR_tbl(ti,ci);
        end
    end

    Sens_tbl(~isfinite(Sens_tbl)) = NaN;
end

function winTbl = tauWindowsFromDerivative(temps, L, dLdT, tau, mode)


    if nargin < 5 || isempty(mode)
        mode = 'leftpoint';
    end

    temps = temps(:);
    L     = L(:);
    dLdT  = dLdT(:);

    above = isfinite(dLdT) & (dLdT > tau);

    segStart = find(diff([false; above]) == 1);
    segEnd   = find(diff([above; false]) == -1);

    if isempty(segStart)
        winTbl = table();
        return;
    end

    switch lower(mode)
        case 'leftpoint'
            T_start = temps(segStart);
            T_end   = temps(segEnd);

            L_start = L(segStart);
            L_end   = L(segEnd);

        case 'interval'
            T_start = temps(segStart);
            T_end   = temps(segEnd + 1);

            L_start = L(segStart);
            L_end   = L(segEnd + 1);

        otherwise
            error('Unknown mode: %s. Use ''leftpoint'' or ''interval''.', mode);
    end

    DeltaT  = T_end - T_start;
    DeltaL  = L_end - L_start;
    MagL    = abs(DeltaL);

    winTbl = table(T_start, T_end, DeltaT, L_start, L_end, DeltaL, MagL, ...
        'VariableNames', {'T_start_C','T_end_C','DeltaT_C','L_start','L_end','DeltaL','MagDeltaL'});
end
