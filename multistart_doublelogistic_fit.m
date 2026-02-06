function ResultsTbl = multistart_doublelogistic_fit(dataDir, outCsv)


close all; clc;
applyPublicationDefaults();

if nargin < 1 || isempty(dataDir)
    dataDir = pwd;
end
if exist(dataDir,'dir') ~= 7
    error('Data directory not found: %s', dataDir);
end
if nargin < 2 || isempty(outCsv)
    outCsv = 'ManualGlobalSearch_DoubleLogistic_FitResults.csv';
end
if exist('lsqcurvefit','file') ~= 2
    error('lsqcurvefit not found. Optimization Toolbox is required.');
end

rng(1);


T_hotplate = [90:2.5:140, 150:10:250];
T_actual   = [50:2.5:100, 110:10:210];


concList = [5 10 12 15 20];
timeList = [5 15 30 60];

D = struct();
for ci = 1:numel(concList)
    conc = concList(ci);
    keyC = sprintf('C%d', conc);
    D.(keyC) = struct();

    for ti = 1:numel(timeList)
        tsec = timeList(ti);
        keyT = sprintf('T%d', tsec);

        fname = fullfile(dataDir, sprintf('%dmg_%ds_Background_Adjusted.xlsx', conc, tsec));D.(keyC).(keyT) = readLabFile(fname);
    end
end


lb_base = [0    0    0.02   65    0     0.001   100];
ub_base = [100  100  3.00   105   100   0.15    300];


budgetSecPerDataset = 180;
coarseBatchN = 25000;
topK         = 60;
nearFrac     = 1.01;

useParallelLocal = true;
maxLocalIters    = 8000;
maxLocalEvals    = 2e5;

useWarmStart = true;
warmP = [];


Fit  = struct();
rows = [];

for ci = 1:numel(concList)
    conc = concList(ci);
    keyC = sprintf('C%d', conc);

    for ti = 1:numel(timeList)
        tsec = timeList(ti);
        keyT = sprintf('T%d', tsec);


        rawTemps = D.(keyC).(keyT).raw(:,1);
        x = correctTempsFromTable(rawTemps, T_hotplate, T_actual);
        y = D.(keyC).(keyT).L(:);
        [x,y] = sanitizeXY(x,y);


        [lb,ub] = boundsForDataset(lb_base, ub_base, x);


        out = fit_double_sigmoid_budgeted( ...
            x, y, lb, ub, ...
            budgetSecPerDataset, coarseBatchN, topK, nearFrac, ...
            warmP, useParallelLocal, maxLocalIters, maxLocalEvals);


        out = plateauBoundaryEscapeIfNeeded( ...
            x, y, lb, ub, out, ...
            budgetSecPerDataset, coarseBatchN, topK, nearFrac, ...
            warmP, useParallelLocal, maxLocalIters, maxLocalEvals);


        out = tempBoundaryExpandIfNeeded( ...
            x, y, lb, ub, out, ...
            budgetSecPerDataset, coarseBatchN, topK, nearFrac, ...
            warmP, useParallelLocal, maxLocalIters, maxLocalEvals);


        Fit.(keyC).(keyT).x = x;
        Fit.(keyC).(keyT).y = y;
        Fit.(keyC).(keyT).p = out.pbest;

        if useWarmStart
            warmP = out.pbest;
        end

        rows = [rows; ...
            conc, tsec, ...
            out.appParams.a, out.appParams.b, out.appParams.c, out.appParams.d, ...
            out.appParams.e, out.appParams.f, out.appParams.g, ...
            out.x01, out.x02, out.SSE, out.R2, out.exitflag, ...
            out.e_minNearBest, out.e_maxNearBest, out.nNearBest];

        fprintf('Done: %2d mg, %2d s | R2=%.5f | a=%.3f b=%.3f e=%.3f | x01=%.2f x02=%.2f | SSE=%.4g\n', ...
            conc, tsec, out.R2, out.appParams.a, out.appParams.b, out.appParams.e, out.x01, out.x02, out.SSE);
    end
end

ResultsTbl = array2table(rows, ...
    'VariableNames', {'Conc_mg_mL','Time_s', ...
                      'a','b','c','d','e','f','g', ...
                      'x01','x02','SSE','R2','exitflag', ...
                      'e_minNearBest','e_maxNearBest','nNearBest'});

disp(ResultsTbl);
writetable(ResultsTbl, fullfile(dataDir, outCsv));


concColors = lines(numel(concList));
expMarker  = 'o';
mSize      = 5;

for ti = 1:numel(timeList)
    tsec = timeList(ti);

    figure('Name', sprintf('%ds: Experimental (solid) vs Model (dashed)', tsec));
    hold on;

    for ci = 1:numel(concList)
        conc = concList(ci);
        keyC = sprintf('C%d', conc);
        keyT = sprintf('T%d', tsec);

        xExp = Fit.(keyC).(keyT).x;
        yExp = Fit.(keyC).(keyT).y;
        p    = Fit.(keyC).(keyT).p;

        plot(xExp, yExp, '-', ...
            'Color', concColors(ci,:), ...
            'Marker', expMarker, ...
            'MarkerSize', mSize, ...
            'MarkerFaceColor', concColors(ci,:), ...
            'MarkerEdgeColor', concColors(ci,:), ...
            'DisplayName', sprintf('%d mg/mL', conc));

        xFit = linspace(min(xExp), max(xExp), 600).';
        yFit = modelOnly(p, xFit);

        plot(xFit, yFit, '--', ...
            'Color', concColors(ci,:), ...
            'HandleVisibility','off');
    end

    xlabel('Temperature ($^{\circ}$C)');
    ylabel('Avg $L^*$');
    title(sprintf('%ds activation: solid = experimental, dashed = model', tsec));
    grid minor;
    ylim([0 100]);
    legend('Location','bestoutside');
end

end


function out = fit_double_sigmoid_budgeted(x, y, lb, ub, budgetSec, batchN, topK, nearFrac, warmP, useParallelLocal, maxIters, maxEvals)
x = x(:); y = y(:);

[b0, a0, e0, x01_0, x02_0, c0, f0] = dataDrivenSeeds(x,y);
seed1 = clampP([a0 b0 c0 x01_0 e0 f0 x02_0], lb, ub);
seed2 = clampP([min(max(a0+5,0),100) b0 0.6  (x01_0+3)  e0 0.01 (x02_0+5)], lb, ub);
seed3 = clampP([a0 b0 0.15 (x01_0-3) min(max(e0+5,0),100) 0.03 (x02_0-5)], lb, ub);

seeds = [seed1; seed2; seed3];
if ~isempty(warmP) && numel(warmP)==7 && all(isfinite(warmP))
    seeds = [clampP(warmP(:).', lb, ub); seeds];
end

Pbest = seeds;
Sbest = batchSSE(Pbest, x, y);

t0 = tic;
while toc(t0) < budgetSec
    Pbatch = generateRandomStarts(lb, ub, batchN);

    Pcand = [Pbest; Pbatch];
    Scand = [Sbest; batchSSE(Pbatch, x, y)];

    [Sbest, idx] = mink(Scand, min(topK, numel(Scand)));
    Pbest = Pcand(idx,:);

    if toc(t0) > 0.92*budgetSec
        break;
    end
end

nearMask = Sbest <= nearFrac*min(Sbest);
eNear = Pbest(nearMask,5);
if isempty(eNear)
    eMin = NaN; eMax = NaN; nNear = 0;
else
    eMin = min(eNear); eMax = max(eNear); nNear = numel(eNear);
end

opts = optimoptions('lsqcurvefit', ...
    'Display','off', ...
    'Algorithm','trust-region-reflective', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'MaxIterations', maxIters, ...
    'MaxFunctionEvaluations', maxEvals, ...
    'TypicalX',[80 50 0.3 median(x) 90 0.02 max(x)-20]);

bestSSE  = inf;
bestP    = nan(1,7);
bestExit = -inf;

if useParallelLocal && license('test','Distrib_Computing_Toolbox')
    K = size(Pbest,1);
    pFits = nan(K,7);
    sse   = inf(K,1);
    ef    = -inf(K,1);

    parfor k = 1:K
        p0 = Pbest(k,:);
        try
            [pfit, resnorm, ~, efk] = lsqcurvefit(@(p,xx) modelAndJacobian(p,xx), p0, x, y, lb, ub, opts);
            pFits(k,:) = pfit;
            sse(k)     = resnorm;
            ef(k)      = efk;
        catch
        end
    end

    [bestSSE, kbest] = min(sse);
    bestP    = pFits(kbest,:);
    bestExit = ef(kbest);
else
    for k = 1:size(Pbest,1)
        p0 = Pbest(k,:);
        try
            [pfit, resnorm, ~, ef] = lsqcurvefit(@(p,xx) modelAndJacobian(p,xx), p0, x, y, lb, ub, opts);
            if resnorm < bestSSE
                bestSSE  = resnorm;
                bestP    = pfit;
                bestExit = ef;
            end
        catch
        end
    end
end

if ~isfinite(bestSSE)
    error('All lsqcurvefit refinements failed. Check x/y finiteness and bounds.');
end

yhat = modelOnly(bestP,x);
SST = sum((y - mean(y)).^2);
R2  = 1 - bestSSE/max(SST, eps);

a   = bestP(1); b = bestP(2); c = bestP(3); x01 = bestP(4);
e   = bestP(5); f = bestP(6); x02 = bestP(7);
d   = c*(x01 - 50);
g   = f*(x02 - 50);

out.pbest = bestP;
out.SSE = bestSSE;
out.R2 = R2;
out.exitflag = bestExit;
out.appParams = struct('a',a,'b',b,'c',c,'d',d,'e',e,'f',f,'g',g);
out.x01 = x01;
out.x02 = x02;
out.e_minNearBest = eMin;
out.e_maxNearBest = eMax;
out.nNearBest = nNear;
end


function out2 = plateauBoundaryEscapeIfNeeded(x, y, lb, ub, out, budgetSec, batchN, topK, nearFrac, warmP, useParallelLocal, maxIters, maxEvals)
p = out.pbest;


tolL = 1e-2;
hitA = isNearBound(p(1), 0, 100, tolL);
hitB = isNearBound(p(2), 0, 100, tolL);
hitE = isNearBound(p(5), 0, 100, tolL);

if ~(hitA || hitB || hitE)
    out2 = out;
    return;
end

base = out;
baseSSE = base.SSE;


relSSEcap = 0.005;


margins = [1e-3 1e-2 0.05 0.1 0.25 0.5 1 2 5];


extraBudget = max(2*budgetSec, 240);
extraBatchN = max(batchN, 60000);
extraTopK   = max(topK, 150);
extraNear   = min(nearFrac, 1.005);

bestInterior = [];
bestInteriorSSE = inf;
bestMargin = NaN;


minSep = 10;
xmax = max(x);
xmin = min(x);


x02RightExtension = 2000;

for m = margins
    lb2 = lb; ub2 = ub;


    if hitA
        lb2(1) = max(lb2(1), 0 + m);
        ub2(1) = min(ub2(1), 100 - m);
    end
    if hitB
        lb2(2) = max(lb2(2), 0 + m);
        ub2(2) = min(ub2(2), 100 - m);
    end
    if hitE
        lb2(5) = max(lb2(5), 0 + m);
        ub2(5) = min(ub2(5), 100 - m);
    end

    if any(lb2 >= ub2)
        continue;
    end


    lb2(3) = min(lb2(3), 0.001);  ub2(3) = max(ub2(3), 8.0);
    lb2(6) = min(lb2(6), 0.0001); ub2(6) = max(ub2(6), 1.0);

    lb2(4) = max(0, min(lb2(4), xmin - 120));
    ub2(4) = max(ub2(4), xmax + 300);

    lb2(7) = max([lb2(7), lb2(4) + minSep]);
    ub2(7) = max(ub2(7), xmax + x02RightExtension);


    [lb2, ub2] = enforceBoundFeasibility(lb2, ub2);

    warm2 = base.pbest;

    try
        cand = fit_double_sigmoid_budgeted( ...
            x, y, lb2, ub2, ...
            extraBudget, extraBatchN, extraTopK, extraNear, ...
            warm2, useParallelLocal, maxIters, maxEvals);
    catch
        continue;
    end

    pc = cand.pbest;
    interiorOK = true;
    if hitA, interiorOK = interiorOK && ~isNearBound(pc(1), 0, 100, tolL); end
    if hitB, interiorOK = interiorOK && ~isNearBound(pc(2), 0, 100, tolL); end
    if hitE, interiorOK = interiorOK && ~isNearBound(pc(5), 0, 100, tolL); end

    if interiorOK && (cand.SSE <= (1 + relSSEcap)*baseSSE)
        out2 = cand;
        out2.plateau_escape = struct('used',true,'margin',m,'relSSE',cand.SSE/baseSSE,'keptInterior',true);
        return;
    end

    if cand.SSE < bestInteriorSSE
        bestInterior = cand; bestInteriorSSE = cand.SSE; bestMargin = m;
    end
end


if ~isempty(bestInterior) && isfinite(bestInteriorSSE) && bestInteriorSSE <= (1 + 3*relSSEcap)*baseSSE
    out2 = bestInterior;
    out2.plateau_escape = struct('used',true,'margin',bestMargin,'relSSE',bestInteriorSSE/baseSSE,'keptInterior',true);
else
    out2 = base;
    out2.plateau_escape = struct('used',true,'margin',NaN,'relSSE',1.0,'keptInterior',false);
end

end

function tf = isNearBound(v, lo, hi, tolAbs)
tf = (v <= lo + tolAbs) || (v >= hi - tolAbs);
end


function out2 = tempBoundaryExpandIfNeeded(x, y, lb, ub, out, budgetSec, batchN, topK, nearFrac, warmP, useParallelLocal, maxIters, maxEvals)
p = out.pbest;

tol = 1e-6;
hitX01 = (abs(p(4)-lb(4))<tol) || (abs(p(4)-ub(4))<tol);
hitX02 = (abs(p(7)-lb(7))<tol) || (abs(p(7)-ub(7))<tol);

if ~(hitX01 || hitX02)
    out2 = out;
    return;
end

lb2 = lb; ub2 = ub;

xmin = min(x);
xmax = max(x);
minSep = 10;


lb2(4) = max(0, xmin - 200);
ub2(4) = xmax + 400;

lb2(7) = max([lb2(7), lb2(4) + minSep, 0]);
ub2(7) = max(ub2(7), xmax + 2000);


lb2(3) = min(lb2(3), 0.001);  ub2(3) = max(ub2(3), 10.0);
lb2(6) = min(lb2(6), 0.0001); ub2(6) = max(ub2(6), 2.0);

[lb2, ub2] = enforceBoundFeasibility(lb2, ub2);

outB = fit_double_sigmoid_budgeted(x, y, lb2, ub2, max(0.75*budgetSec, 120), batchN, max(topK,80), nearFrac, out.pbest, useParallelLocal, maxIters, maxEvals);

if outB.SSE < out.SSE
    out2 = outB;
else
    out2 = out;
end
end


function [lb, ub] = boundsForDataset(lb_base, ub_base, x)
lb = lb_base;
ub = ub_base;

xmin = min(x);
xmax = max(x);
minSep = 10;


lb(4) = max([lb(4), 0, xmin - 50]);
ub(4) = max(ub(4), xmax + 150);

lb(7) = max([lb(7), 0, xmin + minSep]);
ub(7) = max(ub(7), xmax + 600);


lb(7) = max(lb(7), lb(4) + minSep);


[lb, ub] = enforceBoundFeasibility(lb, ub);

end

function [lb, ub] = enforceBoundFeasibility(lb, ub)

for k = 1:numel(lb)
    if lb(k) >= ub(k)
        mid = 0.5*(lb(k)+ub(k));
        lb(k) = mid - 1;
        ub(k) = mid + 1;
    end
end
end


function [yfit, J] = modelAndJacobian(p, x)
a = p(1); b = p(2); c = p(3); x01 = p(4);
e = p(5); f = p(6); x02 = p(7);

z1 = c.*(x - x01);
z2 = f.*(x - x02);

[S1, S1p] = safeSigmoidAndDeriv(z1);
[S2, S2p] = safeSigmoidAndDeriv(z2);

yfit = b + (a-b).*S1 + (e-a).*S2;

dya   = (S1 - S2);
dyb   = (1 - S1);
dyc   = (a-b) .* (x - x01) .* S1p;
dyx01 = -(a-b) .* c .* S1p;
dye   = S2;
dyf   = (e-a) .* (x - x02) .* S2p;
dyx02 = -(e-a) .* f .* S2p;

J = [dya, dyb, dyc, dyx01, dye, dyf, dyx02];
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


function SSE = batchSSE(P, x, y)
a   = P(:,1).';  b   = P(:,2).';  c   = P(:,3).';  x01 = P(:,4).';
e   = P(:,5).';  f   = P(:,6).';  x02 = P(:,7).';

X1 = (x - x01) .* c;
X2 = (x - x02) .* f;

X1 = max(min(X1, 60), -60);
X2 = max(min(X2, 60), -60);

S1 = 1 ./ (1 + exp(-X1));
S2 = 1 ./ (1 + exp(-X2));

Yhat = b + (a-b).*S1 + (e-a).*S2;
R = Yhat - y;
SSE = sum(R.^2, 1).';
SSE(~isfinite(SSE)) = inf;
end


function [b0, a0, e0, x01_0, x02_0, c0, f0] = dataDrivenSeeds(x,y)
b0 = median(y(x <= prctile(x,10)));

midMask = (x >= prctile(x,40)) & (x <= prctile(x,70));
if nnz(midMask) >= 5, a0 = median(y(midMask)); else, a0 = median(y); end

hiMask = (x >= prctile(x,90));
if nnz(hiMask) >= 5, e0 = median(y(hiMask)); else, e0 = max(y); end
e0 = min(max(e0,0),100);

[dL, xmid] = localDerivative(x,y);

x01_0 = pickPeakInWindow(xmid, dL, prctile(x,15), prctile(x,55), median(x));
x02_0 = pickPeakInWindow(xmid, dL, prctile(x,55), prctile(x,95), prctile(x,80));

c0 = 0.3;
f0 = 0.02;
end

function [dL, xmid] = localDerivative(x, y)
dx = diff(x); dy = diff(y);
dL = dy ./ dx;
xmid = 0.5*(x(1:end-1) + x(2:end));
mask = isfinite(dL) & isfinite(xmid);
dL = dL(mask); xmid = xmid(mask);
end

function x0 = pickPeakInWindow(xmid, dL, wlo, whi, fallback)
mask = (xmid >= wlo) & (xmid <= whi) & isfinite(dL);
if nnz(mask) < 3
    x0 = fallback; return;
end
[~, idx] = max(dL(mask));
xm = xmid(mask);
x0 = xm(idx);
end


function starts = generateRandomStarts(lb, ub, n)
U = rand(n,7);
starts = lb + U.*(ub-lb);


c_lb = max(lb(3), eps); c_ub = ub(3);
f_lb = max(lb(6), eps); f_ub = ub(6);
starts(:,3) = exp(log(c_lb) + rand(n,1).*(log(c_ub)-log(c_lb)));
starts(:,6) = exp(log(f_lb) + rand(n,1).*(log(f_ub)-log(f_lb)));


minSep = 10;
x01 = starts(:,4);
x02 = starts(:,7);

bad = x02 < (x01 + minSep);
if any(bad)
    push = minSep + 20*rand(sum(bad),1);
    x02_new = x01(bad) + push;
    starts(bad,7) = min(max(x02_new, lb(7)), ub(7));
end
end

function p = clampP(p, lb, ub)
p = max(p, lb);
p = min(p, ub);


minSep = 10;
if p(7) < p(4) + minSep
    p(7) = min(p(4) + minSep, ub(7));
end
end


function S = readLabFile(fname)
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
