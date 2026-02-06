function ResultsTbl = multistart_singlelogistic_fit(dataDir, outCsv)
if nargin < 1 || isempty(dataDir)
    dataDir = pwd;
end
if nargin < 2 || isempty(outCsv)
    outCsv = 'ManualGlobalSearch_SingleLogistic_FitResults.csv';
end

close all;
clc;
rng(1);

T_hotplate = [90:2.5:140, 150:10:250];
T_actual   = [50:2.5:100, 110:10:210];

concList = [5 10 12 15 20];
timeList = [5 15 30 60];

lb_base = [0    0    0.02   50];
ub_base = [100  100  5.00   250];

budgetCandidates = 30000;
topK = 60;
maxLocalIters = 8000;
maxLocalEvals = 2e5;

rows = [];

for ci = 1:numel(concList)
    conc = concList(ci);

    for ti = 1:numel(timeList)
        tsec = timeList(ti);

        fname = fullfile(dataDir, sprintf('%dmg_%ds_Background_Adjusted.xlsx', conc, tsec));
        S = readLabFile(fname);

        x_raw = S.raw(:,1);
        x = correctTempsFromTable(x_raw, T_hotplate, T_actual);
        y = S.L;

        [x,y] = sanitizeXY(x,y);

        xMin = min(x);
        xMax = max(x);

        lb = lb_base;
        ub = ub_base;
        lb(4) = max(lb(4), xMin - 10);
        ub(4) = min(ub(4), xMax + 10);
        if ub(4) <= lb(4)
            lb(4) = xMin;
            ub(4) = xMax;
        end

        p0 = roughInit(x,y,lb,ub);

        pBest = [];
        sseBest = inf;

        pCand = sampleCandidates(budgetCandidates, lb, ub, p0);
        sseCand = nan(size(pCand,1),1);
        for k = 1:size(pCand,1)
            yhat = logisticModel(pCand(k,:), x);
            r = yhat - y;
            sseCand(k) = sum(r.^2);
        end

        [sseSorted, idx] = sort(sseCand, 'ascend');
        idx = idx(1:min(topK, numel(idx)));

        opts = optimoptions('lsqcurvefit', ...
            'Display','off', ...
            'MaxIterations', maxLocalIters, ...
            'MaxFunctionEvaluations', maxLocalEvals);

        for k = 1:numel(idx)
            pStart = pCand(idx(k),:);
            try
                pFit = lsqcurvefit(@(p,xv) logisticModel(p,xv), pStart, x, y, lb, ub, opts);
                yhat = logisticModel(pFit, x);
                r = yhat - y;
                sse = sum(r.^2);
                if sse < sseBest
                    sseBest = sse;
                    pBest = pFit;
                end
            catch
            end
        end

        if isempty(pBest)
            pBest = p0;
            yhat = logisticModel(pBest, x);
            r = yhat - y;
            sseBest = sum(r.^2);
        end

        rmse = sqrt(sseBest / numel(y));

        rows = [rows; conc, tsec, pBest(:).', rmse, sseBest];
    end
end

ResultsTbl = array2table(rows, 'VariableNames', ...
    {'Conc_mg_mL','Time_s','a','b','c','x0','RMSE','SSE'});

writetable(ResultsTbl, fullfile(dataDir, outCsv));

end

function y = logisticModel(p, x)
a  = p(1);
b  = p(2);
c  = p(3);
x0 = p(4);
y = b + (a - b) ./ (1 + exp(-c .* (x - x0)));
end

function p0 = roughInit(x,y,lb,ub)
a0 = max(y);
b0 = min(y);
yMid = (a0 + b0) / 2;
[~,iMid] = min(abs(y - yMid));
x0 = x(iMid);

if a0 < b0
    tmp = a0; a0 = b0; b0 = tmp;
end

c0 = 0.15;

p0 = [a0 b0 c0 x0];
p0 = max(p0, lb);
p0 = min(p0, ub);
end

function P = sampleCandidates(N, lb, ub, p0)
d = numel(lb);
P = nan(N, d);

P(1,:) = p0;

U = rand(N-1, d);
for j = 1:d
    P(2:end,j) = lb(j) + (ub(j)-lb(j)).*U(:,j);
end

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