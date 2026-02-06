function HeatResponse_Lab_Plots(tauInputsPath)
if nargin < 1 || isempty(tauInputsPath)
    tauInputsPath = 'tau_inputs.mat';
end

close all;
clc;

if exist(tauInputsPath,'file') ~= 2
    error('Missing file: %s', tauInputsPath);
end

S = load(tauInputsPath);

need = {'temps','timeVals','concVals','L_byTime'};
for k = 1:numel(need)
    if ~isfield(S, need{k})
        error('Missing field in %s: %s', tauInputsPath, need{k});
    end
end

temps = S.temps(:);
timeVals = S.timeVals(:);
concVals = S.concVals(:).';
L_byTime = S.L_byTime;

if isfield(S,'timeLabels')
    timeLabels = S.timeLabels;
else
    timeLabels = arrayfun(@(t)sprintf('%g',t), timeVals, 'UniformOutput', false);
end

if isfield(S,'concLabels')
    concLabels = S.concLabels;
else
    concLabels = arrayfun(@(c)sprintf('%g mg/mL',c), concVals, 'UniformOutput', false);
end

A_byTime = [];
B_byTime = [];

if isfield(S,'a_byTime')
    A_byTime = S.a_byTime;
elseif isfield(S,'A_byTime')
    A_byTime = S.A_byTime;
end

if isfield(S,'b_byTime')
    B_byTime = S.b_byTime;
elseif isfield(S,'B_byTime')
    B_byTime = S.B_byTime;
end

if isempty(A_byTime) || isempty(B_byTime)
    error('Missing a_byTime/b_byTime in %s', tauInputsPath);
end

nT = size(L_byTime,1);
nC = size(L_byTime,2);

if ~isequal(size(A_byTime), [nT nC]) || ~isequal(size(B_byTime), [nT nC])
    error('Size mismatch among L_byTime, a_byTime, b_byTime');
end

colors = lines(nC);

for ti = 1:nT
    ttl = sprintf('Heat Response L* (%s)', timeLabels{ti});
    figure('Name', ttl);
    hold on;
    for ci = 1:nC
        y = L_byTime{ti,ci};
        plot(temps, y, 'LineStyle','-', 'Marker','.', 'Color', colors(ci,:));
    end
    xlabel('Temperature (^\circC)');
    ylabel('L*');
    title(ttl);
    legend(concLabels, 'Location','best');
    grid on;
    grid minor;
    hold off;

    ttl = sprintf('Heat Response a* (%s)', timeLabels{ti});
    figure('Name', ttl);
    hold on;
    for ci = 1:nC
        y = A_byTime{ti,ci};
        plot(temps, y, 'LineStyle','-', 'Marker','.', 'Color', colors(ci,:));
    end
    xlabel('Temperature (^\circC)');
    ylabel('a*');
    title(ttl);
    legend(concLabels, 'Location','best');
    grid on;
    grid minor;
    hold off;

    ttl = sprintf('Heat Response b* (%s)', timeLabels{ti});
    figure('Name', ttl);
    hold on;
    for ci = 1:nC
        y = B_byTime{ti,ci};
        plot(temps, y, 'LineStyle','-', 'Marker','.', 'Color', colors(ci,:));
    end
    xlabel('Temperature (^\circC)');
    ylabel('b*');
    title(ttl);
    legend(concLabels, 'Location','best');
    grid on;
    grid minor;
    hold off;
end

end
