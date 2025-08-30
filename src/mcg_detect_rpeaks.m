
function [rLocs, score] = mcg_detect_rpeaks(x_ref, Fs, cfg)
% PDF版 R 波检测（健壮化）：仅用最小峰距，但自动限制到数据长度
% 默认最小峰距 = cfg.minRR_s * Fs（PDF建议 0.5s）

% 目标最小峰距（样本）
minDist_goal = round(cfg.minRR_s * Fs);

% 实际允许的最大值：必须 < length(x_ref)
maxAllowed = max(1, numel(x_ref)-1);
minDist = min(minDist_goal, maxAllowed);

% 若数据极短（<2*minDist），进一步放宽到长度的 1/3
if numel(x_ref) < 2*minDist
    minDist = max(1, floor(maxAllowed/3));
end

[~, rLocs] = findpeaks(x_ref, 'MinPeakDistance', minDist);

score = struct('minDist_goal',minDist_goal,'minDist',minDist,'len',numel(x_ref),'nPeaks',numel(rLocs));
end
