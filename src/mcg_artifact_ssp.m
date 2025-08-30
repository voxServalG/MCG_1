function [X_clean, P, info] = mcg_artifact_ssp(X_main, X_art, varargin)
% MCG_ARTIFACT_SSP  使用“伪影记录”估计噪声子空间，并对主数据做信号空间投影(SSP)
% 输入：
%   X_main [Ns x N] : 主数据（未滤波，仅做过 detrend/去DC 即可）
%   X_art  [Ma x N] : 伪影记录（与主数据同通道数）
%   cfg.artifact.varExplained ∈ (0,1] : 子空间累计方差阈值（默认 0.90）
%   cfg.artifact.kMax         : 上限保留的分量数（默认 8）
%   cfg.artifact.excludeChan  : 需要从子空间估计中排除的通道索引（可选）
%
% 输出：
%   X_clean [Ns x N] : 去伪影后的主数据（仍在原始采样率域）
%   P [N x N]        : 投影矩阵 (I - U_k*U_k^T)
%   info             : 结构体，含奇异值、选取的 k、方差占比等

% 使用 inputParser 处理输入参数（兼容旧版 MATLAB）
p = inputParser;
addRequired(p, 'X_main', @isnumeric);
addRequired(p, 'X_art', @isnumeric);
addParameter(p, 'varExplained', 0.90, @(x) isnumeric(x) && x > 0 && x <= 1);
addParameter(p, 'kMax', 8, @(x) isnumeric(x) && x > 0);
addParameter(p, 'excludeChan', [], @isnumeric);
parse(p, X_main, X_art, varargin{:});

cfg.artifact.varExplained = p.Results.varExplained;
cfg.artifact.kMax = p.Results.kMax;
cfg.artifact.excludeChan = p.Results.excludeChan;

[Ns, N] = size(X_main);
assert(size(X_art, 2) == N, '通道数不一致：主数据=%d, 伪影=%d', N, size(X_art, 2));

% --- 去均值（按通道）---
Xm = X_main - mean(X_main, 1);
Za = X_art - mean(X_art, 1);

% --- 可选：剔除某些通道参与子空间估计（如坏道）---
incl = true(1, N);
if ~isempty(cfg.artifact.excludeChan)
    incl(cfg.artifact.excludeChan) = false;
end
Za_used = Za(:, incl);

% --- 估计伪影的通道协方差 & SVD ---
C = (Za_used.' * Za_used) / max(1, size(Za_used, 1) - 1);
[U, Sv] = svd(C, 'econ');
sing = diag(Sv);

% 方差占比 & 选 k
vr = sing / max(eps, sum(sing));
cvr = cumsum(vr);
k = find(cvr >= cfg.artifact.varExplained, 1, 'first');
if isempty(k), k = 1; end
k = min(k, min(cfg.artifact.kMax, size(U, 2)));

% 构造完整尺寸的 U_full (N x k)
U_full = zeros(N, k);
U_full(incl, :) = U(:, 1:k);

% 投影矩阵（I - U*U^T），对主数据做空间投影
P = eye(N) - (U_full * U_full.');
X_clean = Xm * P;

info = struct('singular', sing, 'var_ratio', vr, 'cumvar', cvr, ...
              'k', k, 'incl', find(incl), 'excluded', find(~incl));
end
