
function qc_grid64(preproc_mat, mode, out_png, varargin)
%QC_GRID64  以 8x8 子图显示 64 通道波形，便于快速体检。
% 用法：
%   % 1) 叠加平均视图（推荐，用 avgWave & t_epoch）
%   qc_grid64('out/preproc_session01.mat', 'avg', 'out/qc_grid64_avg_session01.png');
%
%   % 2) 原始预处理波形的前10秒（X_bp）
%   qc_grid64('out/preproc_session01.mat', 'segment', 'out/qc_grid64_seg10s_session01.png', 'Dur', 10, 'T0', 0);
%
% 选项（Name-Value）：
%   'Dur' (double)  : segment 模式下的时长（秒），默认 10
%   'T0'  (double)  : segment 模式的起始时间（秒），默认 0
%   'CLimPct' (double) : 振幅裁剪百分位（对称），默认 99（稳健缩放）
%
% 输出：保存 PNG（若 out_png 为空，则仅显示不保存）

p = inputParser;
addRequired(p, 'preproc_mat', @(s)ischar(s)||isstring(s));
addRequired(p, 'mode', @(s)ischar(s)||isstring(s));
addRequired(p, 'out_png', @(s)ischar(s)||isstring(s)||isempty(s));
addParameter(p, 'Dur', 10, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p, 'T0', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p, 'CLimPct', 99, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<100);
parse(p, preproc_mat, mode, out_png, varargin{:});
opt = p.Results;

S = load(preproc_mat);
mode = lower(string(mode));

if mode == "avg"
    assert(isfield(S,'avgWave') && isfield(S,'t_epoch'), '缺少 avgWave/t_epoch');
    X = S.avgWave;        % [L x 64]
    t = S.t_epoch(:);     % [L x 1]
    xlab = 'Time (s)  (relative to R)';
    ttl  = '64通道叠加平均 (avgWave)';
elseif mode == "segment"
    assert(isfield(S,'X_bp') && isfield(S,'Fs_ds'), '缺少 X_bp/Fs_ds');
    Fs = S.Fs_ds;
    t0 = opt.T0;
    dur = min(opt.Dur, (size(S.X_bp,1)/Fs) - t0);
    idx0 = max(1, floor(t0*Fs)+1);
    idx1 = min(size(S.X_bp,1), idx0 + floor(dur*Fs) - 1);
    X = S.X_bp(idx0:idx1, :);
    t = ( (idx0:idx1) - 1 )'/Fs;
    xlab = 'Time (s)';
    ttl  = sprintf('64通道预处理波形片段 [%.2f, %.2f] s', t(1), t(end));
else
    error('未知模式：%s（支持 avg / segment）', mode);
end

if size(X,2) ~= 64
    warning('通道数=%d（非64），仍尝试绘图。', size(X,2));
end

% --- 绘图 ---
figure('Color','w','Position',[100 100 1400 1200]);
tiledlayout(8,8,'TileSpacing','compact','Padding','compact');

nChan = size(X,2);
for ch = 1:nChan
    nexttile;
    y = X(:,ch);
    % 稳健缩放：以百分位限定 y 轴（对称）
    pctl = prctile(abs(y), opt.CLimPct);
    if pctl <= 0 || ~isfinite(pctl)
        pctl = max(1e-12, std(y)*3);
    end
    plot(t, y,'-'); grid on
    ylim([-pctl pctl]);
    title(sprintf('Ch%02d', ch), 'FontSize', 8);
    if ch > (nChan-8)
        xlabel(xlab);
    end
end
sgtitle(ttl, 'FontWeight','bold');

if ~isempty(out_png)
    out_dir = fileparts(out_png);
    if ~isempty(out_dir) && ~exist(out_dir,'dir'), mkdir(out_dir); end
    exportgraphics(gcf, out_png, 'Resolution', 160);
    fprintf('[OK] 已输出：%s\n', out_png);
end
end
