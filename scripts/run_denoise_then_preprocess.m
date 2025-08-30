function run_denoise_then_preprocess(main_mat, artifact_mat, varargin)
% RUN_DENOISE_THEN_PREPROCESS 先做伪影子空间投影(SSP)，再按PDF流程滤波/检波/叠加
% 用法：
%   run_denoise_then_preprocess('data/session01.mat','data/session01_artifact.mat', ...
%                               'RefChan', 1, 'BadChannels', [57 60 64]);

addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..','src')));

cfg = mcg_config();
p = inputParser; p.KeepUnmatched=true;
addParameter(p,'RefChan',cfg.refChan,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'BadChannels',[],@(x)isnumeric(x)&&isvector(x)); % 从子空间估计中排除
addParameter(p,'ArtifactWindows',[],@(x)isnumeric(x)&&size(x,2)==2);
parse(p,varargin{:});
if ~isempty(p.Results.RefChan),   cfg.refChan   = p.Results.RefChan; end
artifact_windows = p.Results.ArtifactWindows;

% --- 载入主数据/伪影数据（MAT）---
S_main = load(main_mat, 'X', 'Fs');
S_art  = load(artifact_mat, 'X', 'Fs');
assert(S_main.Fs == cfg.Fs_raw && S_art.Fs == cfg.Fs_raw, 'Fs 不一致或与配置不符');

X_main = double(S_main.X);
X_art  = double(S_art.X);

% 使用去基线后的数据来去伪影
[X_denoise, P, info_ssp] = mcg_artifact_ssp(detrend(X_main, 'constant'), detrend(X_art, 'constant'), varargin{:});

% --- 之后进入原有预处理（bandpass + baseline + R峰 + 叠加平均）---
[X_bp, filtInfo] = mcg_bandpass_baseline(X_denoise, cfg.Fs_raw, cfg);

x_ref = X_bp(:, cfg.refChan);
[rLocs, rScore] = mcg_detect_rpeaks(x_ref, cfg.Fs_raw, cfg);

[avgWave, epochs, t_epoch, keep_idx, qc] = mcg_epoch_average(X_bp, rLocs, cfg.Fs_raw, cfg, artifact_windows);

% --- 保存 ---
[~, b_main, ~] = fileparts(main_mat);
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'out');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
save(fullfile(out_dir, sprintf('denoise_preproc_%s.mat', b_main)), ...
    'X_bp', 'rLocs', 'rScore', 'epochs', 't_epoch', 'keep_idx', 'avgWave', 'filtInfo', 'cfg', 'qc', ...
    'P', 'info_ssp', '-v7.3');

fprintf('\n[OK] 去伪影→预处理完成 → out/denoise_preproc_%s.mat\n', b_main);
end
