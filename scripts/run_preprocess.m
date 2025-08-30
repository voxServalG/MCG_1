function run_preprocess(input_path, varargin)
% 一键运行：原始数据 → 降采样 → 带通/去基线 → R峰 → 心动周期 → 叠加平均 → 保存
% 用法：run_preprocess(fullfile('..','data','session01.mat'))
% 可选：run_preprocess(..., 'RefChan', 3, 'TargetFs', 1000)

addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..','src')));

cfg = mcg_config();
% 覆盖可选参数
p = inputParser; p.KeepUnmatched=true;
addParameter(p,'RefChan',cfg.refChan,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'TargetFs',cfg.targetFs,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'ArtifactWindows',[],@(x)isnumeric(x)&&size(x,2)==2); % [t_start t_end] 秒
parse(p,varargin{:});
if ~isempty(p.Results.RefChan),   cfg.refChan   = p.Results.RefChan; end
if ~isempty(p.Results.TargetFs),  cfg.targetFs  = p.Results.TargetFs; end
artifact_windows = p.Results.ArtifactWindows;

% 1) 读取
[X_raw, Fs_raw, meta] = mcg_io_load(input_path, cfg);  % X_raw: [Nsamples x Nchan]
assert(Fs_raw==cfg.Fs_raw, 'Fs_raw与配置不一致：%.3f vs %.3f', Fs_raw, cfg.Fs_raw);

% 2) 抗混叠降采样到 cfg.targetFs (默认 1000 Hz)
[X_ds, Fs_ds] = mcg_downsample(X_raw, Fs_raw, cfg.targetFs);

% 3) 0.5–45 Hz 零相位带通 + 去基线
[X_bp, filtInfo] = mcg_bandpass_baseline(X_ds, Fs_ds, cfg);

% 4) 参考通道 R 峰检测
x_ref = X_bp(:, cfg.refChan);
[rLocs, rScore] = mcg_detect_rpeaks(x_ref, Fs_ds, cfg);

% 5) 切片、质控、叠加平均
[avgWave, epochs, t_epoch, keep_idx, qc] = mcg_epoch_average(X_bp, rLocs, Fs_ds, cfg, artifact_windows);

% 6) 保存
[~,base,~] = fileparts(input_path);
out_dir = fullfile(fileparts(mfilename('fullpath')),'..','out');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
save(fullfile(out_dir, sprintf('preproc_%s.mat', base)), ...
    'X_bp','Fs_ds','rLocs','rScore','epochs','t_epoch','keep_idx','avgWave','filtInfo','meta','cfg','qc','-v7.3');

fprintf('\n[OK] 预处理完成 → out/preproc_%s.mat\n', base);
end
