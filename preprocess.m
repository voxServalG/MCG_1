# 目录结构（建议）

```
mcg_pipeline/
├─ data/                     % 原始数据放这里（.mat 或自定义 .baseData）
├─ out/                      % 处理结果输出
├─ scripts/
│   └─ run_preprocess.m      % 一键运行脚本（入口）
└─ src/
    ├─ mcg_config.m          % 全局参数配置
    ├─ mcg_io_load.m         % 数据读取（.mat 已支持，.baseData 留接口）
    ├─ mcg_downsample.m      % 多级抗混叠降采样（250kHz→1kHz）
    ├─ mcg_bandpass_baseline.m % 0.5–45 Hz 零相位带通+去基线
    ├─ mcg_detect_rpeaks.m   % R 峰检测（参考通道）
    ├─ mcg_epoch_average.m   % 心动周期切片、质控、叠加平均
    └─ mcg_utils.m           % 工具函数（可选）
```

---

## scripts/run_preprocess.m

```matlab
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
```

---

## src/mcg_config.m

```matlab
function cfg = mcg_config()
% 全局参数（基于 Fs_raw = 250 kHz, Nchan = 64）
cfg = struct();

% 原始数据规格
cfg.Fs_raw    = 250e3;   % 原始采样率（Hz）
cfg.nChan     = 64;      % 通道数

% 目标降采样规格
cfg.targetFs  = 1000;    % 目标采样率（Hz）

% 预处理滤波
cfg.bpBand    = [0.5 45];  % 带通（Hz），兼顾 P/QRS/T
cfg.bpOrder   = 4;         % IIR 阶数（Butterworth via designfilt）

% 参考通道与 R 峰检测
cfg.refChan     = 1;        % 可改为最清晰的通道索引
cfg.rBand       = [5 15];   % R 峰增强带通（Hz）
cfg.rMinRR      = 0.30;     % MinPeakDistance ≈ 0.30 s（对应 >200 bpm 的上限）
cfg.rThreshSD   = 2.5;      % 峰值阈值：均值 + rThreshSD×std
cfg.rSmoothWin  = 0.10;     % 包络平滑窗口（秒）

% 心动周期切片（相对 R 峰）
cfg.epochWin    = [-0.33 0.67]; % 秒（≈ P→T 全覆盖）

% 质控（MAD/阈值）
cfg.qc.maxRRDev = 0.2;     % RR 偏差容忍（相对中位数）
cfg.qc.zMaxAbs  = 6;       % 振幅绝对 Z 分数上限（单通道）

end
```

---

## src/mcg_io_load.m

```matlab
function [X, Fs, meta] = mcg_io_load(input_path, cfg)
% 读取 .mat 或自定义二进制（留接口）
% 输出：X [Nsamp x Nchan]（double），Fs（Hz）

[~,~,ext] = fileparts(input_path);
meta = struct('path',input_path,'note','');

switch lower(ext)
  case '.mat'
    S = load(input_path);
    if isfield(S,'X') && isfield(S,'Fs')
      X = double(S.X); Fs = S.Fs;
    else
      error('.mat 需包含变量 X [Nsamp x Nchan], Fs');
    end

  otherwise
    % TODO: 如果是 .baseData，请按你们的格式在此实现：
    % 例如：header_bytes=512; dtype='single'; 行主序/列主序等。
    error('暂未实现此扩展名的读取：%s（请先转存为 .mat）', ext);
end

% 基本一致性检查
assert(size(X,2)==cfg.nChan, '通道数不符：数据=%d, 期望=%d', size(X,2), cfg.nChan);
X = detrend(X,'constant'); % 去 DC
end
```

---

## src/mcg_downsample.m

```matlab
function [Y, Fs_out] = mcg_downsample(X, Fs_in, Fs_out)
% 多级抗混叠降采样（默认 250k → 1k）。优先使用 resample；否则用 decimate 链 [5 5 5 2]。

ratio = Fs_in / Fs_out; assert(abs(ratio-round(ratio))<1e-6, '需为整数倍率');
Y = X;

if exist('resample','file')
  % 一步降到目标（内置抗混叠 FIR）
  Y = resample(Y, 1, ratio);
else
  % fallback：分级 decimate（每级自带低通）
  stages = factorize_ratio(ratio);
  for k = 1:numel(stages)
    Y = decimate(Y, stages(k));
  end
end

% 小工具：把比率因子化为 [5 5 5 2] 之类（尽量 2/3/5 的乘积）
  function s = factorize_ratio(r)
    s = [];
    for p = [5 3 2]
      while mod(r,p)==0; s(end+1)=p; r=r/p; end %#ok<AGROW>
    end
    assert(r==1, '无法仅用 2/3/5 分解比率，请安装 Signal Processing Toolbox 以使用 resample');
  end
end
```

---

## src/mcg_bandpass_baseline.m

```matlab
function [Y, info] = mcg_bandpass_baseline(X, Fs, cfg)
% 0.5–45 Hz 零相位带通；额外去基线（高通已足够，保留接口）

bp = designfilt('bandpassiir','FilterOrder',cfg.bpOrder, ...
  'HalfPowerFrequency1',cfg.bpBand(1), 'HalfPowerFrequency2',cfg.bpBand(2), ...
  'DesignMethod','butter','SampleRate',Fs);

Y = filtfilt(bp, X);   % 零相位

% 可选：中值/低阶多项式去基线（通常带通已覆盖 <0.5 Hz 漂移）
% for ch=1:size(Y,2)
%   base = medfilt1(Y(:,ch), round(0.3*Fs));
%   Y(:,ch) = Y(:,ch) - base;
% end

info = struct('bpBand',cfg.bpBand,'bpOrder',cfg.bpOrder);
end
```

---

## src/mcg_detect_rpeaks.m

```matlab
function [rLocs, score] = mcg_detect_rpeaks(x_ref, Fs, cfg)
% 简化 Pan–Tompkins 思路：5–15 Hz → 绝对值平方 → 平滑 → 自适应阈值 → findpeaks

rb = designfilt('bandpassiir','FilterOrder',4, ...
  'HalfPowerFrequency1',cfg.rBand(1),'HalfPowerFrequency2',cfg.rBand(2), ...
  'DesignMethod','butter','SampleRate',Fs);
xf = filtfilt(rb, x_ref);

y  = xf.^2;
win = max(1, round(cfg.rSmoothWin*Fs));
y  = movmean(y, win);

mu = mean(y); sd = std(y);
th = mu + cfg.rThreshSD*sd;

[~, rLocs] = findpeaks(y, 'MinPeakDistance', round(cfg.rMinRR*Fs), ...
                           'MinPeakHeight', th);

score = struct('mu',mu,'sd',sd,'th',th,'nPeaks',numel(rLocs));
end
```

---

## src/mcg_epoch_average.m

```matlab
function [avgWave, epochs, t_epoch, keep_idx, qc] = mcg_epoch_average(X, rLocs, Fs, cfg, artifact_windows)
% 以 R 峰为中心切片，质控剔除异常搏动，叠加平均

if nargin<5 || isempty(artifact_windows), artifact_windows = zeros(0,2); end

N = size(X,2);
pre  = round(abs(cfg.epochWin(1))*Fs);
post = round(abs(cfg.epochWin(2))*Fs);
L = pre + post + 1;

t_epoch = (-pre:post)'/Fs;  % 相对 R 峰的时间轴（秒）

% 1) 剔除越界与伪影重叠的心搏
valid = rLocs(rLocs>pre & rLocs+post<=size(X,1));
if ~isempty(artifact_windows)
  t = (0:size(X,1)-1)'/Fs;
  keep = true(size(valid));
  for i=1:numel(valid)
    seg = [t(valid(i)-pre), t(valid(i)+post)];
    if any(overlap(seg, artifact_windows))
      keep(i) = false;
    end
  end
  valid = valid(keep);
end

% 2) 构造 epoch: [L x N x K]
K = numel(valid);
epochs = zeros(L, N, K);
for k=1:K
  idx = (valid(k)-pre):(valid(k)+post);
  epochs(:,:,k) = X(idx,:);
end

% 3) 质控：RR 与幅度 Z 分数
RR = diff(valid)/Fs; medRR = median(RR);
keep = true(1,K);
for k=2:K-1  % 跳过首尾（缺 RR）
  if abs(RR(k-1)-medRR)>cfg.qc.maxRRDev*medRR || abs(RR(k)-medRR)>cfg.qc.maxRRDev*medRR
    keep(k) = false; continue;
  end
  % 幅度 QC：参考通道
  seg = epochs(:,1,k); z = (seg-mean(seg))/std(seg+eps);
  if any(abs(z)>cfg.qc.zMaxAbs)
    keep(k) = false;
  end
end
keep_idx = find(keep);

% 4) 叠加平均（对剩余搏动）
if isempty(keep_idx)
  warning('无有效心搏通过质控');
  avgWave = mean(epochs,3,'omitnan');
else
  avgWave = mean(epochs(:,:,keep_idx), 3);
end

qc = struct('nEpochs',K, 'nKept',numel(keep_idx), 'keepIdx',keep_idx);

% ---- 内部：区间重叠检查 ----
  function tf = overlap(seg, win)
    % seg: 1x2, win: Mx2
    tf = any( seg(1) < win(:,2) & seg(2) > win(:,1) );
  end
end
```

---

## src/mcg_utils.m（可选）

```matlab
function y = zscore_per_channel(X)
mu = mean(X,1); sd = std(X,[],1); y = (X - mu) ./ (sd + eps);
end
```

---

# 运行示例

```matlab
% 1) 将原始数据放到 mcg_pipeline/data/session01.mat
%    变量：X [Nsamp x 64], Fs = 250000
% 2) 打开 MATLAB，cd 到 mcg_pipeline/scripts
% 3) 运行：
run_preprocess(fullfile('..','data','session01.mat'), 'RefChan', 1, 'TargetFs', 1000);
% 结果：mcg_pipeline/out/preproc_session01.mat
%        含 X_bp（预处理波形，1kHz）、rLocs（R 峰）、epochs/t_epoch、avgWave（64 通道叠加平均）、质控统计 qc
