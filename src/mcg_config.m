function cfg = mcg_config()
% 全局参数（基于 Fs_raw = 250 kHz, Nchan = 64）
cfg = struct();

% 原始数据规格
cfg.Fs_raw    = 1000;    % 原始采样率（Hz）
cfg.nChan     = 64;      % 通道数

% 目标降采样规格
cfg.targetFs  = 1000;    % 目标采样率（Hz）（与原始相同将自动跳过降采样）

% 预处理滤波
cfg.bpBand    = [1 40];  % 与PDF一致
cfg.bpOrder   = 4;         % 保留但未使用（bandpass内部自动设计）

% R 峰最小间隔（秒），PDF=0.5，可按需要改小
cfg.minRR_s   = 0.5;

% 参考通道与 R 峰检测
cfg.refChan     = 1;        % 可改为最清晰的通道索引
cfg.rBand       = [5 15];   % R 峰增强带通（Hz）
cfg.rMinRR      = 0.30;     % MinPeakDistance ≈ 0.30 s（对应 >200 bpm 的上限）
cfg.rSmoothWin  = 0.10;     % 包络平滑窗口（秒）
% 鲁棒阈值参数（对单位/幅值不敏感）
cfg.rHeightFrac = 0.35;     % MinPeakHeight = p50 + rHeightFrac*(p99-p50)
cfg.rPromFrac   = 0.25;     % MinPeakProminence = rPromFrac*(p99-p50)

% 心动周期切片（相对 R 峰）
cfg.epochWin    = [-0.33 0.67]; % 秒（≈ P→T 全覆盖）

% 质控开关（为与PDF一致，默认关闭）
cfg.qc.enable  = false;

% 质控（MAD/阈值）
cfg.qc.maxRRDev = 0.30;    % RR 偏差容忍（相对中位数）
cfg.qc.zMaxAbs  = 8;       % 振幅绝对 Z 分数上限（单通道）

end
