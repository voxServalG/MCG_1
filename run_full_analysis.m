function run_full_analysis()
% RUN_FULL_ANALYSIS  从原始数据开始做伪影去除、预处理与 P/QRS/T 特征提取
%
% 用法：
%   run_full_analysis()
%   （数据需位于 data/session01.mat 与 data/session01_artifact.mat）

%% --------- 1. 准备与载入 ---------
addpath(genpath('src'));            % 确保 src 路径已加入
cfg = mcg_config();                 % 读入默认参数

% 载入主信号与伪影信号
S_main = load('data/session01.mat');        % 包含 X (Nsamp x Nchan), Fs
S_art  = load('data/session01_artifact.mat'); 
assert(S_main.Fs == S_art.Fs, 'Fs mismatch');
Fs = S_main.Fs;

X_main = double(S_main.X);
X_art  = double(S_art.X);

%% --------- 2. SSP 伪影去除 ---------
[X_clean, P, info_ssp] = mcg_artifact_ssp(X_main, X_art, ...
    'varExplained', 0.90, 'kMax', 8);

%% --------- 3. 带通滤波 + 去基线 ---------
[X_bp, filtInfo] = mcg_bandpass_baseline(X_clean, Fs, cfg);

%% --------- 4. R 峰检测与心动周期切片 ---------
x_ref = X_bp(:, cfg.refChan);          % 参考通道
[rLocs, rScore] = mcg_detect_rpeaks(x_ref, Fs, cfg);

[avgWave, epochs, t_epoch, keep_idx, qc] = mcg_epoch_average( ...
    X_bp, rLocs, Fs, cfg, []);        % 若有伪影时间窗，可在 [] 处填入

%% --------- 5. P/QRS/T 波特征提取 ---------
t  = t_epoch;        % 相对 R 峰的时间轴（秒）
Ts = 1/Fs;

% 定义窗口：可根据实际需要微调
base_idx = t < -0.25;                        % 基线窗口
p_idx    = (t >= -0.12 & t <= -0.04);        % P 波
qrs_idx  = (t >= -0.04 & t <=  0.04);        % QRS
t_idx    = (t >=  0.10 & t <=  0.25);        % T 波

baseline = mean(avgWave(base_idx,:),1);      % 每通道基线
N        = size(avgWave,2);

feat = struct('P_amp',zeros(1,N),'P_dur',zeros(1,N), ...
              'QRS_amp',zeros(1,N),'QRS_dur',zeros(1,N), ...
              'T_amp',zeros(1,N),'T_dur',zeros(1,N));

for ch = 1:N
    % ----- P 波 -----
    p_w = avgWave(p_idx,ch) - baseline(ch);
    [pmax,~] = max(p_w);
    idx1 = find(p_w>0.1*pmax,1,'first');
    idx2 = find(p_w>0.1*pmax,1,'last');
    feat.P_amp(ch) = pmax;
    feat.P_dur(ch) = (idx2-idx1)*Ts;

    % ----- QRS -----
    qrs_w = avgWave(qrs_idx,ch) - baseline(ch);
    [qmax,~] = max(qrs_w); [qmin,~] = min(qrs_w);
    qpp = qmax - qmin;
    idx1 = find(abs(qrs_w)>0.1*qpp,1,'first');
    idx2 = find(abs(qrs_w)>0.1*qpp,1,'last');
    feat.QRS_amp(ch) = qpp;
    feat.QRS_dur(ch) = (idx2-idx1)*Ts;

    % ----- T 波 -----
    t_w = avgWave(t_idx,ch) - baseline(ch);
    [tmax,~] = max(t_w);
    idx1 = find(t_w>0.1*tmax,1,'first');
    idx2 = find(t_w>0.1*tmax,1,'last');
    feat.T_amp(ch) = tmax;
    feat.T_dur(ch) = (idx2-idx1)*Ts;
end

%% --------- 6. 保存与可视化 ---------
if ~exist('out','dir'), mkdir('out'); end
save('out/analysis_session01.mat', ...
    'X_bp','rLocs','t_epoch','epochs','avgWave','feat', ...
    'cfg','P','info_ssp','qc','filtInfo','rScore');

% 示例：绘制参考通道的平均波形
chan = cfg.refChan;
figure;
plot(t, avgWave(:,chan), 'LineWidth', 1.2); hold on;
yline(baseline(chan),'k--'); xline(0,'r');
title(sprintf('Channel %d Average Waveform', chan));
xlabel('Time (s)'); ylabel('Amplitude');
legend('Average','Baseline','R peak');

disp('Feature summary (channel-wise):');
disp(feat);
end
