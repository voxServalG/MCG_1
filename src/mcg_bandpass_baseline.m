
function [Y, info] = mcg_bandpass_baseline(X, Fs, cfg)
% PDF版预处理：1–40 Hz 带通(零相位) + 300点中值去基线
% 参考：指导书 3.2 中 (a) 与 (b) 步骤

% 1) 1–40 Hz 带通（与 PDF 参数一致：Steepness 0.85, StopbandAttenuation 60 dB）
Y1 = bandpass(X, [1 40], Fs, 'Steepness', 0.85, 'StopbandAttenuation', 60);

% 2) 300 窗口中值滤波估计基线并扣除（沿样本维度）
base = medfilt1(Y1, 300);

Y = Y1 - base;

info = struct('bpBand',[1 40],'steepness',0.85,'stopAtten',60,'medfiltWin',300);
end
