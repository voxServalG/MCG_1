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
  % 幅度 QC：参考通道（索引 1，可按需改）
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
