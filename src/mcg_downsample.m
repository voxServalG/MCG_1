function [Y, Fs_out] = mcg_downsample(X, Fs_in, Fs_out)
% 多级抗混叠降采样（默认 250k → 1k）。优先使用 resample；否则用 decimate 链 [5 5 5 2].

ratio = Fs_in / Fs_out;
Y = X;
if abs(ratio-1) < 1e-12
  Fs_out = Fs_in;  % 采样率相同，直接返回
  return;
end

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
