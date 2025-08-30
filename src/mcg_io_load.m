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
