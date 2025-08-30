function convert_binary_to_mat(bin_path, out_mat_path, Fs, nChan, dtype, header_bytes, endianness)
% 将自定义二进制时序文件转换为 .mat（支持超大文件分块写入）
% 假设连续按通道交错存储（C-order）：s1_ch1, s1_ch2, ... s1_chN, s2_ch1, ...，以此类推
% 参数：
%   bin_path     : 源二进制文件路径（例如 raw.baseData）
%   out_mat_path : 输出 .mat 路径（例如 ../data/session01.mat）
%   Fs           : 采样率（如 250000）
%   nChan        : 通道数（如 64）
%   dtype        : 数据类型（'single' | 'int16' | ...）
%   header_bytes : 文件头长度（字节数，如 512；若无头则 0）
%   endianness   : 字节序（'ieee-le' 小端 | 'ieee-be' 大端）
%
% 用法示例：
%   convert_binary_to_mat('../raw/rec.baseData','../data/session01.mat',250000,64,'single',512,'ieee-le')

arguments
  bin_path (1,:) char
  out_mat_path (1,:) char
  Fs (1,1) double {mustBePositive}
  nChan (1,1) double {mustBeInteger, mustBePositive}
  dtype (1,:) char
  header_bytes (1,1) double {mustBeNonnegative}
  endianness (1,:) char {mustBeMember(endianness,{'ieee-le','ieee-be','l','b'})}
end

% matfile 方式分块写（避免 2GB 限制 & 内存占用过大）
m = matfile(out_mat_path, 'Writable', true);
m.Fs = double(Fs);

% 打开文件并跳过头
fid = fopen(bin_path, 'r', endianness);
assert(fid>0, '无法打开文件：%s', bin_path);
cleanup = onCleanup(@() fclose(fid));
fseek(fid, header_bytes, 'bof');

% 每块样本数（每块大小 ≈ nSampChunk*nChan*bytesPerSample）
bytesPer = bytes_per_sample(dtype);
targetChunkBytes = 64*1024*1024; % 64MB 左右，可按需调整
nSampChunk = floor(targetChunkBytes / (nChan * bytesPer));
nSampChunk = max(nSampChunk, 1e5); % 至少 1e5 个样本/块以提高效率

row = 1;
total_samples = 0;
while true
    % 读一块
    count = nSampChunk * nChan;
    [raw, got] = fread(fid, count, ['*' dtype]);
    if got==0
        break;
    end
    ns = floor(got / nChan);
    if ns==0
        break;
    end
    raw = raw(1:ns*nChan); % 丢弃不足一帧的尾部
    Xc  = reshape(raw, nChan, ns).'; % [ns x nChan]
    % 写入 .mat（double）
    m.X(row:row+ns-1, 1:nChan) = double(Xc);
    row = row + ns;
    total_samples = total_samples + ns;
end

fprintf('[OK] 已生成 %s  (X:[%d x %d], Fs=%.0f)\n', out_mat_path, total_samples, nChan, Fs);

end % function

function n = bytes_per_sample(dtype)
switch lower(dtype)
  case {'single','float32'}, n = 4;
  case {'double','float64'}, n = 8;
  case {'int16','uint16'},   n = 2;
  case {'int32','uint32'},   n = 4;
  case {'int8','uint8','char'}, n = 1;
  otherwise, error('未支持的数据类型：%s', dtype);
end
end
