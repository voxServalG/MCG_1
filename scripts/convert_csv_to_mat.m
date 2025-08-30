function convert_csv_to_mat(csv_path, out_mat_path, Fs, nChan, col_idx)
% 将 64 列 CSV（或含 >=64 列）转换为 .mat：变量 X [Nsamp x nChan], Fs
% 用法示例：
%   convert_csv_to_mat('../raw/data.csv','../data/session01.mat',250000,64,1:64)

if nargin < 5 || isempty(col_idx), col_idx = 1:nChan; end

% 读入（优先 readmatrix）
try
    A = readmatrix(csv_path);
catch
    T = readtable(csv_path);
    A = table2array(T);
end

assert(size(A,2) >= max(col_idx), 'CSV 列数不足：需要至少 %d 列', max(col_idx));

X = double(A(:, col_idx));
assert(size(X,2)==nChan, '选取列数与 nChan 不一致');

Fs = double(Fs); %#ok<NASGU>
save(out_mat_path, 'X', 'Fs', '-v7.3');

fprintf('[OK] 已保存 MAT：%s  （X:[%d x %d], Fs=%.0f)\n', out_mat_path, size(X,1), size(X,2), Fs);
end
