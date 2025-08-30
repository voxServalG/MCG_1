
function quick_view_64()
%QUICK_VIEW_64  一键生成 64 通道视图（叠加平均 + 10s 段）
in  = fullfile('out','preproc_session01.mat');
png1 = fullfile('out','qc_grid64_avg_session01.png');
png2 = fullfile('out','qc_grid64_seg10s_session01.png');

qc_grid64(in, 'avg',     png1);
qc_grid64(in, 'segment', png2, 'Dur', 10, 'T0', 0);

fprintf('已生成：\n  %s\n  %s\n', png1, png2);
end
