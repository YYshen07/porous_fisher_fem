function draw_3d_4()
% draw_3d_5
% 5号脚本：基于 core_0 的结果画3D图
% - 参数用 PARAMS 结构数组
% - 对每个参数 & 时间步画完整3D和半平面3D
% - 每组参数按 L{Lx}x{Ly}_afr{0p***}_flag{iflag} 存图

%% ===== 1. 参数配置 =====
PARAMS = [ ...
    struct('Lx',40,'Ly',40,'afr',0.125,'iflag',1) ...
    % 需要的话可以继续加：
    % , struct('Lx',60,'Ly',60,'afr',0.25,'iflag',1)
];

inum_values = [300, 900, 2100];  % 时间步
h_          = 600;               % 高度缩放
%% ======================

for k = 1:numel(PARAMS)
    p    = PARAMS(k);
    Lx   = p.Lx;
    Ly   = p.Ly;
    afr  = p.afr;
    flag = p.iflag;

    % 参数目录
    outdir = param_dirname(Lx, Ly, afr, flag);
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    fprintf('\n===== draw_3d_5: Lx=%d, Ly=%d, afr=%.3f, iflag=%d =====\n', ...
            Lx, Ly, afr, flag);

    for t = 1:numel(inum_values)
        inum = inum_values(t);

        % 调用 core_0 得到 z,x,y（不再调用 test_one6）
        [~, ~, z, x, y] = core_0(Lx, Ly, afr, inum, flag);

        % 放大坐标，和原来 test_one6 一样 *10
        X = x * 10;
        Y = y * 10;
        [Xg, Yg] = meshgrid(X, Y);   % z 是 n×m，对应 Y 行 X 列

        % ik = 0: 完整3D图
        % ik = 1: 半平面3D图
        for ik = 0:1
            fig = figure('Position', [100, 100, 800, 600]);
            set(fig, 'Color', 'none');
            ah = axes;
            hold(ah, 'on');
            box(ah, 'on');

            % surf 绘制体
            surf(ah, Xg, Yg, (z-1)*h_, 'edgecolor', 'k');

            set(ah,'FontName','Times New Roman','FontWeight','bold');

            if ik == 0
                %% 完整3D图
                view(ah, 3);
                xlim([0 Lx*10])
                ylim([0 Ly*10])
                zlim([-h_ 0])

                % 边框线
                plot3(ah,[0 0],[0 0],[0 -h_],'-k',...
                         [Lx*10 Lx*10],[0 0],[0 -h_],'-k',...
                         [Lx*10 Lx*10],[Ly*10 Ly*10],[0 -h_],'-k',...
                         [0 0],[Ly*10 Ly*10],[0 -h_],'-k',...
                         [0 Lx*10 Lx*10 0 0],[0 0 Ly*10 Ly*10 0],[0 0 0 0 0],'-k',...
                         [0 Lx*10 Lx*10 0 0],[0 0 Ly*10 Ly*10 0],-h_*ones(1,5),'-k');
            else
                %% 半平面3D图（后半部分）
                view(ah, [-30, 20]);
                xlim([0 Lx*10])
                ylim([Ly*5 Ly*10])
                zlim([-h_ 0])

                % 外边框
                plot3(ah,[0 0],[Ly*5 Ly*5],[0 -h_],'-k',...
                         [Lx*10 Lx*10],[Ly*5 Ly*5],[0 -h_],'-k',...
                         [Lx*10 Lx*10],[Ly*10 Ly*10],[0 -h_],'-k',...
                         [0 0],[Ly*10 Ly*10],[0 -h_],'-k',...
                         [0 Lx*10 Lx*10 0 0],[Ly*5 Ly*5 Ly*10 Ly*10 Ly*5],[0 0 0 0 0],'-k',...
                         [0 Lx*10 Lx*10 0 0],[Ly*5 Ly*5 Ly*10 Ly*10 Ly*5],-h_*ones(1,5),'-k');

                % 红色剖面框
                plot3(ah,[0 Lx*10],[Ly*5 Ly*5],[0 0],'-r','LineWidth',2);
                plot3(ah,[0 Lx*10],[Ly*5 Ly*5],[-h_ -h_],'-r','LineWidth',2);
                plot3(ah,[0 0],[Ly*5 Ly*5],[0 -h_],'-r','LineWidth',2);
                plot3(ah,[Lx*10 Lx*10],[Ly*5 Ly*5],[0 -h_],'-r','LineWidth',2);
            end

            % 去掉坐标轴和网格
            axis(ah, 'off');
            grid(ah, 'off');
            set(ah, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
            set(ah, 'XLabel', [], 'YLabel', [], 'ZLabel', []);
            title(ah, '');

            % 文件名
            if ik == 0
                pngfile = fullfile(outdir, sprintf('t%d_3D_full.png', inum));
            else
                pngfile = fullfile(outdir, sprintf('t%d_3D_half.png', inum));
            end

            exportgraphics(fig, pngfile, 'Resolution', 300);
            close(fig);

            fprintf('已保存：%s\n', pngfile);
        end
    end
end

fprintf('\n所有 3D 图绘制完毕。\n');
end


%% 参数目录命名
function outdir = param_dirname(Lx, Ly, afr, flag)
afr_str = strrep(num2str(afr, '%.3f'), '.', 'p');  % 0.125 -> 0p125
outdir  = sprintf('L%dx%d_afr%s_flag%d', Lx, Ly, afr_str, flag);
end