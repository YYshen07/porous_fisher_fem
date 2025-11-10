function feature_time_1()
% feature_time_1
% 计算给定参数组的特征时间（孔内所有单元 >= 0.8 的最早时间步）
% 思路：
%   - 对每一组参数做一次时间推进（最多 max_inum 步）
%   - 每一步检查是否“全场 >= threshold”
%   - 对有随机性的 iflag=4/5/6，会在多组随机种子上重复，直到获得若干个有效样本

%% ===== 1. 配置参数组（按需修改） =====
PARAMS = [ 
    struct('Lx',40,'Ly',40,'afr',0.25,'iflag',1), ...
    struct('Lx',63,'Ly',63,'afr',0.25,'iflag',1), ...
    struct('Lx',80,'Ly',80,'afr',0.25,'iflag',1), ...
    struct('Lx',43,'Ly',43,'afr',0.25,'iflag',1), ...
    struct('Lx',62,'Ly',62,'afr',0.25,'iflag',1), ...
    struct('Lx',79,'Ly',79,'afr',0.25,'iflag',1), ...
    struct('Lx',41,'Ly',41,'afr',0.25,'iflag',1), ...
    struct('Lx',61,'Ly',61,'afr',0.25,'iflag',1), ...
    struct('Lx',80,'Ly',80,'afr',0.25,'iflag',1), ...
    struct('Lx',42,'Ly',42,'afr',0.25,'iflag',1), ...
    struct('Lx',63,'Ly',63,'afr',0.25,'iflag',1), ...
    struct('Lx',82,'Ly',82,'afr',0.25,'iflag',1), ...
    struct('Lx',41,'Ly',41,'afr',0.25,'iflag',1), ...
    struct('Lx',61,'Ly',61,'afr',0.25,'iflag',1), ...
    struct('Lx',84,'Ly',84,'afr',0.25,'iflag',1), ...
    struct('Lx',43,'Ly',43,'afr',0.25,'iflag',1), ...
    struct('Lx',62,'Ly',62,'afr',0.25,'iflag',1), ...
    struct('Lx',80,'Ly',80,'afr',0.25,'iflag',1), ...
    struct('Lx',41,'Ly',41,'afr',0.25,'iflag',1), ...
    struct('Lx',62,'Ly',62,'afr',0.25,'iflag',1), ...
    struct('Lx',83,'Ly',83,'afr',0.25,'iflag',1)  
];


max_inum  = 20000;   % 最大时间步
threshold = 0.8;     % 填充判定阈值
SEEDS     = 42:142;   % 预设的一组随机数种子（100个）
MAX_SUCC  = 3;       % 每组随机边界最多记录 3 个“填充完成”的样本

%% ====================================

fprintf('边长\tafr\tiflag\tseed\t特征时间(步)\t状态\n');
fprintf('----------------------------------------\n');

for k = 1:numel(PARAMS)
    p = PARAMS(k);

    if p.iflag <= 3
        % ===== 1/2/3：无随机性，只跑一个固定种子 =====
        seed = SEEDS(1);  % 比如 42
        [ft, filled] = compute_feature_time_one(p.Lx, p.Ly, p.afr, p.iflag, max_inum, threshold, seed);

        if filled
            fprintf('%d\t%.3f\t%d\t%d\t%d\t填充完成\n', p.Lx, p.afr, p.iflag, seed, ft);
        else
            fprintf('%d\t%.3f\t%d\t%d\t>%d\t未完全填充\n', p.Lx, p.afr, p.iflag, seed, max_inum);
        end

    else
        % ===== 4/5/6：有随机性的边界，多随机种子，直到有 MAX_SUCC 个填充完成 =====
        success_count = 0;

        for seed = SEEDS
            [ft, filled] = compute_feature_time_one(p.Lx, p.Ly, p.afr, p.iflag, max_inum, threshold, seed);

            if filled
                success_count = success_count + 1;
                fprintf('%d\t%.3f\t%d\t%d\t%d\t填充完成\n', p.Lx, p.afr, p.iflag, seed, ft);
            else
                fprintf('%d\t%.3f\t%d\t%d\t>%d\t未完全填充\n', p.Lx, p.afr, p.iflag, seed, max_inum);
            end

            if success_count >= MAX_SUCC
                break;  % 已经拿到足够多的“填充完成”样本
            end
        end

        if success_count == 0
            % 所有种子都没填满，可以在结果里特别标一下（这里只是提示，你已经有上面那几行打印）
            fprintf('※ L=%d, iflag=%d: 所有种子在 %d 步内均未完全填充。\n', p.Lx, p.iflag, max_inum);
        end
    end
end

fprintf('----------------------------------------\n');
fprintf('计算完成。\n');

end


%% ===== 计算单一参数组的特征时间（更快的一次性模拟版本） =====
function [feature_time, filled_status] = compute_feature_time_one(Lx, Ly, afr, iflag, max_inum, threshold, seed)
    % 返回：
    %   feature_time  : 最早满足“全场 >= threshold”的时间步（若没填满则 = max_inum）
    %   filled_status : 是否在 max_inum 以内填满（true/false）

    % 固定随机种子，确保 iflag=4/5/6 的随机边界可复现
    if nargin < 7
        seed = 42;  % 默认种子，防止其他地方调用时忘了传
    end
    rng(seed);

    % 网格尺寸
    dx = 1;
    dy = 1;

    % 单元数量 m*n
    m = Lx / dx;
    n = Ly / dy;

    % 边界数量
    snum1 = (m + 1) * n; % 水平边数
    snum2 = (n + 1) * m; % 竖直边数
    snum  = snum1 + snum2;

    % 差分计算单元的东南西北线编号（和 core_0 保持一致）
    e_s = 1:m*n;                 % 南
    e_n = (n + 1):(m + 1)*n;     % 北

    t_  = (m + 1)*n + (1:m*(n+1));
    t_(n+1:n+1:end) = [];
    e_w = t_;                    % 西

    t_  = (m + 1)*n + (2:m*(n+1));
    t_(n+1:n+1:end) = [];
    e_e = t_;                    % 东

    % 初始化场
    s0 = zeros(1, snum);   % 边上的值
    e0 = zeros(1, m*n);    % 单元内部值

    % 预设返回值
    feature_time  = max_inum;
    filled_status = false;

    %% ===== 按 core_0 逻辑预生成 4/5/6 的边界模式 =====
    % 单元编号 j: 1..m*n, 行优先：
    %   行 r = 1..m, 每行 n 个单元：索引范围 [(r-1)*n+1, r*n]
    idx_bottom = 1:n;                       % 下边界（r = 1）
    idx_top    = (m-1)*n+1 : m*n;           % 上边界（r = m）
    idx_left   = 1:n:m*n;                   % 左边界（c = 1）
    idx_right  = n:n:m*n;                   % 右边界（c = n）

    % 4/6：四边都算边界
    boundary_mask_all4 = false(1, m*n);
    boundary_mask_all4(idx_bottom) = true;
    boundary_mask_all4(idx_top)    = true;
    boundary_mask_all4(idx_left)   = true;
    boundary_mask_all4(idx_right)  = true;

    % 2/5：三边边界（上 + 左 + 右），下边自由
    boundary_mask_3side = false(1, m*n);
    boundary_mask_3side(idx_top)   = true;
    boundary_mask_3side(idx_left)  = true;
    boundary_mask_3side(idx_right) = true;

    boundary_mask = false(1, m*n);
    boundary_vals = zeros(1, m*n);

    switch iflag
        case 4
            % flag4：四边随机轻微起伏 0.9–1.1
            boundary_mask = boundary_mask_all4;
            num_b = sum(boundary_mask);
            rand_vec = 0.9 + 0.2 * rand(1, num_b);  % U(0.9, 1.1)
            boundary_vals(boundary_mask) = rand_vec;

        case 5
            % flag5：三边边界梯度版，所有 Dirichlet 边界值 >= 0.8
            boundary_mask = boundary_mask_3side;   % 上 + 左 + 右 三条边是 Dirichlet

            Vmin = 0.8;
            Vmax = 1.2;

            % 上边整条 = Vmax
            boundary_vals(idx_top) = Vmax;

            % 左边：从下到上 Vmin -> Vmax
            if numel(idx_left) > 1
                boundary_vals(idx_left) = linspace(Vmin, Vmax, numel(idx_left));
            else
                boundary_vals(idx_left) = (Vmin + Vmax) / 2;
            end

            % 右边：从下到上 Vmin -> Vmax
            if numel(idx_right) > 1
                boundary_vals(idx_right) = linspace(Vmin, Vmax, numel(idx_right));
            else
                boundary_vals(idx_right) = (Vmin + Vmax) / 2;
            end

        case 6
            % flag6：四边“平滑随机团块”
            boundary_mask = false(1, m*n);   % 先全部当“非边界”
            boundary_vals = zeros(1, m*n);   % 默认 0

            % 下边
            vals_b = make_smooth_patches(numel(idx_bottom));
            boundary_vals(idx_bottom) = vals_b;
            boundary_mask(idx_bottom(vals_b > 0)) = true;

            % 上边
            vals_t = make_smooth_patches(numel(idx_top));
            boundary_vals(idx_top) = vals_t;
            boundary_mask(idx_top(vals_t > 0)) = true;

            % 左边
            vals_l = make_smooth_patches(numel(idx_left));
            boundary_vals(idx_left) = vals_l;
            boundary_mask(idx_left(vals_l > 0)) = true;

            % 右边
            vals_r = make_smooth_patches(numel(idx_right));
            boundary_vals(idx_right) = vals_r;
            boundary_mask(idx_right(vals_r > 0)) = true;

        otherwise
            % iflag = 1 / 2 / 3：边界逻辑由下面的 switch 直接控制，不在这里预设
    end

    %% ===== 一次性时间推进，从 1 跑到 max_inum =====
    for step = 1:max_inum
        % 1) 更新单元值 e0
        for j = 1:m*n
            switch iflag
                case 1 % 四边 Dirichlet=1
                    if (j <= n) || (j > (m - 1) * n) || any(j == 1:n:m*n) || any(j == n:n:m*n)
                        e0(j) = 1;
                    else
                        e2 = ((s0(e_e(j)) - 2 * e0(j) + s0(e_w(j))) / dx^2 + ...
                              (s0(e_n(j)) - 2 * e0(j) + s0(e_s(j))) / dy^2);
                        e0(j) = afr * e2 / dx / dy + e0(j);
                    end

                case 2 % 三边 Dirichlet=1
                    if (j > (m - 1) * n) || any(j == 1:n:m*n) || any(j == n:n:m*n)
                        e0(j) = 1;
                    else
                        e2 = ((s0(e_e(j)) - 2 * e0(j) + s0(e_w(j))) / dx^2 + ...
                              (s0(e_n(j)) - 2 * e0(j) + s0(e_s(j))) / dy^2);
                        e0(j) = afr * e2 / dx / dy + e0(j);
                    end

                case 3 % 四边间隔 Dirichlet=1
                    if any(j == (8:8:72)) || ...
                       any(j == 8 + ((m-1)*n:8:(m*n)-9)) || ...
                       any(j == 1 + (n*8:n*8:(m-1)*n)) || ...
                       any(j == (n*8:n*8:(m-8)*n))
                        e0(j) = 1;
                    else
                        e2 = ((s0(e_e(j)) - 2 * e0(j) + s0(e_w(j))) / dx^2 + ...
                              (s0(e_n(j)) - 2 * e0(j) + s0(e_s(j))) / dy^2);
                        e0(j) = afr * e2 / dx / dy + e0(j);
                    end

                case {4,5,6} % 随机/梯度/团块边界
                    if boundary_mask(j)
                        e0(j) = boundary_vals(j);
                    else
                        e2 = ((s0(e_e(j)) - 2 * e0(j) + s0(e_w(j))) / dx^2 + ...
                              (s0(e_n(j)) - 2 * e0(j) + s0(e_s(j))) / dy^2);
                        e0(j) = afr * e2 / dx / dy + e0(j);
                    end

                otherwise
                    % 兜底：当纯扩散处理
                    e2 = ((s0(e_e(j)) - 2 * e0(j) + s0(e_w(j))) / dx^2 + ...
                          (s0(e_n(j)) - 2 * e0(j) + s0(e_s(j))) / dy^2);
                    e0(j) = afr * e2 / dx / dy + e0(j);
            end
        end

        % 2) 单元值 -> 水平交接线 s0(1:snum1)
        for j = 1:snum1
            if j <= n            % 第一行
                s0(j) = e0(j);
           	elseif j > snum1 - n % 最后一行
                s0(j) = e0(j - n);
            else                 % 中间行
                s0(j) = 0.5 * (e0(j) + e0(j - n));
            end
        end

        % 3) 单元值 -> 竖直交接线 s0(snum1+1:end)
        for j = 1:snum2
            if mod(j, n+1) == 1
                s0(snum1 + j) = e0(j + 1 - floor(j / (n+1)) + 1);
            elseif mod(j, n+1) == 0
                s0(snum1 + j) = e0(j - floor(j / (n+1)) - 1);
            else
                id = floor(j / (n+1));
                s0(snum1 + j) = 0.5 * (e0(j-1 - id) + e0(j - id));
            end
        end

        % 4) 检查是否“全场 >= threshold”
        if all(e0 >= threshold)
            feature_time  = step;
            filled_status = true;
            break;
        end
    end
end


function vals = make_smooth_patches(len)
% 生成一条边上的“平滑随机团块”分布：
% - 团块之间间隔大致固定，只做轻微抖动；
% - 团块内部的值在 0.9–1.1 轻微波动；
% - 风格上类似 flag4，只是有“块”和“空”。

    vals = zeros(1, len);
    if len <= 3
        % 太短的话就整体给一个接近常数的随机值
        vals(:) = 0.95 + 0.1 * rand(1, len);
        return;
    end

    % 基准间隔 & 基准半长度（控制团块大小与覆盖率）
    base_period = max(8, round(len/6));          % 基本间隔
    base_half   = max(1, round(base_period/10)); % 团块半长度（覆盖率 ~ 20%）

    pos = round(base_period/2);  % 从靠前的位置开始铺团块

    while pos <= len
        % 轻微抖动中心位置（±2）
        jitter_center = randi([-2, 2]);
        c = pos + jitter_center;
        c = max(1, min(len, c));

        % 轻微抖动半长度（±1）
        jitter_half = randi([-1, 1]);
        h = base_half + jitter_half;
        h = max(1, min(h, floor(len/4)));

        s = max(1, c - h);
        e = min(len, c + h);

        % 团块内部值在 0.9 ~ 1.1 之间
        amp = 0.9 + 0.2 * rand();
        vals(s:e) = amp;

        % 下一块的位置，间隔也轻微抖动（±2）
        pos = pos + base_period + randi([-2, 2]);
    end
end