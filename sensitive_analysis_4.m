function sensitive_analysis_4()
% 敏感性分析（两遍计算版）：
% 第一遍：扫描 Lx/Ly、afr、threshold 组合，记录 50% 与完全铺满步长
% 第二遍：在 50% 步、T50-20、T50+20、T50-2%full、T50+2%full 保存阈值等值线点数据

    L_values   = [40, 60, 80];
    thresholds = [0.65, 0.8, 0.95];
    afr_values = [0.1, 0.25, 0.4];
    iflag      = 1;
    max_steps  = 60000;

    output_dir = fullfile(pwd, 'sensitive_analysis');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fprintf('Lx\tLy\tafr\tthreshold\t50%%步\t完全步\tT50-20\tT50+20\tT50-2%%full\tT50+2%%full\n');
    fprintf('--------------------------------------------------------------------------------------------\n');

    summary_data = struct('Lx', {}, 'Ly', {}, 'afr', {}, 'threshold', {}, ...
        'fifty_step', {}, 'full_step', {}, 'step_minus20', {}, 'step_plus20', {}, ...
        'step_minus_2pct_full', {}, 'step_plus_2pct_full', {});

    for Lx = L_values
        Ly = Lx;
        for afr = afr_values
            for threshold = thresholds

                % ===== 第一遍：只计算 T50 和 Tfull =====
                [fifty_step, full_step] = run_single_case_first_pass( ...
                    Lx, Ly, afr, iflag, threshold, max_steps);

                % ===== 根据 T50 和 Tfull 计算 5 个目标步长 =====
                if isnan(fifty_step)
                    step_minus20        = NaN;
                    step_plus20         = NaN;
                    step_minus_2pct     = NaN;
                    step_plus_2pct      = NaN;
                else
                    % T50-20 / T50+20
                    step_minus20 = fifty_step - 20;
                    if step_minus20 < 1 || step_minus20 > max_steps
                        step_minus20 = NaN;
                    end

                    step_plus20 = fifty_step + 20;
                    if step_plus20 < 1 || step_plus20 > max_steps
                        step_plus20 = NaN;
                    end

                    % T50 ± 2%·Tfull
                    if isnan(full_step)
                        step_minus_2pct = NaN;
                        step_plus_2pct  = NaN;
                    else
                        delta_2pct = floor(0.02 * full_step);

                        step_minus_2pct = fifty_step - delta_2pct;
                        if step_minus_2pct < 1 || step_minus_2pct > max_steps
                            step_minus_2pct = NaN;
                        end

                        step_plus_2pct = fifty_step + delta_2pct;
                        if step_plus_2pct < 1 || step_plus_2pct > max_steps
                            step_plus_2pct = NaN;
                        end
                    end
                end

                % ===== 打印一行结果 =====
                fifty_str    = value_or_na(fifty_step);
                full_str     = value_or_na(full_step);
                minus20_str  = value_or_na(step_minus20);
                plus20_str   = value_or_na(step_plus20);
                minus2pct_str= value_or_na(step_minus_2pct);
                plus2pct_str = value_or_na(step_plus_2pct);

                fprintf('%d\t%d\t%.2f\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
                    Lx, Ly, afr, threshold, fifty_str, full_str, ...
                    minus20_str, plus20_str, minus2pct_str, plus2pct_str);

                % ===== 记录到 summary_data =====
                summary_data(end+1) = struct( ... %#ok<AGROW>
                    'Lx', Lx, ...
                    'Ly', Ly, ...
                    'afr', afr, ...
                    'threshold', threshold, ...
                    'fifty_step', fifty_step, ...
                    'full_step', full_step, ...
                    'step_minus20', step_minus20, ...
                    'step_plus20', step_plus20, ...
                    'step_minus_2pct_full', step_minus_2pct, ...
                    'step_plus_2pct_full', step_plus_2pct);

                % ===== 第二遍：在 5 个目标步长保存等值线 csv =====
                target_steps = [fifty_step, step_minus20, step_plus20, ...
                                step_minus_2pct, step_plus_2pct];
                run_single_case_second_pass( ...
                    Lx, Ly, afr, iflag, threshold, max_steps, ...
                    target_steps, output_dir);
            end
        end
    end

    fprintf('----------------------------------------------\n');
    fprintf('敏感性分析完成。\n');

    if ~isempty(summary_data)
        write_summary_csv(output_dir, summary_data);
    end
end

%% ====================== 第一遍：只算 T50 和 Tfull ======================
function [fifty_step, full_step] = run_single_case_first_pass(Lx, Ly, afr, iflag, threshold, max_steps)
    dx = 1;
    dy = 1;
    m  = Lx / dx;
    n  = Ly / dy;

    snum1 = (m + 1) * n;
    snum2 = (n + 1) * m;

    % 单元与边对应
    e_s = 1:m*n;
    e_n = (n + 1):(m + 1) * n;

    t_ = (m + 1) * n + (1:m*(n + 1));
    t_(n + 1:n + 1:end) = [];
    e_w = t_;

    t_ = (m + 1) * n + (2:m*(n + 1));
    t_(n + 1:n + 1:end) = [];
    e_e = t_;

    s0 = zeros(1, snum1 + snum2);
    e0 = zeros(1, m * n);

    fifty_step = NaN;
    full_step  = NaN;

    for step = 1:max_steps
        % --- 更新单元量 e0 ---
        for j = 1:m*n
            switch iflag
                case 1
                    % 边界单元保持为 1
                    if (j <= n) || (j > (m - 1) * n) || any(j == 1:n:m*n) || any(j == n:n:m*n)
                        e0(j) = 1;
                    else
                        e2 = ((s0(e_e(j)) - 2 * e0(j) + s0(e_w(j))) / dx^2 + ...
                              (s0(e_n(j)) - 2 * e0(j) + s0(e_s(j))) / dy^2);
                        e0(j) = afr * e2 / dx / dy + e0(j);
                    end
                otherwise
                    error('当前敏感性分析只支持 iflag = 1。');
            end
        end

        % --- 更新边量 s0（水平方向） ---
        for j = 1:snum1
            if j <= n
                % 下边界
                s0(j) = e0(j);
            elseif j > snum1 - n
                % 上边界
                s0(j) = e0(j - n);
            else
                % 内部水平边
                s0(j) = 0.5 * (e0(j) + e0(j - n));
            end
        end

        % --- 更新边量 s0（竖直方向） ---
        for j = 1:snum2
            if mod(j, n + 1) == 1
                % 左边界
                s0(snum1 + j) = e0(j + 1 - floor(j / (n + 1)) + 1);
            elseif mod(j, n + 1) == 0
                % 右边界
                s0(snum1 + j) = e0(j - floor(j / (n + 1)) - 1);
            else
                % 内部竖直边
                id = floor(j / (n + 1));
                s0(snum1 + j) = 0.5 * (e0(j - 1 - id) + e0(j - id));
            end
        end

        % --- 计算铺满比例 ---
        filled_ratio = sum(e0 >= threshold) / numel(e0);

        if isnan(fifty_step) && filled_ratio >= 0.5
            fifty_step = step;
        end

        if isnan(full_step) && all(e0 >= threshold)
            full_step = step;
        end

        % T50 和 Tfull 都算出来了就可以停
        if ~isnan(fifty_step) && ~isnan(full_step)
            break;
        end
    end
end

%% ====================== 第二遍：固定步长取状态并保存 csv ======================
function run_single_case_second_pass(Lx, Ly, afr, iflag, threshold, max_steps, target_steps, output_dir)
    dx = 1;
    dy = 1;
    m  = Lx / dx;
    n  = Ly / dy;

    snum1 = (m + 1) * n;
    snum2 = (n + 1) * m;

    % 单元与边对应
    e_s = 1:m*n;
    e_n = (n + 1):(m + 1) * n;

    t_ = (m + 1) * n + (1:m*(n + 1));
    t_(n + 1:n + 1:end) = [];
    e_w = t_;

    t_ = (m + 1) * n + (2:m*(n + 1));
    t_(n + 1:n + 1:end) = [];
    e_e = t_;

    s0 = zeros(1, snum1 + snum2);
    e0 = zeros(1, m * n);

    % 只保留合法的目标步长
    valid_steps = target_steps(~isnan(target_steps) & ...
                               target_steps >= 1 & ...
                               target_steps <= max_steps);
    valid_steps = unique(valid_steps);  % 去重

    if isempty(valid_steps)
        return;
    end

    max_target_step = max(valid_steps);

    for step = 1:max_steps
        % --- 更新单元量 e0 ---
        for j = 1:m*n
            switch iflag
                case 1
                    if (j <= n) || (j > (m - 1) * n) || any(j == 1:n:m*n) || any(j == n:n:m*n)
                        e0(j) = 1;
                    else
                        e2 = ((s0(e_e(j)) - 2 * e0(j) + s0(e_w(j))) / dx^2 + ...
                              (s0(e_n(j)) - 2 * e0(j) + s0(e_s(j))) / dy^2);
                        e0(j) = afr * e2 / dx / dy + e0(j);
                    end
                otherwise
                    error('当前敏感性分析只支持 iflag = 1。');
            end
        end

        % --- 更新边量 s0（水平方向） ---
        for j = 1:snum1
            if j <= n
                s0(j) = e0(j);
            elseif j > snum1 - n
                s0(j) = e0(j - n);
            else
                s0(j) = 0.5 * (e0(j) + e0(j - n));
            end
        end

        % --- 更新边量 s0（竖直方向） ---
        for j = 1:snum2
            if mod(j, n + 1) == 1
                s0(snum1 + j) = e0(j + 1 - floor(j / (n + 1)) + 1);
            elseif mod(j, n + 1) == 0
                s0(snum1 + j) = e0(j - floor(j / (n + 1)) - 1);
            else
                id = floor(j / (n + 1));
                s0(snum1 + j) = 0.5 * (e0(j - 1 - id) + e0(j - id));
            end
        end

        % --- 在目标步长上保存 csv ---
        if any(step == valid_steps)
            save_contour_state(output_dir, Lx, Ly, afr, iflag, ...
                step, threshold, e0, m, n);
        end

        % 所有目标步长都已经到达了就可以停
        if step >= max_target_step
            break;
        end
    end
end

%% ====================== 保存等值线 csv（带阈值到文件名） ======================
function save_contour_state(output_dir, Lx, Ly, afr, iflag, step, threshold, state, m, n)
    if isempty(state) || isnan(step)
        return;
    end

    z = reshape(state, [n, m]);
    curve_data = sample_interface_curve(z, Lx, Ly, threshold);
    if isempty(curve_data)
        return;
    end

    afr_str = trim_numeric_str(afr);
    thr_str = trim_numeric_str(threshold);

    filename = sprintf('Lx%d_Ly%d_afr%s_thr%s_flag%d_timestep%d.csv', ...
        Lx, Ly, afr_str, thr_str, iflag, step);
    file_path = fullfile(output_dir, filename);

    T = array2table(curve_data, 'VariableNames', {'x', 'y', 'angle', 'curvature'});
    writetable(T, file_path);
end

%% ====================== 等值线采样为 0–90° 曲线 ======================
function data = sample_interface_curve(z, Lx, Ly, threshold)
    [n, m] = size(z);
    x_grid = linspace(0, Lx, m + 1);
    y_grid = linspace(0, Ly, n + 1);

    [X, Y] = meshgrid(x_grid(1:m), y_grid(1:n));
    [Xi, Yi] = meshgrid(linspace(0, Lx, m * 10), linspace(0, Ly, n * 10));
    Zi = interp2(X, Y, z, Xi, Yi, 'cubic');

    C = contourc(Xi(1, :), Yi(:, 1), Zi, [threshold threshold]);
    if isempty(C) || size(C, 2) <= 1
        data = [];
        return;
    end

    idx = 1;
    contour_points = [];
    while idx < size(C, 2)
        level = C(1, idx);
        count = C(2, idx);
        if level == threshold && count > 0 && (idx + count) <= size(C, 2)
            segment = C(:, idx + 1:idx + count);
            contour_points = [contour_points; segment.']; %#ok<AGROW>
        end
        idx = idx + count + 1;
    end

    if isempty(contour_points)
        data = [];
        return;
    end

    center_x = Lx / 2;
    center_y = Ly / 2;
    right_top_mask   = contour_points(:, 1) >= center_x & contour_points(:, 2) >= center_y;
    right_top_points = contour_points(right_top_mask, :);

    if isempty(right_top_points)
        data = [];
        return;
    end

    angles = atan2d(right_top_points(:, 2) - center_y, ...
                    right_top_points(:, 1) - center_x);
    angles = mod(angles, 360);

    valid_mask   = angles >= 0 & angles <= 90;
    valid_points = right_top_points(valid_mask, :);
    valid_angles = angles(valid_mask);

    if isempty(valid_points)
        data = [];
        return;
    end

    [sorted_angles, order] = sort(valid_angles);
    sorted_points = valid_points(order, :);

    target_angles  = 0:1:90;
    sampled_points = nan(numel(target_angles), 2);

    for k = 1:numel(target_angles)
        target_angle = target_angles(k);
        [min_diff, closest_idx] = min(abs(sorted_angles - target_angle));

        if min_diff < 0.1
            sampled_points(k, :) = sorted_points(closest_idx, :);
            continue;
        end

        if sorted_angles(closest_idx) > target_angle && closest_idx > 1
            idx1 = closest_idx - 1;
            idx2 = closest_idx;
        elseif closest_idx < numel(sorted_angles)
            idx1 = closest_idx;
            idx2 = closest_idx + 1;
        else
            sampled_points(k, :) = sorted_points(closest_idx, :);
            continue;
        end

        angle1 = sorted_angles(idx1);
        angle2 = sorted_angles(idx2);
        point1 = sorted_points(idx1, :);
        point2 = sorted_points(idx2, :);

        if angle2 ~= angle1
            weight = (target_angle - angle1) / (angle2 - angle1);
            sampled_points(k, :) = point1 + weight * (point2 - point1);
        else
            sampled_points(k, :) = point1;
        end
    end

    if any(isnan(sampled_points(:, 1)))
        data = [];
        return;
    end

    curvatures = compute_curvature(sampled_points);
    data = [sampled_points, target_angles.', curvatures];
end

%% ====================== 小工具函数 ======================
function out = value_or_na(val)
    if isnan(val)
        out = 'NA';
    else
        out = sprintf('%d', val);
    end
end

function txt = trim_numeric_str(num)
    txt = sprintf('%.2f', num);
    txt = regexprep(txt, '0+$', '');
    txt = regexprep(txt, '\.$', '');
end

function curvatures = compute_curvature(points)
    num_points = size(points, 1);
    curvatures = zeros(num_points, 1);

    for idx = 1:num_points
        prev_idx = mod(idx - 2, num_points) + 1;
        next_idx = mod(idx, num_points) + 1;

        P1 = points(prev_idx, :);
        P2 = points(idx, :);
        P3 = points(next_idx, :);

        V1 = P1 - P2;
        V2 = P3 - P2;

        L1 = norm(V1);
        L2 = norm(V2);

        if L1 == 0 || L2 == 0
            curvatures(idx) = 0;
            continue;
        end

        V1 = V1 / L1;
        V2 = V2 / L2;

        dx1 = P1(1) - P2(1);
        dy1 = P1(2) - P2(2);
        dx2 = P3(1) - P2(1);
        dy2 = P3(2) - P2(2);

        dx_prime  = (dx2 - dx1) / 2;
        dy_prime  = (dy2 - dy1) / 2;
        dx2prime  = dx2 + dx1;
        dy2prime  = dy2 + dy1;

        denominator = (dx_prime^2 + dy_prime^2)^(1.5);
        if denominator > 0
            curvatures(idx) = abs(dx_prime * dy2prime - dy_prime * dx2prime) / denominator;
        else
            curvatures(idx) = 0;
        end
    end
end

function write_summary_csv(output_dir, summary_data)
    summary_table = struct2table(summary_data);
    summary_table = summary_table(:, {'Lx', 'Ly', 'afr', 'threshold', 'fifty_step', 'full_step', ...
        'step_minus20', 'step_plus20', 'step_minus_2pct_full', 'step_plus_2pct_full'});
    file_path = fullfile(output_dir, 'sensitive_analysis_summary.csv');
    writetable(summary_table, file_path);
end