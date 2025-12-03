function sensitive_analysis_3d_4()
% 3D图像最低点敏感性分析
% 分析不同参数组合下，3D表面最低点（侵蚀最深点）的特征：
%   - 最低点深度（min_z）
%   - 最低点位置（x, y）
%   - 最低点数量
%   - 平均深度、深度标准差
%   - 最低点到中心的距离

%% ===== 1. 参数配置 =====
L_values     = [40, 60, 80];                  % 区域大小
afr_values   = [0.1, 0.25, 0.4];              % afr参数
thresholds   = [0.65, 0.8, 0.95];             % 阈值（虽然不直接用于最低点，但影响计算）
time_steps   = [200, 400, 800, 1600, 3200, 6400];  % 时间步
iflag        = 1;                             % 边界条件标志

output_dir = fullfile(pwd, 'sensitive_analysis_3d');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% ===== 2. 初始化结果数组 =====
results = struct('Lx', {}, 'Ly', {}, 'afr', {}, 'threshold', {}, ...
    'timestep', {}, 'min_z', {}, 'min_z_scaled', {}, ...
    'min_x', {}, 'min_y', {}, 'num_min_points', {}, ...
    'avg_z', {}, 'std_z', {}, 'min_dist_to_center', {}, 'valid_ratio', {});

fprintf('\n========== 3D 最低点敏感性分析（并行版）==========\n\n');

%% ===== 3. 生成所有参数组合 =====
param_combinations = [];
for Lx = L_values
    for afr = afr_values
        for threshold = thresholds
            for inum = time_steps
                param_combinations(end+1, :) = [Lx, afr, threshold, inum]; %#ok<AGROW>
            end
        end
    end
end

total_cases = size(param_combinations, 1);
fprintf('总共 %d 个参数组合\n', total_cases);

%% ===== 4. 启动并行池 =====
try
    pool = gcp('nocreate');
    if isempty(pool)
        fprintf('启动并行池...\n');
        pool = parpool('local');
        fprintf('并行池已启动，使用 %d 个工作进程\n', pool.NumWorkers);
    else
        fprintf('使用现有并行池，%d 个工作进程\n', pool.NumWorkers);
    end
catch
    fprintf('警告：无法启动并行池，将使用串行计算\n');
    pool = [];
end

%% ===== 5. 并行计算所有案例 =====
fprintf('\n开始并行计算...\n\n');

% 预分配结果数组
results_cell = cell(total_cases, 1);

% 使用 parfor 并行计算
parfor idx = 1:total_cases
    Lx = param_combinations(idx, 1);
    Ly = Lx;  % 保持正方形区域
    afr = param_combinations(idx, 2);
    threshold = param_combinations(idx, 3);
    inum = param_combinations(idx, 4);
    h_ = Lx * 10;  % 高度缩放因子
    
    try
        % 调用 core_0 获取 z, x, y
        [~, ~, z, x, y] = core_0(Lx, Ly, afr, inum, iflag, threshold);
        
        % 分析最低点（传入threshold参数）
        stats = analyze_min_points(z, x, y, Lx, Ly, h_, threshold);
        
        % 保存结果到cell
        results_cell{idx} = struct( ...
            'Lx', Lx, ...
            'Ly', Ly, ...
            'afr', afr, ...
            'threshold', threshold, ...
            'timestep', inum, ...
            'min_z', stats.min_z, ...
            'min_z_scaled', stats.min_z_scaled, ...
            'min_x', stats.min_x, ...
            'min_y', stats.min_y, ...
            'num_min_points', stats.num_min_points, ...
            'avg_z', stats.avg_z, ...
            'std_z', stats.std_z, ...
            'min_dist_to_center', stats.min_dist_to_center, ...
            'valid_ratio', stats.valid_ratio);
        
        fprintf('[%3d/%3d] Lx=%d, afr=%.2f, thr=%.2f, t=%d 完成\n', ...
            idx, total_cases, Lx, afr, threshold, inum);
        
    catch ME
        fprintf('[%3d/%3d] Lx=%d, afr=%.2f, thr=%.2f, t=%d 失败: %s\n', ...
            idx, total_cases, Lx, afr, threshold, inum, ME.message);
        results_cell{idx} = [];
    end
end

%% ===== 6. 整理结果 =====
fprintf('\n整理结果...\n');
results = struct('Lx', {}, 'Ly', {}, 'afr', {}, 'threshold', {}, ...
    'timestep', {}, 'min_z', {}, 'min_z_scaled', {}, ...
    'min_x', {}, 'min_y', {}, 'num_min_points', {}, ...
    'avg_z', {}, 'std_z', {}, 'min_dist_to_center', {}, 'valid_ratio', {});

for idx = 1:total_cases
    if ~isempty(results_cell{idx})
        results(end+1) = results_cell{idx}; %#ok<AGROW>
    end
end

%% ===== 7. 保存结果到CSV =====
if ~isempty(results)
    save_results_to_csv(output_dir, results);
    fprintf('\n所有结果已保存到: %s\n', output_dir);
    
    % 生成汇总统计
    generate_summary_stats(output_dir, results);
else
    fprintf('\n警告：没有成功计算的案例。\n');
end

fprintf('\n========== 分析完成 ==========\n');
end


%% ====================== 分析最低点函数 ======================
function stats = analyze_min_points(z, x, y, Lx, Ly, h_, threshold)
% 分析3D表面的最低点特征（只在 z >= threshold 的区域内分析）
%
% 输入：
%   z         : n×m 的场值矩阵
%   x         : 长度为m的x坐标向量
%   y         : 长度为n的y坐标向量
%   Lx        : x方向区域大小
%   Ly        : y方向区域大小
%   h_        : 高度缩放因子
%   threshold : 阈值（只分析 z >= threshold 的区域）
%
% 输出：
%   stats : 包含最低点统计信息的结构体

    [n, m] = size(z);
    
    % 1. 创建掩码：只考虑 z >= threshold 的区域
    valid_mask = z >= threshold;
    
    if ~any(valid_mask(:))
        % 如果没有满足条件的点，返回NaN
        stats = struct();
        stats.min_z = NaN;
        stats.min_z_scaled = NaN;
        stats.min_x = NaN;
        stats.min_y = NaN;
        stats.num_min_points = 0;
        stats.avg_z = mean(z(:));
        stats.std_z = std(z(:));
        stats.min_dist_to_center = NaN;
        return;
    end
    
    % 2. 在有效区域内找到最小值
    z_valid = z(valid_mask);
    min_z = min(z_valid);
    
    % 3. 找到所有最小值点的位置
    [row_idx, col_idx] = find(z == min_z);
    num_min_points = numel(row_idx);
    
    % 3. 将索引转换为实际坐标
    % z是n×m矩阵，对应 Y行 X列
    x_coords = x(col_idx) * 10;  % 放大10倍（与draw_3d_4一致）
    y_coords = y(row_idx) * 10;
    
    % 4. 计算第一个最低点的坐标（如果有多个，取第一个）
    min_x = x_coords(1);
    min_y = y_coords(1);
    
    % 5. 计算缩放后的深度
    min_z_scaled = (min_z - 1) * h_;  % 负值表示深度
    
    % 6. 计算平均值和标准差（只在有效区域内）
    avg_z = mean(z_valid);
    std_z = std(z_valid);
    
    % 7. 计算最低点到区域中心的距离
    center_x = Lx * 10 / 2;
    center_y = Ly * 10 / 2;
    
    % 如果有多个最低点，计算平均距离
    distances = sqrt((x_coords - center_x).^2 + (y_coords - center_y).^2);
    min_dist_to_center = mean(distances);
    
    % 8. 计算有效区域的覆盖率
    valid_ratio = sum(valid_mask(:)) / numel(z);
    
    % 9. 组装输出
    stats = struct();
    stats.min_z = min_z;
    stats.min_z_scaled = min_z_scaled;
    stats.min_x = min_x;
    stats.min_y = min_y;
    stats.num_min_points = num_min_points;
    stats.avg_z = avg_z;
    stats.std_z = std_z;
    stats.min_dist_to_center = min_dist_to_center;
    stats.valid_ratio = valid_ratio;  % 添加有效区域覆盖率
end


%% ====================== 保存结果到CSV ======================
function save_results_to_csv(output_dir, results)
    % 将结果保存为CSV文件
    results_table = struct2table(results);
    
    % 添加额外的计算列
    results_table.depth_ratio = abs(results_table.min_z_scaled) ./ (results_table.Lx * 10);
    results_table.z_range = results_table.avg_z - results_table.min_z;
    
    % 保存详细结果
    detailed_file = fullfile(output_dir, '3d_minpoint_detailed_results.csv');
    writetable(results_table, detailed_file);
    fprintf('  详细结果已保存: %s\n', detailed_file);
end


%% ====================== 生成汇总统计 ======================
function generate_summary_stats(output_dir, results)
    % 按不同维度生成汇总统计
    
    results_table = struct2table(results);
    
    % 1. 按 L 值汇总
    fprintf('\n----- 按区域大小汇总 -----\n');
    L_values = unique(results_table.Lx);
    for L = L_values'
        mask = results_table.Lx == L;
        fprintf('  L=%d: 平均最小深度=%.2f, 标准差=%.2f\n', ...
            L, mean(results_table.min_z_scaled(mask)), std(results_table.min_z_scaled(mask)));
    end
    
    % 2. 按 afr 汇总
    fprintf('\n----- 按 afr 汇总 -----\n');
    afr_values = unique(results_table.afr);
    for afr = afr_values'
        mask = results_table.afr == afr;
        fprintf('  afr=%.2f: 平均最小深度=%.2f, 标准差=%.2f\n', ...
            afr, mean(results_table.min_z_scaled(mask)), std(results_table.min_z_scaled(mask)));
    end
    
    % 3. 按 threshold 汇总
    fprintf('\n----- 按 threshold 汇总 -----\n');
    thr_values = unique(results_table.threshold);
    for thr = thr_values'
        mask = results_table.threshold == thr;
        fprintf('  threshold=%.2f: 平均最小深度=%.2f, 标准差=%.2f\n', ...
            thr, mean(results_table.min_z_scaled(mask)), std(results_table.min_z_scaled(mask)));
    end
    
    % 4. 按时间步汇总
    fprintf('\n----- 按时间步汇总 -----\n');
    time_values = unique(results_table.timestep);
    for t = time_values'
        mask = results_table.timestep == t;
        fprintf('  t=%d: 平均最小深度=%.2f, 标准差=%.2f\n', ...
            t, mean(results_table.min_z_scaled(mask)), std(results_table.min_z_scaled(mask)));
    end
    
    % 5. 保存汇总统计到CSV
    summary_file = fullfile(output_dir, '3d_minpoint_summary.csv');
    
    % 创建汇总表（修复：直接构建二维cell数组）
    summary_rows = cell(0, 5);  % 初始化为0行5列的cell数组
    
    % 按L汇总
    for L = L_values'
        mask = results_table.Lx == L;
        summary_rows(end+1, :) = {'L', L, mean(results_table.min_z_scaled(mask)), ...
            std(results_table.min_z_scaled(mask)), mean(results_table.min_dist_to_center(mask))};
    end
    
    % 按afr汇总
    for afr = afr_values'
        mask = results_table.afr == afr;
        summary_rows(end+1, :) = {'afr', afr, mean(results_table.min_z_scaled(mask)), ...
            std(results_table.min_z_scaled(mask)), mean(results_table.min_dist_to_center(mask))};
    end
    
    % 按threshold汇总
    for thr = thr_values'
        mask = results_table.threshold == thr;
        summary_rows(end+1, :) = {'threshold', thr, mean(results_table.min_z_scaled(mask)), ...
            std(results_table.min_z_scaled(mask)), mean(results_table.min_dist_to_center(mask))};
    end
    
    % 按时间步汇总
    for t = time_values'
        mask = results_table.timestep == t;
        summary_rows(end+1, :) = {'timestep', t, mean(results_table.min_z_scaled(mask)), ...
            std(results_table.min_z_scaled(mask)), mean(results_table.min_dist_to_center(mask))};
    end
    
    % 创建表格
    summary_table = cell2table(summary_rows, ...
        'VariableNames', {'parameter', 'value', 'avg_min_depth', 'std_min_depth', 'avg_dist_to_center'});
    writetable(summary_table, summary_file);
    fprintf('\n  汇总统计已保存: %s\n', summary_file);
end

