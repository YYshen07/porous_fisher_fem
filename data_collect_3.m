function results = data_collect_3()
% run_interface_analysis_5b
% 功能：使用“老版” collect_interface_data_5 的完整逻辑，
%       在若干时间步上提取 u=0.8 界面点（右上 1/4，按 0–90° 采样 91 点），
%       计算曲率，并以 CSV / MAT 保存。
%
% 说明（微小数值差异）：
% 本函数沿用最初用于论文分析的 collect_interface_data_5 流水线，
% 包括在细网格上插值后再提取 0.8 等值线等数值实现细节。
% 与后续基于 core_0 + data_collect 的实现相比，界面提取方式略有不同，
% 因而在特征时间、曲率等量上会出现少量步数级 / 数百分点级的差异，
% 这些差异属于数值实现层面的“微瑕”，不会改变 T_fill ∝ L^2 等主尺度律
% 以及不同孔径之间的相对比较结论。
%
% 调用形式仿照新版 data_collect：
%   - 不带输入参数
%   - 通过 PARAMS 结构数组配置多组参数并批量处理
%
% 目录结构与旧版保持一致：
%   interface_data_Lx%d_Ly%d_afr%.2f_csv/
%       timestep_0001.csv  (x, y, angle, curvature)
%       timestep_0061.csv
%       ...
%       index.csv          (time_step, filename, num_points)
%       parameters.csv     (Lx, Ly, afr, total_steps, step_interval)
%   interface_data_Lx%d_Ly%d_afr%.2f_steps%d_interval%d.mat  (备份)

%% ===== 1. 批量参数配置（按需修改） =====
PARAMS = [ ...
    struct('Lx',60,'Ly',60,'afr',0.25,'total_steps',1801,'step_interval',60) ...
    % 如果需要更多参数组，在这里继续加：
    % , struct('Lx',80,'Ly',80,'afr',0.25,'total_steps',2001,'step_interval',60)
];
%% =====================================

n_case = numel(PARAMS);
results = cell(n_case, 1);

for k = 1:n_case
    p = PARAMS(k);
    fprintf('\n===== 运行第 %d 组参数 =====\n', k);
    results{k} = run_one_case_legacy(p);
end

% 如果只有一组参数，直接返回结构体，兼容原来的用法
if n_case == 1
    results = results{1};
end

fprintf('\n所有参数组的数据收集与保存已完成。\n');

end

function result_data = run_one_case_legacy(p)
    % 从结构体中取出参数
    Lx           = p.Lx;
    Ly           = p.Ly;
    afr          = p.afr;
    total_steps  = p.total_steps;
    step_interval = p.step_interval;

    % 显示运行参数
    fprintf('运行参数: Lx=%d, Ly=%d, afr=%.2f, 总步数=%d, 收集间隔=%d\n', ...
            Lx, Ly, afr, total_steps, step_interval);

    % 调用 collect_interface_data_5 收集数据（老代码逻辑）
    fprintf('开始收集界面数据和计算曲率（collect_interface_data_3）...\n');
    [interface_points, curvatures] = collect_interface_data_3(Lx, Ly, afr, total_steps, step_interval);
    fprintf('数据收集完成。\n');

    % 计算实际收集的时间点数
    num_collect = length(interface_points);

    % 显示每个收集时间点的数据点数
    fprintf('\n每个收集时间点的界面点数:\n');
    for i = 1:num_collect
        time_step = (i - 1) * step_interval + 1;
        if ~isempty(interface_points{i})
            fprintf('收集点 %d (时间步 %d): %d 个界面点\n', ...
                    i, time_step, size(interface_points{i}, 1));
        else
            fprintf('收集点 %d (时间步 %d): 无界面点\n', i, time_step);
        end
    end

    % 将结果保存到结构体
    result_data = struct( ...
        'interface_points', {interface_points}, ...
        'curvatures',       {curvatures}, ...
        'parameters',       struct('Lx', Lx, 'Ly', Ly, 'afr', afr, ...
                                   'total_steps', total_steps, ...
                                   'step_interval', step_interval) );

    % ========= 保存 CSV =========
    fprintf('正在保存 CSV 格式数据...\n');

    % 创建保存目录（与老版本一致）
    csv_dir = sprintf('interface_data_Lx%d_Ly%d_afr%.2f_csv', Lx, Ly, afr);
    if ~exist(csv_dir, 'dir')
        mkdir(csv_dir);
    end

    % 保存每个时间步的界面点到单独的 CSV 文件
    for i = 1:num_collect
        time_step = (i - 1) * step_interval + 1;
        if ~isempty(interface_points{i}) && ~isempty(curvatures{i})
            points = interface_points{i};  % 91×2, [x,y]
            curv   = curvatures{i};        % 91×1

            % 角度：老逻辑中就是固定 0–90°
            angles = (0:1:90)';            % 91×1

            % [x, y, angle, curvature]
            data = [points, angles, curv];

            csv_filename = sprintf('%s/timestep_%04d.csv', csv_dir, time_step);

            T = array2table(data, 'VariableNames', {'x', 'y', 'angle', 'curvature'});
            writetable(T, csv_filename);
        end
    end

    % 创建索引文件 index.csv，方便 Python 批量处理
    index_filename = sprintf('%s/index.csv', csv_dir);
    time_steps = zeros(num_collect, 1);
    filenames  = cell(num_collect, 1);
    num_points = zeros(num_collect, 1);

    for i = 1:num_collect
        time_steps(i) = (i - 1) * step_interval + 1;
        if ~isempty(interface_points{i})
            filenames{i}  = sprintf('timestep_%04d.csv', time_steps(i));
            num_points(i) = size(interface_points{i}, 1);
        else
            filenames{i}  = 'NA';
            num_points(i) = 0;
        end
    end

    index_table = table(time_steps, filenames, num_points, ...
                        'VariableNames', {'time_step', 'filename', 'num_points'});
    writetable(index_table, index_filename);

    % 另外保存参数文件 parameters.csv
    params_table = table([Lx; Ly; afr; total_steps; step_interval], ...
                         'RowNames', {'Lx', 'Ly', 'afr', 'total_steps', 'step_interval'}, ...
                         'VariableNames', {'Value'});
    params_filename = sprintf('%s/parameters.csv', csv_dir);
    writetable(params_table, params_filename, 'WriteRowNames', true);

    fprintf('CSV 数据已保存到目录: %s\n', csv_dir);

    % ========= 保存 MAT 备份 =========
    % 保存到 CSV 目录内部，方便一起管理
    data_filename = sprintf('%s/interface_data_Lx%d_Ly%d_afr%.2f_steps%d_interval%d.mat', ...
        csv_dir, Lx, Ly, afr, total_steps, step_interval);
    save(data_filename, 'interface_points', 'curvatures', 'result_data');
    fprintf('MAT 文件已保存为: %s\n', data_filename);

end