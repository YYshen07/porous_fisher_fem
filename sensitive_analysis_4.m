function sensitive_analysis_4()
% 敏感性分析：扫描 Lx/Ly、afr、threshold 组合，记录 50% 与完全铺满步长
% 并在 50% 铺满时刻及其 ±20 步保存阈值等值线点数据

    L_values = [40, 60, 80];
    thresholds = [0.65, 0.8, 0.95];
    afr_values = [0.1, 0.25, 0.4];
    iflag = 1;
    max_steps = 60000;

    output_dir = fullfile(pwd, 'sensitive_analysis');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fprintf('Lx\tLy\tafr\tthreshold\t50%%步\t完全步\n');
    fprintf('----------------------------------------------\n');

    for Lx = L_values
        Ly = Lx;
        for afr = afr_values
            for threshold = thresholds
                result = run_single_case(Lx, Ly, afr, iflag, threshold, max_steps);

                fifty_str = value_or_na(result.fifty_step);
                full_str = value_or_na(result.full_step);
                fprintf('%d\t%d\t%.2f\t%.2f\t%s\t%s\n', ...
                    Lx, Ly, afr, threshold, fifty_str, full_str);

                save_contour_state(output_dir, Lx, Ly, afr, iflag, ...
                    result.fifty_step_minus, threshold, result.state_minus, result.m, result.n);
                save_contour_state(output_dir, Lx, Ly, afr, iflag, ...
                    result.fifty_step, threshold, result.state_fifty, result.m, result.n);
                save_contour_state(output_dir, Lx, Ly, afr, iflag, ...
                    result.fifty_step_plus, threshold, result.state_plus, result.m, result.n);
            end
        end
    end

    fprintf('----------------------------------------------\n');
    fprintf('敏感性分析完成。\n');
end

function result = run_single_case(Lx, Ly, afr, iflag, threshold, max_steps)
    dx = 1;
    dy = 1;
    m = Lx / dx;
    n = Ly / dy;

    snum1 = (m + 1) * n;
    snum2 = (n + 1) * m;

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

    buffer_len = 64;
    recent_steps = nan(1, buffer_len);
    recent_states = cell(1, buffer_len);

    fifty_step = [];
    full_step = [];
    state_minus = [];
    state_fifty = [];
    state_plus = [];

    target_minus = [];
    target_plus = [];

    for step = 1:max_steps
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

        for j = 1:snum1
            if j <= n
                s0(j) = e0(j);
            elseif j > snum1 - n
                s0(j) = e0(j - n);
            else
                s0(j) = 0.5 * (e0(j) + e0(j - n));
            end
        end

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

        idx = mod(step - 1, buffer_len) + 1;
        recent_steps(idx) = step;
        recent_states{idx} = e0;

        filled_ratio = sum(e0 >= threshold) / numel(e0);
        if isempty(fifty_step) && filled_ratio >= 0.5
            fifty_step = step;
            state_fifty = e0;
            target_minus = step - 20;
            if target_minus >= 1
                minus_idx = find(recent_steps == target_minus, 1);
                if ~isempty(minus_idx)
                    state_minus = recent_states{minus_idx};
                end
            end
            target_plus = step + 20;
        end

        if isempty(full_step) && all(e0 >= threshold)
            full_step = step;
        end

        if ~isempty(target_plus) && isempty(state_plus) && step >= target_plus && target_plus <= max_steps
            state_plus = e0;
        end

        need_plus = false;
        if ~isempty(target_plus) && target_plus <= max_steps
            need_plus = true;
        end

        if ~isempty(fifty_step) && (~need_plus || ~isempty(state_plus)) && ~isempty(full_step)
            break;
        end
    end

    if isempty(fifty_step)
        fifty_step = NaN;
    end
    if isempty(full_step)
        full_step = NaN;
    end
    if isempty(target_minus) || target_minus < 1
        fifty_step_minus = NaN;
    else
        fifty_step_minus = target_minus;
    end
    if isempty(target_plus) || target_plus > max_steps
        fifty_step_plus = NaN;
    else
        fifty_step_plus = target_plus;
    end

    result = struct( ...
        'fifty_step', fifty_step, ...
        'full_step', full_step, ...
        'fifty_step_minus', fifty_step_minus, ...
        'fifty_step_plus', fifty_step_plus, ...
        'state_minus', state_minus, ...
        'state_fifty', state_fifty, ...
        'state_plus', state_plus, ...
        'm', m, ...
        'n', n ...
    );
end

function save_contour_state(output_dir, Lx, Ly, afr, iflag, step, threshold, state, m, n)
    if isempty(state) || isnan(step)
        return;
    end

    z = reshape(state, [n, m]);
    points = extract_contour_points(z, Lx, Ly, threshold);
    if isempty(points)
        return;
    end

    filename = sprintf('Lx%d_Ly%d_afr%s_flag%d_timestep%d.csv', ...
        Lx, Ly, trim_numeric_str(afr), iflag, step);
    file_path = fullfile(output_dir, filename);

    T = array2table(points, 'VariableNames', {'x', 'y'});
    writetable(T, file_path);
end

function points = extract_contour_points(z, Lx, Ly, threshold)
    [n, m] = size(z);
    x = linspace(0, Lx, m);
    y = linspace(0, Ly, n);

    C = contourc(x, y, z, [threshold threshold]);
    if isempty(C)
        points = [];
        return;
    end

    idx = 1;
    segments = [];
    while idx < size(C, 2)
        level = C(1, idx);
        count = C(2, idx);
        if level == threshold && count > 0
            segment = C(:, idx + 1:idx + count);
            segments = [segments; segment.']; %#ok<AGROW>
        end
        idx = idx + count + 1;
    end

    points = segments;
end

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

