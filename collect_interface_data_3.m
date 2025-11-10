function [interface_points, curvatures] = collect_interface_data_3(varargin)
% collect_interface_data_5 - 高精度收集扩散界面点，计算曲率
% 在右上1/4区域 (x ≥ Lx/2 且 y ≥ Ly/2) 的0-90度范围内均匀采样91个点
% 总共模拟steps步，每step_interval步收集一次数据

if nargin==5 % 当被调用有输入参数时
    Lx=varargin{1};
    Ly=varargin{2};
    afr=varargin{3};
    total_steps=varargin{4}; % 总时间步数
    step_interval=varargin{5}; % 数据收集间隔
else % 缺省时为以下参数
    Lx=60;
    Ly=60;
    afr=0.25;
    total_steps=3000; % 默认总步数为1500
    step_interval=60; % 默认每10步收集一次
end

% 只使用iflag=1的情况
iflag = 1;

% 计算需要收集数据的时间步
collect_steps = 1:step_interval:total_steps;
num_collect = length(collect_steps);

fprintf('总时间步: %d, 收集间隔: %d, 将收集 %d 个时间点的数据\n', ...
    total_steps, step_interval, num_collect);

%% 网格尺寸
dx=1;
dy=1;

%% 单元数量m*n
m=Lx/dx;
n=Ly/dy;

% 创建正确的网格坐标
x_grid = linspace(0, Lx, m+1);
y_grid = linspace(0, Ly, n+1);

% 检查网格维度
fprintf('网格维度: m=%d (x方向单元数), n=%d (y方向单元数)\n', m, n);

%% 线编号，单元编号
% 边界数量
snum1=(m+1)*n; % 水平边数，左到右，下到上编号
snum2=(n+1)*m; % 竖直边数，左到右，下到上编号
snum=snum1+snum2;

% 差分计算单元的东南西北线编号
e_s=1:m*n; % 南
e_n=(n+1):(m+1)*n;% 北
t_=(m+1)*n+(1:m*(n+1));
t_(n+1:n+1:end)=[];
e_w=t_;% 西
t_=(m+1)*n+(2:m*(n+1));
t_(n+1:n+1:end)=[];
e_e=t_;% 东

%% 用于存储每个时间步的界面点和曲率
interface_points = cell(num_collect, 1);  % 存储每个收集时间步的界面点
curvatures = cell(num_collect, 1);        % 存储每个收集时间步的曲率值

%% 差分计算
s0=zeros(1,snum);
e0=zeros(1,m*n);

% 用于记录当前收集时间步的索引
curr_collect_idx = 0;

% 循环计算所有时间步
for i=1:total_steps
    % 根据偏微分公式差分法，计算下一时刻的值。
    % 根据边界线+上一时刻值=>当前时刻值
    for j=1:m*n        
        % 只保留iflag=1的情况：四边边界
        if (j<=n)||(j>(m-1)*n)||sum(find(j==[1:n:m*n]))||sum(find(j==[n:n:m*n]))
            e0(j)=1;
        else
            % 根据论文公式，只取第二项
            e2=((s0(e_e(j))-2*e0(j)+s0(e_w(j)))/dx^2+(s0(e_n(j))-2*e0(j)+s0(e_s(j)))/dy^2);
            e0(j)=afr*e2/dx/dy+e0(j);
        end
    end
    
    % 根据单元值=>计算相邻水平交接线的值
    for j=1:snum1
        if j<=n % 第一行，相邻只有一个单元，线编号与单元编号相同
            s0(j)=e0(j);
        elseif j>snum1-n % 最后一行，相邻只有一个单元，线编号与单元编号相差n
           s0(j)=e0(j-n);
        else  % 其余行，线编号为相邻单元编号
            s0(j)=mean([e0(j)  e0(j-n)]);
        end
    end
    
    % 根据单元值=>计算相邻竖直交接线的值
    for j=(1:snum2)        
        if mod(j,n+1)==1 % 使用mod函数获取第一列，并获取相邻单元编号
            s0(snum1+j)=e0(j+1-fix(j/(n+1))+1);
        elseif mod(j,n+1)==0 % 使用mod函数获取最后一列，并获取相邻单元编号
            s0(snum1+j)=e0(j-fix(j/(n+1))-1);
        else % 其余列，线编号为相邻单元编号
            id=fix(j/(n+1));
            s0(snum1+j)=mean(e0([j-1:j]-id));
        end
    end
    
    % 检查当前步是否需要收集数据
    is_collect_step = ismember(i, collect_steps);
    
    if is_collect_step
        % 更新收集索引
        curr_collect_idx = curr_collect_idx + 1;
        
        % 将当前时刻的单元值重新排列为二维数组
        z = reshape(e0, [n, m]);
        
        % 使用高分辨率插值进行精确的等值线提取
        % 插值到更密集的网格上（10倍分辨率）
        [X, Y] = meshgrid(x_grid(1:m), y_grid(1:n));
        [Xi, Yi] = meshgrid(linspace(0, Lx, m*10), linspace(0, Ly, n*10));
        Zi = interp2(X, Y, z, Xi, Yi, 'cubic');
        
        % 提取0.8密度等值线
        [C, h] = contour(Xi, Yi, Zi, [0.8 0.8]);
        delete(h); % 删除等值线对象，不显示图形
        
        % 从等值线数据中提取点
        if ~isempty(C) && size(C, 2) > 1
            % 获取等值线点
            level_start = 1;
            contour_points = [];
            
            while level_start < size(C, 2)
                num_points = C(2, level_start);
                segment = C(:, level_start+1:level_start+num_points);
                contour_points = [contour_points, segment];
                level_start = level_start + num_points + 1;
                if level_start >= size(C, 2)
                    break;
                end
            end
            
            % 转置为N×2的点坐标
            if ~isempty(contour_points)
                points = contour_points';
                
                % 只保留右上1/4区域的点 (x ≥ Lx/2 且 y ≥ Ly/2)
                right_top_points = points(points(:,1) >= Lx/2 & points(:,2) >= Ly/2, :);
                
                % 如果有点，按角度排序并均匀提取91个点
                if ~isempty(right_top_points)
                    % 中心点
                    center_x = Lx/2;
                    center_y = Ly/2;
                    
                    % 计算每个点相对于中心的角度
                    angles = atan2d(right_top_points(:,2) - center_y, right_top_points(:,1) - center_x);
                    
                    % 规范化角度到0-90度范围
                    angles = mod(angles, 360);
                    valid_points = right_top_points(angles >= 0 & angles <= 90, :);
                    valid_angles = angles(angles >= 0 & angles <= 90);
                    
                    % 如果有效点足够多，则均匀采样91个点
                    if ~isempty(valid_points)
                        % 按角度排序
                        [sorted_angles, idx] = sort(valid_angles);
                        sorted_points = valid_points(idx, :);
                        
                        % 目标角度：0到90度，步进1度
                        target_angles = 0:1:90;  % 91个点
                        sampled_points = zeros(91, 2);
                        
                        % 对每个目标角度，找到最接近的点或插值
                        for j = 1:length(target_angles)
                            target_angle = target_angles(j);
                            
                            % 找到最接近的两个点
                            [~, closest_idx] = min(abs(sorted_angles - target_angle));
                            
                            if abs(sorted_angles(closest_idx) - target_angle) < 0.1
                                % 如果有非常接近的点，直接使用
                                sampled_points(j, :) = sorted_points(closest_idx, :);
                            else
                                % 找到角度两侧的点进行线性插值
                                if sorted_angles(closest_idx) > target_angle && closest_idx > 1
                                    idx1 = closest_idx - 1;
                                    idx2 = closest_idx;
                                elseif closest_idx < length(sorted_angles)
                                    idx1 = closest_idx;
                                    idx2 = closest_idx + 1;
                                else
                                    % 边界情况处理
                                    sampled_points(j, :) = sorted_points(closest_idx, :);
                                    continue;
                                end
                                
                                % 线性插值
                                angle1 = sorted_angles(idx1);
                                angle2 = sorted_angles(idx2);
                                point1 = sorted_points(idx1, :);
                                point2 = sorted_points(idx2, :);
                                
                                if angle2 ~= angle1
                                    weight = (target_angle - angle1) / (angle2 - angle1);
                                    sampled_points(j, :) = point1 + weight * (point2 - point1);
                                else
                                    sampled_points(j, :) = point1;
                                end
                            end
                        end
                        
                        % 计算曲率
                        curve_curvatures = zeros(91, 1);
                        
                        % 使用有限差分法计算曲率
                        for j = 1:91
                            % 对于边界点，使用循环边界
                            prev_idx = mod(j-2, 91) + 1;
                            next_idx = mod(j, 91) + 1;
                            
                            % 获取三个点坐标
                            P1 = sampled_points(prev_idx, :);
                            P2 = sampled_points(j, :);
                            P3 = sampled_points(next_idx, :);
                            
                            % 计算向量及其模
                            V1 = P1 - P2;
                            V2 = P3 - P2;
                            L1 = norm(V1);
                            L2 = norm(V2);
                            
                            if L1 > 0 && L2 > 0
                                % 归一化向量
                                V1 = V1 / L1;
                                V2 = V2 / L2;
                                
                                % 计算单位切向量（平均方向）
                                T = (V1 + V2) / norm(V1 + V2);
                                
                                % 计算曲率（使用有限差分近似）
                                dx1 = P1(1) - P2(1);
                                dy1 = P1(2) - P2(2);
                                dx2 = P3(1) - P2(1);
                                dy2 = P3(2) - P2(2);
                                
                                % 一阶导数（中心差分）
                                dx_prime = (dx2 - dx1) / 2;
                                dy_prime = (dy2 - dy1) / 2;
                                
                                % 二阶导数
                                dx_double_prime = dx2 + dx1;
                                dy_double_prime = dy2 + dy1;
                                
                                % 计算曲率
                                denominator = (dx_prime^2 + dy_prime^2)^1.5;
                                if denominator > 0
                                    curve_curvatures(j) = abs(dx_prime * dy_double_prime - dy_prime * dx_double_prime) / denominator;
                                else
                                    curve_curvatures(j) = 0;
                                end
                            else
                                curve_curvatures(j) = 0;
                            end
                        end
                        
                        % 存储采样后的91个点和对应的曲率
                        interface_points{curr_collect_idx} = sampled_points;
                        curvatures{curr_collect_idx} = curve_curvatures;
                    else
                        interface_points{curr_collect_idx} = [];
                        curvatures{curr_collect_idx} = [];
                    end
                else
                    interface_points{curr_collect_idx} = [];
                    curvatures{curr_collect_idx} = [];
                end
            else
                interface_points{curr_collect_idx} = [];
                curvatures{curr_collect_idx} = [];
            end
        else
            interface_points{curr_collect_idx} = [];
            curvatures{curr_collect_idx} = [];
        end
        
        % 显示进度
        fprintf('完成时间步 %d/%d (收集点 %d/%d): %d个界面点\n', ...
            i, total_steps, curr_collect_idx, num_collect, ...
            size(interface_points{curr_collect_idx}, 1));
    elseif mod(i, 100) == 0
        % 每100步显示一次进度
        fprintf('计算时间步 %d/%d\n', i, total_steps);
    end
end

end
