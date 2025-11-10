function [fig, C_08, z, x, y] = core_0(varargin)
% 计算差分结果并绘图（不保存），返回：
% fig : 图句柄
% C_08: u=0.8 等值线的 contour 矩阵
% z   : n×m 的场值矩阵
% x,y : 网格坐标向量

% —— 输入参数（兼容无参默认）
if nargin>=5
    Lx=varargin{1}; Ly=varargin{2}; afr=varargin{3}; inum=varargin{4}; iflag=varargin{5};
else
    Lx=20; Ly=20; afr=0.5; inum=20; iflag=1;
end
if nargin>=6
    threshold=varargin{6};
else
    threshold=0.8;
end

% 固定随机数种子，保证所有含随机性的边界设置可复现
rng(42);

% —— 网格/单元
dx=1; dy=1;
m=Lx/dx; n=Ly/dy;

% —— 线编号
snum1=(m+1)*n;              % 水平边
snum2=(n+1)*m;              % 竖直边
snum=snum1+snum2;           % 总边数

% —— 单元四邻边索引
e_s=1:m*n;                                  % 南（这里没用到也无妨）
e_n=(n+1):(m+1)*n;                          % 北
t_=(m+1)*n+(1:m*(n+1)); t_(n+1:n+1:end)=[]; e_w=t_;  % 西
t_=(m+1)*n+(2:m*(n+1)); t_(n+1:n+1:end)=[]; e_e=t_;  % 东

% —— 场变量
s0=zeros(1,snum);           % 边上的场值
e0=zeros(1,m*n);            % 单元内部场值

%% ===== 为 4/5/6 准备边界索引与预生成分布 =====
% 单元编号 j: 1..m*n, 行优先：
%   行 r = 1..m, 每行 n 个单元：索引范围 [(r-1)*n+1, r*n]
idx_bottom = 1:n;                       % 下边界（r = 1）
idx_top    = (m-1)*n+1 : m*n;           % 上边界（r = m）
idx_left   = 1:n:m*n;                   % 左边界（c = 1）
idx_right  = n:n:m*n;                   % 右边界（c = n）

% 4 和 6：四边都是边界
boundary_mask_all4 = false(1, m*n);
boundary_mask_all4(idx_bottom) = true;
boundary_mask_all4(idx_top)    = true;
boundary_mask_all4(idx_left)   = true;
boundary_mask_all4(idx_right)  = true;

% 2 / 5：三边边界（上 + 左 + 右），下边自由
boundary_mask_3side = false(1, m*n);
boundary_mask_3side(idx_top)   = true;
boundary_mask_3side(idx_left)  = true;
boundary_mask_3side(idx_right) = true;

% 当前 flag 实际使用的掩码和值（只对 4/5/6 有效）
boundary_mask = false(1, m*n);
boundary_vals = zeros(1, m*n);

switch iflag
    case 4
        % flag4：四边随机轻微起伏 0.9–1.1，类似真实边缘不均匀贴附
        boundary_mask = boundary_mask_all4;
        num_b = sum(boundary_mask);
        rand_vec = 0.9 + 0.2 * rand(1, num_b);  % U(0.9, 1.1)
        boundary_vals(boundary_mask) = rand_vec;

    case 5
        % flag5：三边边界的梯度版，但所有 Dirichlet 边界值都 >= 0.8
        % 设计：
        %   top    = Vmax = 1.2
        %   left   : 从下到上 Vmin -> Vmax （0.8 -> 1.2）
        %   right  : 从下到上 Vmin -> Vmax
        %   bottom : 仍然自由演化（和 flag2 一致，只是三边的值不同）
        boundary_mask = boundary_mask_3side;   % 上 + 左 + 右 三条边是 Dirichlet
    
        Vmin = 0.8;
        Vmax = 1.2;
    
        % 上边：整条 = Vmax
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
        
        % flag6：四边随机团块，但只有团块位置是 Dirichlet，
        % 其他边界点仍然走 PDE，类似 flag3 的“随机升级版”
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
        % iflag = 1,2,3 或其他：这里不预设边界，由后面 case 1/2/3 控制
        % boundary_mask 保持全 false
end
%% =====================================================

% —— 时间推进
for i=1:inum
    % 1) 更新单元值
    for j=1:m*n
        switch iflag
            case 1 % 四边边界（原版）
                if (j<=n)||(j>(m-1)*n)||sum(find(j==[1:n:m*n]))||sum(find(j==[n:n:m*n]))
                    e0(j)=1;
                else
                    e2=((s0(e_e(j))-2*e0(j)+s0(e_w(j)))/dx^2 + ...
                        (s0(e_n(j))-2*e0(j)+s0(e_s(j)))/dy^2);
                    e0(j)=afr*e2/dx/dy+e0(j);
                end

            case 2 % 三边边界（原版）
                % 上 + 左 + 右 = 1，下边自由
                if (j>(m-1)*n)||sum(find(j==[1:n:m*n]))||sum(find(j==[n:n:m*n]))
                    e0(j)=1;
                else
                    e2=((s0(e_e(j))-2*e0(j)+s0(e_w(j)))/dx^2 + ...
                        (s0(e_n(j))-2*e0(j)+s0(e_s(j)))/dy^2);
                    e0(j)=afr*e2/dx/dy+e0(j);
                end

            case 3 % 四边间隔边界（原逻辑，保留）
                if sum([find(j==[8:8:72])  find(j==8+[(m-1)*n:8:(m*n)-9])  ...
                        find(j==1+[n*8:n*8:(m-1)*n])  find(j==[n*8:n*8:(m-8)*n])])
                    e0(j)=1;
                else
                    e2=((s0(e_e(j))-2*e0(j)+s0(e_w(j)))/dx^2 + ...
                        (s0(e_n(j))-2*e0(j)+s0(e_s(j)))/dy^2);
                    e0(j)=afr*e2/dx/dy+e0(j);
                end

            case {4,5,6} % 升级版 flag：统一的 PDE 内部更新 + 预生成边界值
                if boundary_mask(j)
                    % 边界单元：使用预生成的随机 / 梯度 / 团块值
                    e0(j) = boundary_vals(j);
                else
                    % 内部单元：统一走同一个扩散 PDE 更新
                    e2=((s0(e_e(j))-2*e0(j)+s0(e_w(j)))/dx^2 + ...
                        (s0(e_n(j))-2*e0(j)+s0(e_s(j)))/dy^2);
                    e0(j)=afr*e2/dx/dy+e0(j);
                end

            otherwise
                % 其他非法 iflag，就当纯扩散处理（不建议用）
                e2=((s0(e_e(j))-2*e0(j)+s0(e_w(j)))/dx^2 + ...
                    (s0(e_n(j))-2*e0(j)+s0(e_s(j)))/dy^2);
                e0(j)=afr*e2/dx/dy+e0(j);
        end
    end

    % 2) 单元 -> 水平线（与原版一致）
    for j=1:snum1
        if j<=n
            s0(j)=e0(j);
        elseif j>snum1-n
            s0(j)=e0(j-n);
        else
            s0(j)=mean([e0(j) e0(j-n)]);
        end
    end

    % 3) 单元 -> 竖直线（与原版一致）
    for j=1:snum2
        if mod(j,n+1)==1
            s0(snum1+j)=e0(j+1-fix(j/(n+1))+1);
        elseif mod(j,n+1)==0
            s0(snum1+j)=e0(j-fix(j/(n+1))-1);
        else
            id=fix(j/(n+1));
            s0(snum1+j)=mean(e0([j-1:j]-id));
        end
    end
end

% —— 坐标与二维矩阵
x=0:Lx/(m-1):Lx;
y=0:Ly/(n-1):Ly;
z=reshape(e0,[n,m]);
minz=min(z,[],'all'); 
maxz=max(z,[],'all');

% —— 绘图（不保存）
isoLevel = threshold;                 % 只画给定阈值的红色等值线
fig = figure('Visible','off');
hold on; box on;
contourf(x,y,z,[minz:0.2:maxz],'edgecolor','k');
[C_08,~] = contour(x,y,z,[isoLevel isoLevel],'color','r','linewidth',4);
colorbar; caxis([0 1]); axis equal; axis off;

% —— 不保存、不关闭；由外层统一保存与关闭
end

%% ========= 本地函数：生成“平滑随机团块”边界（用于 flag6） =========
function vals = make_smooth_patches(len)
% 生成一条边上的随机团块分布：
% - 团块之间间隔大致稳定，每段长度轻微抖动；
% - 团块内部的值在 0.9–1.1 轻微波动；
% - 风格上类似 flag4，只是有“块”和“空”。

    vals = zeros(1, len);
    if len <= 3
        % 太短就直接给一块近似常数
        vals(:) = 0.95 + 0.1 * rand(1, len);
        return;
    end

    % 基准间隔 & 基准半长度（控制团块大小）
    base_period = max(8, round(len/6));       % 基本间隔
    base_half   = max(1, round(base_period/10)); % 团块半长度（覆盖率大约 ~20% 左右）

    % 从大概 base_period/2 位置开始，往前排团块
    pos = round(base_period/2);

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