%非定常薄翼绕流
% 参数列表如下，(*)表示可设置
% n(*)          翼型上标识点的数目，翼型标示点和点涡取均匀间隔
% vortex        翼型上点涡的横纵坐标，2 x (n+1)
% point         翼型上点涡的横纵坐标，2 x n
% alpha         攻角为0，不可设置
% normal        翼型在标识点处的法向量，由翼型和标识点决定，经归一化处理，2 x n
% v_inf(*)      无穷远出来流速度，1 x 2
% delta_t(*)    每次迭代的时间间隔
% Gamma         翼型上点涡的环量，(itermax+1) x (n+1)
% vortex_tail   尾涡的横纵坐标，2 x itermax
% Gamma_tail    尾涡的环量，(itermax+1) x itermax

% 1. 因为影响因素较多，故设置较大的迭代次数强迫程序收敛，而非设阈值判断是否达到收敛
% 2. 画图过程利用pause以求动态效果，请耐心等待
% 3. 出于美观考虑，figure的参数均已设置，相关代码可根据实际情况在main.m修改
%% 
clc
clear
close all
%% 初始化
% 取点数，设所取点间隔均匀
n = 100; 
% 翼型上取 n+1 个点涡，2 x n
vortex_x = linspace(0, 1, n+1); 
vortex_y = wing(vortex_x);
vortex = [vortex_x; vortex_y];
% 翼型上去取 n 个标识点，2 x (n-1)
point_x = linspace(0+1/n/2, 1-1/n/2, n);
point_y = wing(point_x);
point = [point_x; point_y];
% 攻角
% alpha = 0;
% 无穷远来流，1 x 2
v_inf = [5, 0];
% 法向量
normal = normal_vector(point_x);
% 最大迭代次数
itermax = 500;
% 每次迭代时间间隔
delta_t = 0.001;
% 翼型上点涡的环量
Gamma = zeros(itermax+1, n+1);
% 尾涡
Gamma_tail = zeros(itermax+1, itermax);
% 记录每次迭代的尾涡位置，每次迭代记录2行数据
vortex_tail = zeros(2*(itermax+1), itermax);
%% t=0
% 方程组的矩阵形式
% A每次都是相同的，算一次就行
A = ones(n+1, n+1);
for i = 1:n
    [vx, vy] = induced_v(vortex, point(:, i));   
    A(i, :) = vx * normal(1, i) + vy * normal(2, i);
end
b = vector(vortex, point, vortex_tail(1:2, :), Gamma_tail(1, :), normal, v_inf);
Gamma(1, :) = -(A\b)';
%% 开始产生尾涡 + 可视化处理
for m = 1:itermax  % m 即此次迭代完尾涡个数    
    % 计算新尾涡位置，由上一时刻涡的位置和强度决定
    % 先算位置再算强度是因为强度对位置有影响
    vortex_tail((2*m+1):(2*m+2), :) = tail_pos(vortex, vortex_tail((2*m-1):(2*m), :), Gamma, Gamma_tail, v_inf, delta_t, m);
    % 计算新尾涡涡量强度
    Gamma_tail(m+1, 2:m) = Gamma_tail(m, 1:m-1);
    Gamma_tail(m+1, 1) = Gamma(m, n+1);
    % 可视化，画出尾涡位置，调整figure参数    
    set(gcf, 'Position', get(0, 'ScreenSize')) 
    % 圈的大小表征涡量强度，不分正负
    sz = abs(Gamma_tail(m+1, 1:m))*500;
    scatter(vortex_tail(2*m+1, 1:m), vortex_tail(2*m+2, 1:m), sz, 'b', 'filled')
    axis equal
    xlabel('x(m)'), ylabel('y(m)')
    set(gca, 'FontSize', 12);
    box on
    title('涡位置变化动态图')
    % 给定坐标边界看图比较方便，但是若迭代次数和时间间隔改变，边界需调整
    xlim([0 3.6])
        
    % 求方程组
    b = vector(vortex, point, vortex_tail((2*m+1):(2*m+2), :), Gamma_tail(m+1, :), normal, v_inf);
    Gamma(m+1, :) = -(A\b)';
    
    % 可视化，画出翼型上点涡位置
    hold on 
    % 圈的大小表征涡量强度，不分正负
    sz = abs(Gamma(m+1, :))*500;
    scatter(vortex(1, :), vortex(2, :), sz, 'r', 'filled')
    hold off
    legend('尾涡', '附着涡', 'Location', 'northwest', 'FontSize', 16)
    pause(0.05)
end
pause(0.1)
%% 第一个分离涡轨迹图
figure
% 翼型
set(gcf, 'Position', get(0, 'ScreenSize'))
box on
plot(vortex_x, vortex_y, 'r-.', 'LineWidth', 1)
% 第一个分离涡轨迹
hold on
plot([vortex_x(end), vortex_tail(3, 1)], [vortex_y(end), vortex_tail(4, 1)], ...
    'b-', 'MarkerSize', 2, 'LineWidth', 1)
for i = 2:itermax
    hold on
    plot([vortex_tail(2*i-1, i-1), vortex_tail(2*i+1, i)], ...
        [vortex_tail(2*i, i-1), vortex_tail(2*i+2, i)], 'b-', 'MarkerSize', 2, 'LineWidth', 1)
    hold off
end
% axis equal
xlabel('x(m)'), ylabel('y(m)')
set(gca, 'FontSize', 14);
title('第一个分离涡的轨迹图')
legend('翼型', '第一个分离涡轨迹', 'Location', 'northwest', 'FontSize', 16)
% 给定坐标边界看图比较方便，但是若迭代次数和时间间隔改变，边界需调整
xlim([0 3.6])
saveas(gcf, '第一个分离涡轨迹图.png')
pause(0.1)
%% 后缘涡量强度变化，作图
figure
set(gcf, 'Position', get(0, 'ScreenSize'))
box on
plot(1:itermax+1, Gamma(:, n+1), 'r', 'linewidth', 2);
set(gca, 'FontSize', 12);
xlabel('迭代次数'), ylabel('\Gamma(m^2/s)')
title('后缘\Gamma随迭代次数变化趋势图')
saveas(gcf, '后缘涡量强度变化趋势图.png')

