function vortex_tail = tail_pos(vortex, vortex_tail, Gamma, Gamma_tail, v_inf, delta_t, m)
%计算新尾涡的位置
% vortex        input: 翼型上点涡的坐标矩阵
% vortex_tail   input: 尾涡的坐标矩阵
% Gamma         input: 翼型上点涡的涡量强度
% Gamma_tail    input: 尾涡的涡量强度
% v_inf         input: 无穷远处来流
% delta_t       input: 每次迭代的时间间隔
% m             input: 当前迭代次数，迭代几次即当前有几个尾涡

n = size(vortex, 2) - 1;
%% 先计算上一次迭代第1个到最后一个尾涡的诱导速度
u1 = zeros(2, m-1);
for k = 1:m-1
    % 上一次迭代已生成尾涡处的诱导速度
    % 翼型点涡的诱导速度
    [vx, vy] = induced_v(vortex, vortex_tail(:, k));
    % 其他尾涡在该尾涡的诱导速度        
    index = 1:m-1;
    index = find(index~=k);
    [vx_tail, vy_tail] = induced_v(vortex_tail(:, index), vortex_tail(:, k));
    u1(:, k) = [vx; vy] * Gamma(m, :)' + [vx_tail; vy_tail] * Gamma_tail(m, index)' + v_inf';
end
%% 再计算第一个尾涡
% 上一次迭代在翼型后缘处的诱导速度
% 翼型点涡的诱导速度
[vx, vy] = induced_v(vortex(:, 1:n), vortex(:, n+1));
% 其他尾涡在该尾涡的诱导速度            
[vx_tail, vy_tail] = induced_v(vortex_tail, vortex(:, n+1));
u2 = [vx; vy] * Gamma(m, 1:n)' + [vx_tail; vy_tail] * Gamma_tail(m, :)' + v_inf';
%% 
% 先求出各自速度，再同时更新位置，否则新更新的位置会对前面的点造成影响   
vortex_tail(:, 2:m) = vortex_tail(:, 1:m-1) + u1 * delta_t;
vortex_tail(:, 1) = vortex(:, n+1) + u2 * delta_t;
