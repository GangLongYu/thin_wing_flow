function b = vector(vortex, point, vortex_tail, Gamma_tail, normal, v_inf)
%返回方程组 Ax+b=0 的 b 向量
% vortex        input:  翼型上点涡的坐标矩阵，2x(n+1)
% point         input:  标识点坐标矩阵，2xn
% vortex_tail   input:  尾涡坐标矩阵，2 x itermax
% Gamma_tail    input:  尾涡涡量强度向量，1 x itermax
% normal        input:  翼型法向量
% v_inf         input:  来流速度
% b             output: Ax + b = 0

n = size(vortex, 2) - 1;
b = zeros(n+1, 1);
for i = 1:n
    [vx_tail, vy_tail] = induced_v(vortex_tail, point(:, i));
    b(i) = (v_inf(1) + vx_tail * Gamma_tail') * normal(1, i) + ...
        (v_inf(2) + vy_tail * Gamma_tail') * normal(2, i);
end
b(n+1) = sum(Gamma_tail);