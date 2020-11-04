function normal = normal_vector(point_x)
%求曲线法向量
% point input:  标识点横纵坐标
% nx    output: 法向量横坐标
% ny    output: 法向量纵坐标

ny = ones(size(point_x));
[~, dy] = wing(point_x);
nx = -dy;
% 向量单位化，行向量，加快收敛
n = size(nx, 2);
normal = [nx; ny];
for i = 1:n
    normal(:, i) = normal(:, i) / norm(normal(:, i));
end
