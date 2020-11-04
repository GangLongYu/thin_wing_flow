function [y, dy] = wing(x)
%翼型函数及其导数
% x  input:  横坐标
% y  output: 对应翼型纵坐标
% dy output: 翼型导数纵坐标

y = 1/8 * x .^ 2;
dy = 1 / 4 * x;
