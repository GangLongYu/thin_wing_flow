function [vx, vy] = induced_v(vortex, point_single)
%多个点涡在某一点的诱导速度的相关项，便于构造矩阵
% vortex        input:  点涡的横纵坐标
% point_single  input:  某一点的横纵坐标
% vx            output: 某一点x方向诱导速度的相关项
% vy            output: 某一点y方向诱导速度的相关项

vx = -(point_single(2)-vortex(2, :))./((vortex(1, :)-point_single(1)).^2+...
    (vortex(2, :)-point_single(2)).^2)/2/pi;
vy = (point_single(1)-vortex(1, :))./((vortex(1, :)-point_single(1)).^2+...
    (vortex(2, :)-point_single(2)).^2)/2/pi;


