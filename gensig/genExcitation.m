function [t, vt, ut] = genExcitation(systems, step, samples, vinit)
%GENEXCITATION 产生激励信号
% Linear Nonhyperbolic Zero Dynamic ODE & SDE
% ODE:  dv(t) = Sv(t)dt
%        u(t) = Uv(t) + w(t) (?)

    % 使用持久变量加速代码执行
    persistent Sv

    % 函数变量初始化
    Sv = systems.Sv;

    % 时间轴
    span = 0:step:(samples-1)*step;  % 求解长度
    t = span.';  % 返回时间轴

    % 定义嵌套函数
    function ret = vfunc(~, vt), ret = Sv*vt; end  % 返回$$\dot{v} = S \cdot v$$对应的函数向量
    % spanlimit = @(x) min(max(x, 1), samples);  % 时间区间饱和函数

    % 解算$$\dot{v} = S \cdot v$$ (仅一次实验)
    % odeset附注: 'Stats', 'on', 'OutputFcn', 'odephas3'
    odeOption = odeset('reltol',1e-6, 'abstol',1e-6);  % 求解选项
    [~, vt] = ode45(@vfunc, span, vinit, odeOption);  % 解算ODE
    
    % 计算u(t)
    ut = vt * (systems.SimUv.');

end

