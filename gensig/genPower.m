function pwmat = genPower(xt, xsize)
%GENPOWER 得到指定信号的功率
% 假定shape = [row col], 其中row为时间轴(samples), col为信号维数轴(xsize)
% 输出shape = [    col], 其中col为信号维数轴(xsize)
    
    % 按每个信号维度分别计算信号功率
    pwmat = bandpower(xt);

    % 防止忽略最后一维
    pwmat = reshape(pwmat, [xsize 1]);

end

