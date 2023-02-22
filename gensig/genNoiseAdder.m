function xnt = genNoiseAdder(xt, snr, samples, xsize)
%GENNOISEADDER 返回xt和高斯白噪音的叠加
% 假定shape = [row col], 其中row为时间轴(samples), col为信号维数轴(xsize)

    xt = xt.';  % 转置为(xsize samples)
    xnt = adder(xt, snr);  % 加入噪声
    xnt = reshape(xnt, [xsize samples]);  % 防止一维变量导致的维数错误
    xnt = xnt.';  % 还原为(samples xsize)

end

function retmat = adder(mat, snr)
    
    % 若snr为nan则无噪音
    if sum(isnan(snr)), retmat = mat; return; end

    % 准备参数
    length = size(mat, 1);
    samples = size(mat, 2);
    rng('shuffle', 'twister')
    seed = randi(intmax('uint32'));
    randStream = RandStream.create('mrg32k3a', ...
        'NumStreams', length, 'Seed', seed, 'CellOutput', true);

    % 按每个数据维度分别产生噪声, 保证每个数据维各自满足指定的snr, 且互不相关
    retmat = zeros(length, samples);
    for iter = 1:length
        if isnan(snr), retmat(iter, :) = mat(iter, :);
        else, retmat(iter, :) = awgn(mat(iter, :), snr(iter), 'measured', randStream{iter}); end
    end

end


