function retStruct = genSystemInfo(varargin)
% GENSYSTEMINFO 返回激励信号生成系统的基本信息
% 说明: 
% 1. UInit为一向量, 可取1或0: 
%   a) 若为1, 表示该端口同时用于仿真和辨识, 正常初始化;
%   b) 若为0, 表示该端口仅用于仿真, 令K矩阵对应对角线处为1, 其中size(K) = (usize, usize).

    % 输入解析
    parser = inputParser;
    % parameter
    addParameter(parser, 'Step', 0, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'Period', 0, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'YSize', 0, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'USize', 0, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'Freqs', 0, @(i)(isvector(i)));
    addParameter(parser, 'USat', 0, @(i)(ismatrix(i)));
    addParameter(parser, 'USelect', 0, @(i)(isvector(i)));
    addParameter(parser, 'USnr', 0, @(i)(isvector(i)));
    addParameter(parser, 'YSnr', 0, @(i)(isvector(i)));
    % parser
    parse(parser, varargin{:});
    signal_step = parser.Results.Step;  % 采样频率
    signal_period = parser.Results.Period;  % 采样周期
    y_size = parser.Results.YSize;  % y(t)维数
    u_size = parser.Results.USize;  % u(t)维数
    freqs = parser.Results.Freqs;  % v(t)频率组成
    u_sat = parser.Results.USat;    % u(t)上下限
    u_select = parser.Results.USelect;  % u(t)仿真/辨识端口选择
    u_snr = parser.Results.USnr;    % u(t)信噪比
    y_snr = parser.Results.YSnr;    % y(t)信噪比
    % 处理
    u_size_sim = length(u_select);
    u_select = logical(u_select);
    % 检查
    if sum(u_select) ~= u_size, return; end
    if ~isequal(size(u_sat), [u_size 2]), return; end

    % 计算周期
    % period_sample = fix(2*pi*period/signal_step);
    % 计算频率
    period_sample = fix(signal_period/signal_step);
    freqs = frequenciesConverter(freqs, signal_step, period_sample);

    % 激励信号系统结构
    % v维数
    if min(freqs) == 0, q = length(freqs)-1; v_size = 2*q+1;
    else, q = length(freqs); v_size = 2*q; end
    % v初值
    v_init = genIniter(v_size);
    % S参数
    sv = SvMatBuilder(v_size, freqs(end-q+1:end));
    % K参数
    k_sim = KwMatBuilder(u_select);

    % U参数
    % 用于仿真的SimUv
    u_amp = genAmpBuilder(v_size, u_size, u_sat, u_select);
    u_phase = genPhaseBuilder(v_size, u_size, u_select);
    uv_sim = UvMatBuilder(v_size, u_size_sim, u_amp, u_phase, v_init);
    % 用于辨识的Uv
    uv = uv_sim(u_select, :);
    
    % return
    retStruct = struct( ...
        'YSize', y_size, 'USize', u_size, 'VSize', v_size, ...
        'Freqs', freqs, 'Usnr', u_snr, 'Ysnr', y_snr, ...
        'PeriodSamples', period_sample, ...
        'VInit', v_init, 'Sv', sv, 'Uv', uv, ...
        'SimUSize', u_size_sim, 'SimUv', uv_sim, 'SimKw', k_sim);
end

function return_frequencies = frequenciesConverter(continuous_frequencies, signal_step, period_sample)
    
    % 将连续时域频率转换为对应采样时间下的离散频率值
    inverse_convertor = (signal_step*period_sample)/(2*pi);
    return_frequencies = continuous_frequencies .* inverse_convertor;
    % 取整
    return_frequencies = round(return_frequencies);
    % 去掉过大的频率点
    return_frequencies = return_frequencies(return_frequencies <= fix(period_sample/2));
    % 转换回连续时域频率
    return_frequencies = return_frequencies ./ inverse_convertor;

end

function vinit = genIniter(vsize, ~)
%GENINITER 初始化v(t)状态变量

    persistent randSeed
    % rng('shuffle', 'twister');
    rng(142857,'twister')
    randSeed = randi(intmax('uint32'), 'uint32');

    % 固定幅值(均为1), 随机相位初始化
    vinit = zeros(vsize, 1);
    randStreams = RandStream('mt19937ar', 'Seed', randSeed);
    % if mod(vsize, 2) ~= 1, return, end
    % if mod(vsize, 2) ~= 0, return, end
    % 常量初始化
    if mod(vsize, 2) ~= 0  % 带直流分量
        % 直流分量初始化
        vinit(1) = 1;
        
        % 正弦量初始化
        for iter_loc = 2:2:vsize
            randPhi = 2*pi*rand(randStreams)-pi;
            temp_sin = sin(randPhi);
            temp_cos = cos(randPhi);
            
            vinit(iter_loc) = temp_sin;
            vinit(iter_loc+1) = temp_cos;
        end
    else  % 不带直流分量
        for iter_loc = 1:2:vsize
            randPhi = 2*pi*rand(randStreams)-pi;
            temp_sin = sin(randPhi);
            temp_cos = cos(randPhi);
            
            vinit(iter_loc) = temp_sin;
            vinit(iter_loc+1) = temp_cos;
        end
    end

end

function retmat = genPhaseBuilder(vsize, usize, uselect)
%GENPHASEBUILDER 按照给定的规则简单生成一个相位参数矩阵
% 输出size @ (usize q+1) / (usize, q)

    % 准备返回值
    q = fix(vsize / 2);
    hasDC = vsize - 2*q;
    retmat = zeros(usize, q+1);

    % 按照指定角度初始化
    retmat(:, 2:end) = 45.*ones(usize, q);
    % 转换为弧度制
    retmat = deg2rad(retmat);

    % vsize为偶数时处理输出
    if ~hasDC, retmat = retmat(:, 2:end); end

    % 添加仿真输入(非受控输入)
    temp = retmat;
    retmat = zeros(length(uselect), q+1);
    retmat(uselect, :) = temp;
    
end

function retmat = genAmpBuilder(vsize, usize, usat, uselect)
%GENAMPBUILDER 按照USat给出的输出信号上下限确定幅值, 要求必须有直流信号(可以置零)
% 方法:
% 将直流信号幅值定为(max+min)/2, 
% 余下正弦信号(正)幅值组合, 使得总和为(max-min)/2, 然后依据公式还原所需的Amp
% 由于是欠定方程, 且信号幅值不能过低(否则不易辨识), 采取如下算法:
% 1. 全部频率(除常值外)的幅值按顺序排成一个数组:[a1, a2, ..., ai, ..., aq], 设Σai = K
% 2. 为所有幅值ai分配最低值K/(2q): ai = K/(2q);
% 3. 为所有幅值ai分配剩下的一半: 将全部频率分为二半(若为奇数则少者在低端),
%    随机分配剩下的部分给低端, 剩下的给高端.
% 4. 递归运行.
% 注意: 由于幅值的正负可由相位决定, 因此幅值均假定为正
% 输出size @ (usize q+1) / (usize, q)
% 补充: 当vsize为偶数时, 只使用usat中绝对值较大者a, 范围为[-a, a]

    % 准备返回值
    q = fix(vsize / 2);
    hasDC = vsize - 2*q;
    retmat = zeros(usize, q+1);
    
    % 每行第一为下限, 第二为上限
    usat = sort(usat, 2, 'ascend');
    % vsize为偶数时处理usat
    if ~hasDC
        usat = max(abs(usat), [], 2);
        usat = [-usat usat];
    end

    % 准备随机数种子
    rng('shuffle', 'twister')
    seed = randi(intmax('uint32'));
    randStream = RandStream.create('mrg32k3a', ...
        'NumStreams', usize, 'Seed', seed, 'CellOutput', true);

    % 确定常数幅值
    retmat(:, 1) = (usat(:, 2) + usat(:, 1)) ./ 2;

    % 确定谐振幅值
    cap = (usat(:, 2) - usat(:, 1)) ./ 2;
    % 确定谐振幅值最低值
    temp = cap ./ (2*q);
    retmat(:, 2:end) = repmat(temp, [1 q]);
    % 确定谐振幅值剩余值
    % 递归运算
    temp = recDistribute(q, cap/2);
    retmat(:, 2:end) = retmat(:, 2:end) + temp;
    % 递归函数
    function ampmat = recDistribute(qsize, capacity)
        % 递归到底
        if qsize == 1, ampmat = capacity; return; end
        
        % 递归未到底
        % 计算分配的capacity
        lowcap = zeros(usize, 1);
        for iter = 1:usize, lowcap(iter) = capacity(iter)*rand(randStream{iter}); end
        highcap = capacity - lowcap;
        % 计算分配的lenth
        lowlen = fix(qsize/2);
        highlen = qsize - lowlen;
        % 递归运算
        ampmat = [recDistribute(lowlen, lowcap) recDistribute(highlen, highcap)];
    end

    % vsize为偶数时处理输出
    if ~hasDC, retmat = retmat(:, 2:end); end
    
    % 添加仿真输入(非受控输入)
    temp = retmat;
    retmat = zeros(length(uselect), q+1);
    retmat(uselect, :) = temp;
    
end


function retMatrix = SvMatBuilder(VSize, Freqs)
    q = fix(VSize / 2);
    retMatrix = zeros(VSize);

    % 检查变量size
    if ~(isvector(Freqs) && length(Freqs) == q), return; end

    if mod(VSize, 2) ~= 0  % 带直流分量
        % 计算
        for n = 1:q
            retMatrix(2*n, 2*n+1) =  Freqs(n);
            retMatrix(2*n+1, 2*n) = -Freqs(n);
        end
    else  % 不带直流分量
        % 计算
        for n = 1:q
            retMatrix(2*n-1, 2*n) =  Freqs(n);
            retMatrix(2*n, 2*n-1) = -Freqs(n);
        end
    end
end

function retMatrix = UvMatBuilder(VSize, USize, UAmp, UPhase, VInit)
    % 准备
    q = fix(VSize / 2);
    retMatrix = zeros(USize, VSize);  % 准备返回值
    % if mod(VSize, 2) ~= 1, return; end
    % if mod(VSize, 2) ~= 0, return; end

    if mod(VSize, 2) ~= 0  % 带直流分量
        % 检查变量size
        if ~isequal(size(UAmp), [USize q+1]), return; end
        if ~isequal(size(UPhase), [USize q+1]), return; end
        % 计算
        for nu = 1:USize  % 逐行
            AmpBlock = UAmp(nu, 1)/VInit(1);
            retMatrix(nu, 1) = AmpBlock;
            for nv = 2:(VSize+1)/2  % 每行的逐元素
                VInitBlock = VInit(2*(nv-1):2*(nv-1)+1);
                AmpBlock = UAmp(nu, nv)/(norm(VInitBlock)^2);
                PhaseBlock = phaseBuilder(UPhase(nu, nv));
                retMatrix(nu, 2*(nv-1):2*(nv-1)+1) = AmpBlock .* PhaseBlock * VInitBlock;
            end
        end
    else  % 不带直流分量
        % 检查变量size
        if ~isequal(size(UAmp), [USize q]), return; end
        if ~isequal(size(UPhase), [USize q]), return; end
        % 计算
        for nu = 1:USize  % 逐行
            for nv = 1:VSize/2  % 每行的逐元素
                VInitBlock = VInit(2*nv-1:2*nv);
                AmpBlock = UAmp(nu, nv)/(norm(VInitBlock)^2);
                PhaseBlock = phaseBuilder(UPhase(nu, nv));
                retMatrix(nu, 2*nv-1:2*nv) = AmpBlock .* PhaseBlock * VInitBlock;
            end
        end
    end

    % 嵌套函数
    function retMat = phaseBuilder(phi), retMat = [cos(phi) -sin(phi); sin(phi) cos(phi)]; end
end

function retMatrix = KwMatBuilder(USele)
    USele = double(USele);
    USele = abs(USele-1);
    retMatrix = diag(USele);
end

