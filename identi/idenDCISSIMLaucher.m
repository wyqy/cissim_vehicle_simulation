function ret_struct = idenDCISSIMLaucher(un_test, varargin)
% IDENDCISSIMLAUCHER discrete-cISSIM系统辨识 - 初始化

    % 输入解析
    parser = inputParser;
    % 输入定义
    addParameter(parser, 'y_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'u_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'period_samples', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'cutted_periods', 0, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'dcissim_type', 'offline', @(i)(ischar(i)));
    addParameter(parser, 'isim_excitation_type', 'full', @(i)(ischar(i)));
    addParameter(parser, 'x_size_upbound', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'sim_ss_bdx_type', 'analytical', @(i)(ischar(i)));
    addParameter(parser, 'sim_ss_d_type', 'null', @(i)(ischar(i)));
    addParameter(parser, 'cov_cross_type', 'null', @(i)(ischar(i)));
    addParameter(parser, 'sim_x_size_type', 'estimate', @(i)(ischar(i)));
    addParameter(parser, 'cov_order_type', 'estimate', @(i)(ischar(i)));
    addParameter(parser, 'sim_x_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'cov_order', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'cov_est_type', 'simple', @(i)(ischar(i)));
    % 输入提取
    parse(parser, varargin{:});
    y_size = parser.Results.y_size;  % 输出信号维数
    u_size = parser.Results.u_size;  % 输出信号维数
    period_samples = parser.Results.period_samples;  % 单周期采样点数
    cutted_periods = parser.Results.cutted_periods;  % 切除的非稳态部分, 按周期计数 (only in 离线和方差估计)
    dcissim_type = parser.Results.dcissim_type;  % 算法运行方式 - 离线or在线
    isim_excitation_type = parser.Results.isim_excitation_type;  % ISIM激励方式
    x_size_upbound = parser.Results.x_size_upbound;  % 状态变量上界
    sim_ss_bdx_type = parser.Results.sim_ss_bdx_type;  % SIM求解方式
    sim_ss_d_type = parser.Results.sim_ss_d_type;  % D矩阵是否辨识
    cov_cross_type = parser.Results.cov_cross_type;  % 是否具有协方差
    sim_x_size_type = parser.Results.sim_x_size_type;  % 状态变量估计方式
    cov_order_type = parser.Results.cov_order_type;  % 协方差估计方式
    sim_x_size = parser.Results.sim_x_size;  % 预定义状态变量大小
    cov_order = parser.Results.cov_order;  % 预定义自协方差估计阶数
    cov_est_type = parser.Results.cov_est_type;  % 定义方差估计方法

    % 运行时准备
    % 清除临时存储
    clear idenISIM idenDCISSIMRunner
    % 准备辨识频率
    frequencies_bound = persistentExcitationCondition(x_size_upbound, u_size);
    switch isim_excitation_type
        case 'full'  % 全激励频率
            frequencies = fullExcitationFrequencies(period_samples);
        case 'reduced'  % 部分激励频率
            un_test = un_test(:, cutted_periods*period_samples+1:(cutted_periods+1)*period_samples);
            frequencies = reducedExcitationFrequencies(un_test, period_samples, frequencies_bound);
        otherwise, frequencies = 0;
    end
    % 初始化回归元
    mat_s = matSIniter(period_samples, frequencies);
    regressor = idenRegressor(period_samples, frequencies, 1, 'recursive');
    % 对于在线辨识 - 初始化参数
    if strcmp(dcissim_type, 'online') || strcmp(dcissim_type, 'online-test')
        v_size = size(mat_s, 1);
        idenISIM(zeros(y_size, 1), zeros(u_size, 1), zeros(v_size, 1), 'recursive');
    end
    % 返回值
    ret_struct = struct('y_size', y_size, 'u_size', u_size, 'period_samples', period_samples, 'cutted_periods', cutted_periods, ...
        'dcissim_type', dcissim_type, 'isim_excitation_type', isim_excitation_type, 'x_size_upbound', x_size_upbound, 'sim_ss_bdx_type', sim_ss_bdx_type, 'sim_ss_d_type', sim_ss_d_type, 'cov_cross_type', cov_cross_type, ...
        'frequencies', frequencies, 'mat_s', mat_s, ...
        'regressor', regressor, 'sim_x_size_type', sim_x_size_type, 'cov_order_type', cov_order_type, 'sim_x_size', sim_x_size, 'cov_order', cov_order, ...
        'cov_est_type', cov_est_type);

end

function frequencies_bound = persistentExcitationCondition(x_size_upbound, u_size)
% 持续激励条件(未知信号阶数, 用上界近似)
    frequencies_bound = u_size*2*x_size_upbound;
end

function frequencies = fullExcitationFrequencies(period_samples)
% 全频率
    frequencies = 0:fix(period_samples/2);
end

function frequencies = reducedExcitationFrequencies(un_test, period_samples, frequencies_bound)
% 选择频率
% 方法: 满足持续激励条件(再多几个频率点), 按照输入信号对应fft的幅值从高到低依次选择
% 最后总是加上零频率

    % 参数计算
    harmonic_size = fix(period_samples/2)+1;

    % 升维度
    frequencies_bound = min(round(frequencies_bound*3), harmonic_size);  % 经验参数
    % FFT
    ufreq = fft(un_test, period_samples, 2);
    ufreq_abs = abs(ufreq(:, 1:harmonic_size));

    % 局部极值样条包络交点法提取频率点
    frequencies_selection = peakFinder(ufreq_abs, frequencies_bound);
    
    % plot
    % figure; plot(ufreq_abs); hold on; stem(max(ufreq_abs, [], 'all')*frequencies_selection);
    % 返回值
    frequencies = 0:harmonic_size-1;
    frequencies = frequencies(logical(frequencies_selection));
    % 加上零频率
    if frequencies(1) ~= 0, frequencies = [0 frequencies]; end

end

function selector = peakFinder(mat_curve, peak_number)
% 直接从离群值点中取(保证数目足够, 从高到低依次取)
    
    % 参数计算 & 准备
    % channel_size = size(mat_curve, 1);
    curve_size = size(mat_curve, 2);
    selector = zeros(1, curve_size);
    thresold = 100 * (1 - (peak_number/curve_size)) - 3;

    % 离群值
    tf_ori = isoutlier(mat_curve, 'percentiles', [0 thresold], 2);
    % 取所有通道的离群值所在位置的最大值点
    curve_max = max(mat_curve, [], 1);
    curve_max(tf_ori == 0) = 0;
    % 将联合所有通道的点排序
    [~, sort_idx] = sort(curve_max, 'descend');
    selector(sort_idx(1:peak_number)) = 1;

    % % 样条插值包络 + 峰值提取方法寻找采样频率点, 默认曲线非负
    % % 准备数组
    % select_number = 0;
    % pre_selector = zeros(1, curve_size);  % 临时存储点
    % selector = zeros(1, curve_size);  % 最终选择点
    % 
    % % 按照预定义的尺寸提取包络峰值, 若最后点数不够再缩小查找尺度
    % envelope_np = ceil(curve_size/(3*peak_number));  % 经验参数
    % while select_number < peak_number
    %     % 按照当前查找尺度查找并记录peak点
    %     pre_selector_value_mat = zeros(size(mat_curve));  % 预记录数组
    %     for iter_channel = 1:channel_size
    %         % 局部极值点的样条插值包络
    %         [iter_curve_envelope, ~] = envelope(mat_curve(iter_channel, :), envelope_np, 'peak');
    %         % 查找交点
    %         iter_curve_envelope_peak_location = iter_curve_envelope == mat_curve(iter_channel, :);
    %         iter_curve_envelope_peak_value = iter_curve_envelope(iter_curve_envelope_peak_location);
    %         pre_selector_value_mat(iter_channel, iter_curve_envelope_peak_location) = iter_curve_envelope_peak_value;
    %         % [iter_curve_envelope_peak_value, iter_curve_envelope_peak_location] = findpeaks(iter_curve_envelope);  % 峰值在包络线的np过大时存在问题, 考虑np是由大到小变化的, 因此使用交点
    %     end
    %     % 从所有通道中选择最大值点
    %     pre_selector_value_vec = max(pre_selector_value_mat, [], 1);
    %     [pre_selector_value_vec_sort, pre_selector_value_vec_sort_location] = sort(pre_selector_value_vec, 'descend');
    %     pre_selector(pre_selector_value_vec_sort_location(pre_selector_value_vec_sort > 0)) = pre_selector_value_vec_sort(pre_selector_value_vec_sort > 0);
    %     % 更新
    %     envelope_np = ceil(envelope_np/5);  % 5为经验参数
    %     select_number = sum(pre_selector > 0);
    % end
    % 
    % % 若点过多, 则根据记录的幅值重选择
    % [~, select_sort_location] = sort(pre_selector, 'descend');
    % selector(select_sort_location(1:peak_number)) = 1;
    
end

function mat_s = matSIniter(period_samples, frequencies)
% S矩阵计算

    % 参数准备
    frequencies = unique(frequencies, 'sorted');  % 默认升序排序

    % 参数计算
    omega = (2*pi)/period_samples;
    frequencies_size = length(frequencies);
    if frequencies(1) == 0
        v_size = 2*frequencies_size-1;
        harmonic_size = frequencies_size-1;
        harmonic_frequencies = frequencies(2:end);
    else
        v_size = 2*frequencies_size;
        harmonic_size = frequencies_size;
        harmonic_frequencies = frequencies;
    end
    % 初始化返回值
    mat_s = zeros(v_size, v_size);

    % 矩阵计算
    % 直流部分(如有)
    if frequencies(1) == 0
        mat_s(1, 1) = 1;
        location_base = 1;
    else
        location_base = 0;
    end
    % 谐波部分
    for iter = 1:harmonic_size
        omega_iter = omega*harmonic_frequencies(iter);
        mat_s(location_base+1:location_base+2, location_base+1:location_base+2) = ...
            [cos(omega_iter) sin(omega_iter); -sin(omega_iter) cos(omega_iter)];
        location_base = location_base + 2;
    end
end
