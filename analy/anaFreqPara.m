function [f_amp, f_phase, f_freq] = anaFreqPara(t_data, t_srate, f_freq)
%ANAFREQPARA 通过FFT获得谐波系数
% t_data的形状为(t_length, sizes)
% f_freq为一向量, 长度为f_length, 有效值为[0, 1/2T)
% f_amp和f_phase的形状均为(vx_length, sizes)

    % 计算参数
    t_slength = size(t_data, 1);  % 采样点数
    f_arate = 1/t_srate;  % (2倍)频带上限

    % 验证提取频率参数
    f_freq = reshape(f_freq, [], 1);  % 重构为列向量
    f_freq = sort(f_freq, 'ascend');  % 按升序排序
    % 取有效区间且符合采样率的值
    f_freqloc = freqApproximator(f_arate, t_slength, f_freq);
    f_freq = freqRestore(f_arate, t_slength, f_freqloc);
    
    % FFT变换
    f_data = fft(t_data);
    f_data = f_data(f_freqloc, :);  % 取有效值

    % 变换参数
    f_an =  2*real(f_data);
    f_bn = -2*imag(f_data);
    
    % 变换为参数形式
    f_amp = sqrt(f_an.^2 + f_bn.^2);
    f_phase = atand(f_bn ./ f_an);
    % 处理零频率
    if f_freqloc(1) == 1, f_amp(1) = f_amp(1) / 2; end

end

%   找到vx_freq在频率轴上最接近的近似点
function f_freqloc = freqApproximator(f_arate, t_slength, f_freq)
    % 变换频率点为[0, t_length)范围内的整数
    f_freq = f_freq./(f_arate/t_slength);

    % 四舍五入到最接近的近似点
    f_freqapprox = round(f_freq);
    % 取有效区间内的值
    f_freqapprox = f_freqapprox((f_freqapprox >= 0) & (f_freqapprox < fix(0.5*t_slength)));
    % 求每个频率上的唯一值
    f_freqapprox = unique(f_freqapprox, 'stable', 'rows');

    % 变换为坐标并提取对应的近似频率
    f_freqloc = f_freqapprox + 1;  % 从0开始的"整数化频率"变为从1开始的坐标
end

function f_freq = freqRestore(f_arate, t_slength, f_freqloc)
    % 构造频率轴
    f_axis = (0:t_slength-1).*(f_arate/t_slength);
    % 注: 时间轴为 t_axis = (0:t_slength-1).*t_srate;  % 时间轴
    % 返回值
    f_freq = f_axis(f_freqloc);
end

