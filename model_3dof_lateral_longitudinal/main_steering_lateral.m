%% 准备仿真参数
% 系统参数
% 模型名称
modelname = 'model_3dof_steering_lateral_selfmade.slx';
% 仿真时间
para_dura = 3000;  % 实验时间
para_sim_step = 5e-3;  % 仿真最大时间步
para_resamples = 5;  % 数据采样倍数
% 信号说明
% y(t): "steering rate", "shaft rate", "steering angle", "shaft angle", "ydot", "yaw rate"
% u(t): "WhlAngF"
para_yselect = [0 0 1 1 1 1];
para_uselect = [1]; %#ok<NBRAK2> 
% 输入信号上下限
para_usat = [-5e2 5e2];
% 噪声参数 (nan表示不加噪声)
para_ysnr = 50*ones(sum(para_yselect), 1);
para_usnr = 50*para_uselect.';

% 辨识参数
% 激励周期
para_period = 95;  % 一周期真实秒数
para_sim_period_samples = para_period*fix(1/(para_sim_step*para_resamples));  % 一周期的采样数
% 激励频率
para_freqs = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
    1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 11];
% 算法参数
iden_cutstart = 1;  % 切除最开始的周期
iden_upxsize = 8;   % X的变量数上界
iden_xsize = 4;  % 固定情况下的X的变量数
iden_dtype = 'null';  % 是否假定D矩阵为0

% 时间计算
para_sim_samples = para_dura*fix(1/para_sim_step);  % 仿真采样数
para_sim_span = 0:para_sim_step:(para_sim_samples-1)*para_sim_step;  % 仿真时刻
% 构造参数结构体
para_ysize = sum(para_yselect);
para_usize = sum(para_uselect);
paraStruct = genSystemInfo( ...
    'Step', para_sim_step, 'Period', para_period, ...
    'YSize', para_ysize, 'USize', para_usize, 'Freqs', para_freqs, ...
    'USat', para_usat, 'USelect', para_uselect, ...
    'USnr', para_usnr, 'YSnr', para_ysnr);
% 随机初始参数
para_vinit = paraStruct.VInit;
% 仿真/辨识系统参数
para_vsize = paraStruct.VSize;
para_mats = paraStruct.Sv;
para_matu = paraStruct.Uv;
para_sim_usize = paraStruct.SimUSize;
para_sim_matu = paraStruct.SimUv;
para_sim_matk = paraStruct.SimKw;

%% 准备模型仿真
% 确定激励信号功率, 用于添加激励噪声
% 仿真激励系统
[~, ~, sim_ut] = genExcitation(paraStruct, para_sim_step, fix(para_sim_samples/2), para_vinit);
% 计算激励信号功率
ut_pwr = genPower(sim_ut, para_sim_usize);
% 调整运行参数
if sum(isnan(para_usnr)), sim_usnr = zeros(para_sim_usize, 1); 
else, sim_usnr = para_usnr; end
% 仿真控制信号
sim_iun = Simulink.Parameter(sum(isnan(para_usnr))); 

%% 运行仿真
sim_outputs = sim(modelname);

% 确定时刻和采样数
t = sim_outputs.vt.Time;
iden_samples = sim_outputs.vt.Length;
% 提取数据
unt = anaGetSimout(sim_outputs.ut, para_uselect);
yt = anaGetSimout(sim_outputs.yt, para_yselect);

% 添加输出噪声
ynt = genNoiseAdder(yt, para_ysnr, iden_samples, para_ysize);
% plot(t, ynt)

%% 辨识系统参数

% 辨识
est_Init = idenDCISSIMLaucher(unt.', 'y_size', size(ynt, 2), 'u_size', size(unt, 2), 'period_samples', para_sim_period_samples, 'cutted_periods', iden_cutstart, ...
'dcissim_type', 'offline', 'isim_excitation_type', 'reduced', 'x_size_upbound', iden_upxsize, 'sim_ss_bdx_type', 'analytical', 'sim_ss_d_type', iden_dtype, 'cov_cross_type', 'null', ...
'sim_x_size_type', 'fixed', 'sim_x_size', 5, 'cov_order_type', 'null');
est_Struct = idenDCISSIMRunner(est_Init, ynt.', unt.');


% 读取结果 (注意时域转换)
% 系统参数
est_SSModel = idss(est_Struct.A, est_Struct.B, est_Struct.C, est_Struct.D, ...
    zeros(size(est_Struct.A, 1), size(est_Struct.C, 1)), zeros(size(est_Struct.A, 1), 1), para_sim_step*para_resamples);
est_SSModel_conti = d2c(est_SSModel, 'zoh');
est_A = est_SSModel_conti.A; est_B = est_SSModel_conti.B; est_C = est_SSModel_conti.C; est_D = est_SSModel_conti.D;
% 系统初值
est_xinit = est_SSModel.x0;

% 输出A矩阵特征值
disp(['The eigenvalues of A are: ', mat2str(eig(est_A), 4)]);

%% 分析用无噪声输出
% 节约时间, 只分析几个周期
para_ana_dura = 3 * para_period;
para_ana_samples = para_ana_dura*fix(1/(para_sim_step*para_resamples));

% 原系统无扰动仿真
if ~sum(isnan(para_usnr))
    sim_usnr = zeros(para_sim_usize, 1); sim_iun = Simulink.Parameter(1);
    ana_ori_out = sim(modelname);
    ana_ori_yt = anaGetSimout(ana_ori_out.yt, para_yselect);
else
    ana_ori_yt = yt;
end
ana_ori_yt = ana_ori_yt(iden_cutstart*para_sim_period_samples+1:para_ana_samples, :);

%% 时域分析
% 辨识系统仿真 (可能差一位数据)
ana_t = 0:para_sim_step*para_resamples:(para_ana_samples-1)*para_sim_step*para_resamples;
[~, ~, ana_ut] = genExcitation(paraStruct, para_sim_step*para_resamples, para_ana_samples, para_vinit);
ana_ut = ana_ut(1:para_ana_samples, :);
ana_ind_yt = sim(est_SSModel, ana_ut, simOptions('InitialCondition', 'z'));  % 默认无噪声, 零初始值
ana_ind_yt = ana_ind_yt(iden_cutstart*para_sim_period_samples+1:para_ana_samples, :);

% 绘图
% anaPlot(ana_t, {ana_ori_yt}, {ana_ind_yt}, ...
%     'Time Response of Original and Indentified System (cISSIM)', 'plot');
% 计算拟合指标
[ana_time_FF, ana_time_VAF] = anaTimeMetric(ana_ori_yt, ana_ind_yt, para_ysize);

%% 输出信号频域分解对比
% 频率响应分析
% 原系统频率响应模型(MIMO对应每个输入的单独响应)
ana_freq_sampled = {0.1 11};  % 采样频率范围, 默认为30个频率点
ana_freq_orifrd = anaFreqFRD(modelname(1:end-4), para_sim_step, ana_freq_sampled, para_yselect, para_uselect, para_usat);
% 辨识系统状态空间模型
ana_freq_indss = ss(est_A, est_B, est_C, est_D);
% 绘图
figure;
bode(ana_freq_orifrd, 'b--', ana_freq_indss, 'r-');
title('Bode Plot of Original and Indentified System (cISSIM)');

%% N4SID辨识系统参数
% load('..\saved_results\steering_lateral_selfmade.mat')
sid_Data = iddata(ynt, unt, para_sim_step*para_resamples);
sid_Option = n4sidOptions('InitialState', 'zero', ...  % 初态为0
            'Focus', 'prediction');  % 最小化(一步)预测误差
sid_Model = n4sid(sid_Data, 'best', sid_Option);
% 输出A矩阵特征值
disp(['The eigenvalues of A in SID are: ', mat2str(eig(sid_Model.A), 4)]);

%% N4SID时域分析
% 无扰动输入输出信号(直接使用前面的信号)
ana_sid_yt = sim(sid_Model, ana_ut, simOptions('InitialCondition', 'z'));  % 默认无噪声, 零初始值
ana_sid_yt = ana_sid_yt(iden_cutstart*para_sim_period_samples+1:para_ana_samples, :);
% 绘图
% anaPlot(t, {ana_ori_yt}, {ana_sid_yt}, ...
%     'Time Response of Original and Indentified System (SID)', 'plot');
% 计算拟合指标
[ana_time_FF_sid, ana_time_VAF_sid] = anaTimeMetric(ana_ori_yt, ana_sid_yt, para_ysize);

%% N4SID频域分析
% 同样直接使用前面的数据
figure;
bode(ana_freq_orifrd, 'b--', sid_Model, 'r-');
title('Bode Plot of Original and Indentified System (SID)');

