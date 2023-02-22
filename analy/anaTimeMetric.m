function [mcFF, mcVAF] = anaTimeMetric(ori_data, ind_data, datasize)
%ANATIMEMETRIC 分别给定原系统和辨识系统对同一输入的响应数据, 输出指标
% 指标:
% 1. Normalized fit Metric: F = 100 (1 - ||y-\hat{y}||_2 / ||y-\bar{y}||_2)
%    [L. Ljung, System Identification Theory for the User, 2nd ed.]
% 2. Variance Account for Metric: VAF = 100 (1 - var(y-\hat{y} / var(y)));
% (Akaike's Final Prediction Error?)

    % 计算参数
    sampleSize = size(ori_data, 1);
    
    % 指标
    % F
    mcFF_meanOri = mean(ori_data, 1);
    mcFF_Num = ori_data - ind_data;
    mcFF_Den = ori_data - repmat(mcFF_meanOri, [sampleSize, 1]);
    mcFF = vecnorm(mcFF_Num, 2, 1) ./ vecnorm(mcFF_Den, 2, 1);
    mcFF = 100 .* (1 - mcFF);
    % VAF (equivalent to r^2)
    mcVAF_Num = ori_data - ind_data;
    mcVAF = var(mcVAF_Num, 0, 1) ./ var(ori_data, 0, 1);
    mcVAF = 100 .* (1 - mcVAF);
    
    % 输出
    fprintf('Model Fit Metric:\n');
    specs = specBuilder("F: ", datasize);
    fprintf(specs, mcFF);
    specs = specBuilder("VAF: ", datasize);
    fprintf(specs, mcVAF);

end

function formatSpec = specBuilder(beginStr, vecSize)
    formatSpec = beginStr;
    for iter = 1:vecSize, formatSpec = formatSpec + "%.4f, "; end
    formatSpec = formatSpec + "\n";
end

