function frd_est = anaFreqFRD(mdl_name, sigstep, sigfreqs, yselect, uselect, usaturation)
%ANAFREQFRD 求原系统(非线性)对不同频率信号的响应, 输出频率响应
% 参考:
% 频率响应概述: https://ww2.mathworks.cn/discovery/frequency-response.html
% frestimate: https://www.mathworks.com/help/releases/R2022a/slcontrol/ug/frestimate.html
    
    % 打开模型
    mdl_handle = load_system(mdl_name);
    
    % 参数计算
    % 输入
    sim_usize = length(uselect);
    frd_usize = sum(uselect);
    % 输出
    sim_ysize = length(yselect);
    frd_ysize = sum(yselect);
    % blk_name = [mdl_name '/' 'yt Bus'];
    % sim_yname = get_param(blk_name, 'OutputSignals');
    % sim_yname = split(sim_yname, ',');
    % sim_yname = extract(sim_yname, regexpPattern('(?=.)[(a-z|A-Z|0-9)]+$'));
    
    % 满偏
    usaturation = sort(usaturation, 2, 'ascend');

    % 定义输入输出
    ioLoc = 1;
    mdl_io(frd_usize+frd_ysize, 1) = linearize.IOPoint();
    % 输入
    blk_name = [mdl_name '/' 'ut Demux'];
    for iter_io = 1:sim_usize
        if uselect(iter_io)
            % mdl_io(iter_io) = linio(blockName, iter_io, 'input');  % 在当前multisine信号下的平衡点
            mdl_io(ioLoc) = linio(blk_name, iter_io, 'openinput');  % 在零输入信号下的平衡点
            ioLoc = ioLoc + 1;
        end
    end
    % 输出
    blk_name = [mdl_name '/' 'vec2elem'];
    for iter_io = 1:sim_ysize
        if yselect(iter_io)
            mdl_io(ioLoc) = linio(blk_name, iter_io, 'output');
            ioLoc = ioLoc + 1;
        end
    end

    % 寻找平衡点
    % 定义系统状态和输出
    mdl_spec = operspec(mdl_name);
    blk_name = [mdl_name '/' 'ut Demux'];
    for iter_io = 1:frd_usize
        mdl_spec = addoutputspec(mdl_spec, blk_name, mdl_io(iter_io).PortNumber);
        mdl_spec.Outputs(iter_io).Min = usaturation(iter_io, 1);
        mdl_spec.Outputs(iter_io).Max = usaturation(iter_io, 2);
    end
    % 计算平衡点
    opOpts = findopOptions('DisplayReport', 'off');
    mdl_op = findop(mdl_name, mdl_spec, opOpts);

    % 获取频率响应
    frd_sig = frest.createFixedTsSinestream(sigstep, sigfreqs);
    frd_est = frestimate(mdl_name, mdl_op, mdl_io, frd_sig);

    % 关闭模型
    close_system(mdl_handle);
    
end

