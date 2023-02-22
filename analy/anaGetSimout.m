function retmat = anaGetSimout(simout, xselect)
%ANAGETSIMOUT 获取simulink输出数据
    
    retmat = simout.Data;
    if ~isempty(xselect), retmat = retmat(:, logical(xselect)); end

end

