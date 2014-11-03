function rm = GetResidueTissue(t,params,extra_params,value)
%Select the AIF function
    switch value
        case 1
            rm = GetrmBuxton(t,params);
        case 2   % SPA Model
            rm = GetrmSPAtissue(t,params,extra_params);
        case 3
            rm = GetrmKSL(t,params,extra_params);
    end
end

