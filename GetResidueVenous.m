function rm = GetResidueVenous(t,params,extra_params,value)
    switch value
        case 1   % SPA Model
            rm = GetrmSPAvenous(t,params,extra_params);
        case 2
            rm = GetrmBuxtonVen(t,params);
    end
end

