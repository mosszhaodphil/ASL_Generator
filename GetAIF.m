function AIF = GetAIF(t,params,extra_params,XASL,value)
%Select the AIF function

    switch value
    case 1   % None
        AIF = GetaBuxton(t,params,XASL);
    case 2   % Gaussian
        AIF = GetaGaussian(t,params,extra_params(1),XASL);
    case 3   % Gamma
        AIF = GetaGamma(t,params,extra_params(2:3),XASL);
    case 4   % GV
        AIF = GetaGV(t,params,extra_params(4:5),XASL);
    case 5   % Gallichan
        AIF = GetaGallichan(t,params,extra_params(6:7),XASL);
    end
end

