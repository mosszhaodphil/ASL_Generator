function AIF = GetaBuxton(t,params,XASL)
%% 

%XASL = 1, PASL
%XASL = 2, CASL

%Define global parameters
deltaTa = params(1);
taua = params(2);
T_1t = params(3);
T_1b = params(4);
fcalib = params(5);
alfa = params(6);

%PASL
if XASL==1
    for j=1:length(t)
    if (deltaTa)>t(j)
        AIF(j) = 0;
    elseif (deltaTa)<=t(j) && t(j)<(deltaTa+taua)
        AIF(j) = alfa*(exp(-t(j)/T_1b));
    elseif t(j)>=(deltaTa+taua);
        AIF(j) = 0;
    end
    end
    
%CASL    
elseif XASL==2
    for j=1:length(t)
    if (deltaTa)>t(j)
        AIF(j) = 0;
    elseif (deltaTa)<=t(j) && t(j)<(deltaTa+taua)
        AIF(j) = alfa*exp(-deltaTa/T_1b);
    elseif t(j)>=(deltaTa+taua);
        AIF(j) = 0;
    end
    end
end
end

        



        
    


