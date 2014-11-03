function yhat = BuxtonAnalytical(t,params,XASL)
% Buxton model (using Hrabe boxcar formululation)
% XASL 1 = PASL
% XASL 2 = CASL

deltaTa = params(1);
taua = params(2);
T_1t = params(3);
T_1b = params(4);
fcalib = params(5);
alpha = params(6);
lambda = params(7);

tau1 = deltaTa; tau2 = deltaTa + taua;

epoch1 = t < tau1;
epoch2 = (t >= tau1) & (t <= tau2);
epoch3 = t > tau2;

yhat = zeros(size(t));

T_1app = 1/((1/T_1t) + (fcalib/lambda));

if XASL == 2   %CASL
    %F = 2*alpha;
    F = 2*fcalib;
    yhat(epoch1) = 0;
    
    yhat(epoch2) = F*T_1app*exp(-tau1/T_1b) .* ( 1 - exp(-(t(epoch2)-tau1)/T_1app) );
    yhat(epoch3) = F*T_1app*exp(-tau1/T_1b) .* exp(-(t(epoch3)-tau2)/T_1app) .* ( 1 - exp(-taua/T_1app) );
    
elseif XASL == 1   %PASL
    
    F = 2*fcalib/lambda*exp(-t/T_1b);
    R = (1/T_1app)-(1/T_1b);
    
    yhat(epoch1) = 0;
    
    yhat(epoch2) = (F(epoch2)/R).*(exp(R*t(epoch2))).*(exp(-R*tau1) - exp(-R*t(epoch2)));
    yhat(epoch3) = (F(epoch3)/R).*(exp(R*t(epoch3))).*(exp(-R*tau1) - exp(-R*tau2));
end

