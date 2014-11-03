function AIF = GetaGaussian(t,params,sigma1,XASL)
%% Hrabe, J., Lewis, D.P., 2004. Two analytical solutions for a model of
%% pulsed arterial spin labeling with randomized blood arrival times. J. Magn. Reson. 167, 49–55.

%XASL = 1, PASL
%XASL = 2, CASL

%Define global parameters
deltaTa = params(1);
taua = params(2);
T_1t = params(3);
T_1b = params(4);
fcalib = params(5);
alfa = params(6);

%Define local parameters

% Low 
%sigma1 = 0.02;
% Medium 
%sigma1 = 0.05;
% High 
%sigma1 = 0.1;
sigma2 = (sqrt(((deltaTa)+taua)/(deltaTa)))*sigma1;

if XASL==1       %PASL
    
    AIF = 0.5*exp(-t/T_1b).*( erf((t-deltaTa)/(sqrt(2)*sigma1)) - erf( (t-deltaTa-taua)/(sqrt(2)*sigma2) ) );
    
elseif XASL==2    %CASL
    
    a = 1/sigma1^2;
    b = 1/T_1b;

    Q = (2*a.*(deltaTa-t)-b)./(2*a);
    S = (a.*(deltaTa-t).^2+b*t)./a;

    pt = t;
    pt(t>taua) = taua;

    AIF = 0.5*exp(-deltaTa/T_1b)*exp(Q.^2-S).*( erf( sqrt(a).*(pt+Q) ) - erf( sqrt(a).*Q ) );
end

