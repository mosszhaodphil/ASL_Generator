function out = GetaGallichan(t,params,extra_params,XASL)
%% Gallichan, D., Jezzard, P., 2008. Modeling the effects of dispersion and
%% pulsatility of blood flow in pulsed arterial spin labeling. Magn Reson Med 60, 53�63.

%XASL = 1, PASL
%XASL = 2, CASL

deltaTa = params(1);
taua = params(2);
T_1t = params(3);
T_1b = params(4);
fcalib = params(5);
alfa = params(6);
xVm = extra_params(1);
gam = extra_params(2);

% x=0.04; %m (distance along pipe in longitudinal direction)
% V_m = 0.4; %m/s (max velocity through pipe)
% X = 0.1; %m (width along the tube)

if XASL == 1   %% PASL
    
    epoch1 = t<deltaTa;
    epoch2 = (t>=deltaTa) & (t<=(deltaTa+taua));
    t2 = t(epoch2);
    epoch3 = t>deltaTa+taua;
    t3 = t(epoch3);

    out(epoch1) = 0;
    out(epoch2) = 2 .* exp(-t2/T_1t) .* ( 1 - xVm./(t2+xVm-deltaTa) ); %% 
    out(epoch3) = 2 .* exp(-t3/T_1t) .* ( taua./(t3+xVm-deltaTa) ); %  


else  %% CASL
    
   out = zeros(1,length(t));

    dt = 0.001;
    for i=1:length(t)
        if t(i)>deltaTa % this implements the 'arrival timeshift'
            tshifted=t(i)-deltaTa;
            lambda = max(tshifted,gam):-dt:max(gam,max(0,tshifted-taua)); % note negative lambda doesn't make sense as this would be including label before it is created
            integrand = 1./lambda.^2.*exp(-lambda/T_1b);
            integrand(lambda==0)=0;
            integral = trapz(integrand)*dt;
            out(i) = gam*integral; % 
        end
    end

    % note the arrival time shift also needs a T1v decay with it
    out = out*exp(-deltaTa/T_1b); 
    
end

