function rm = GetrmSPAvenous(t,params,extra_params)
%%St Lawrence, K.S., Frank, J.A., McLaughlin, A.C., 2000. Effect of
%%restricted water exchange on cerebral blood flow values calculated with arterial spin tagging: a theoretical investigation. Magn Reson Med 44, 440–449.

deltaTa = params(1);
taua = params(2);
T_1 = params(3);
T_1b = params(4);
fcalib = params(5);
alfa = params(6);
lambda = params(7);
taupc = params(8);
tauc   = params(9);
tauv = params(10);


%Grey matter assumptions from KSL paper
% PS = 200; %cc/100g/min
% V_c =2; %cc/100g
% V_b = 90; %cc/100g
% F = 75; %cc/100g.min
PS = extra_params(1);
F = fcalib;

R1b = T_1b; %/s (1)
R1c = T_1; %/s (0.83)

for j=1:length(t)
    if t(j)<=tauc
        rm(j) = 0;
    elseif t(j)>tauc && t(j)<=(tauc+tauv)
        rm(j) = exp(-PS/F)*exp(-R1c*t(j));
    elseif t(j)>(tauc+tauv)
        rm(j) = 0;
    end
end

end
