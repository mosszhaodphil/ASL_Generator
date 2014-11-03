function rm = GetrmSPA(t,params)
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
PS = 200; %cc/100g/min  (80 for white matter)
V_c =2; %cc/100g (1 for white matter)
V_b = 90; %cc/100g
F = 75; %cc/100g.min
R1b = 1; %/s  (1.7 for white matter) = T_1b
R1c = 0.83; %/s (0.83) = T_1

dT = T_1b-T_1;

k_c = PS/V_c;
a_c = k_c + T_1b;
beta = 1/(1+((dT*V_c)/PS));
E_r = 1-(exp((-PS/fcalib)-(dT*tauc)));
E = 1-(exp(-PS/fcalib));
rm=[];
    for j=1:length(t)
        if t(j)<=tauc
            rm(j) = beta*(exp(-T_1*t(j)))+((1-beta)*(exp(-a_c*t(j))));
        elseif t(j)>tauc && t(j)<=(tauc+tauv)
            rm(j) = E_r*beta*(exp(-T_1*t(j)))+((1-E)*(exp(-T_1b*t(j))));
        elseif t(j)>(tauc+tauv)
            rm(j) = E_r*beta*(exp(-T_1*t(j)));
        end
    end
end
