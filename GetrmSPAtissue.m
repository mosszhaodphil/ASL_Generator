function rm = GetrmSPAtissue(t,params,extra_params)
%%St Lawrence, K.S., Frank, J.A., McLaughlin, A.C., 2000. Effect of
%%restricted water exchange on cerebral blood flow values calculated with arterial spin tagging: a theoretical investigation. Magn Reson Med 44, 440–449.

deltaTa = params(1);
taua = params(2);
T_1 = params(3);
T_1b = params(4);
fcalib = params(5);
taupc = params(6);
tauc   = params(7);
tauv = params(8);


%Grey matter assumptions from KSL paper
%PS = 200; %cc/100g/min
%V_c =2; %cc/100g
%V_b = 90; %cc/100g
%F = 75; %cc/100g.min
PS = extra_params(1);
V_c = extra_params(2);
V_b = extra_params(3);
F = fcalib;

R1b = T_1b; %/s (1)
R1c = T_1; %/s (0.83)

dR = R1c-R1b;

k_c = PS/V_c;
a_c = k_c + R1c;
beta = 1/(1+(dR/k_c));
E_r = 1-(exp((-PS/F)-(dR*tauc)));
E = 1-(exp(-PS/F));
    for j=1:length(t)
        if t(j)<=tauc
            rm(j) = beta*(exp(-R1b*t(j)))+((1-beta)*(exp(-a_c*t(j))));
        elseif t(j)>tauc 
            rm(j) = E_r*beta*(exp(-R1b*t(j)));
        end
    end
end
