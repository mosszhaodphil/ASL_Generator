function KSLresidue = GetrmKSL(t,params,extra_params)
%% St. Lawrence, K.S., Lee, T.-Y., 1998. An Adiabatic Approximation to the
%% Tissue Homogeneity Model for Water Exchange in the Brain: I. Theoretical Derivation: Journal of Cerebral Blood Flow & Metabolism 1365–1377.

deltaTa = params(1);
taua = params(2);
T_1 = params(3);
T_1b = params(4);
fcalib = params(5)*6000;
alfa = params(6);
lambda = params(7);
taupc = params(8);
tauc   = params(9);
tauv = params(10);

%E = 0.5; %fraction of tracer extracted into EVS during single capillary transit
PS = extra_params(1); 
V_i = extra_params(2); %mL/100g distribution volume of tracer in IVS
V_e = extra_params(3); %mL/100g distribution volume of tracer in EVS
F = fcalib;

E = 1-(exp(-PS/F));
k_adb = (E*fcalib)/V_e;
tc=V_i/fcalib;
for j=1:length(t)
    if t<tc
       KSLresidue = 1;
    else
       KSLresidue = E*exp(-k_adb*(t-tc));
    end
end


end

    