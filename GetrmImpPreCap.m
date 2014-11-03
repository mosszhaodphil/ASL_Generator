%Residue function for impermeable component

function rmImp = GetrmImpPreCap(t,params)

deltaTa = params(1);
taua = params(2);
T_1 = params(3);
T_1b = params(4);
fcalib = params(5);
alpha = params(6);
lambda = params(7);
taupc = params(8);
tauc   = params(9);
tauv = params(10);

for j=1:length(t)
    if t(j)<taupc
        rmImp(j) = exp((-t(j))/T_1b);
    else
        rmImp(j) = 0;
    end
end
end
