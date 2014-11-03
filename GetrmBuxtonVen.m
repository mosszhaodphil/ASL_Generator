%Residue and magnetisation decay Buxton

function rm = GetrmBuxtonVen(t,params)

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

for j=1:length(t)
    r(j) = 1-(exp((-(fcalib*t(j)))/lambda));
    m(j) = exp(-t(j)/T_1);
    rm(j) = r(j)*m(j);
end
