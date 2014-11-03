function y=Convolution(f1,f2,t)

for i = 2:numel(t)
    y(i) = trapz(t(1:i), f1(1:i).*f2(i-(1:i)+1));
end

y = y + f1(1)*f2(1);

