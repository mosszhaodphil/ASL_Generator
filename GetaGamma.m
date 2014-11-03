function AIF = GetaGamma(t,params,extra_params,XASL)
%% Chappell, M.A., Woolrich, M.W., Kazan, S., Jezzard, P., Payne, S.J.,
%% MacIntosh, B.J., 2013. Modeling dispersion in arterial spin labeling:
%% validation using dynamic angiographic measurements. Magn Reson Med 69, 563–570.

%XASL = 1, PASL
%XASL = 2, CASL

%Define global parameters
deltaTa = params(1);
taua = params(2);
T_1t = params(3);
T_1b = params(4);
fcalib = params(5);
alfa = params(6);

%Define local parameter

%sharpness/s^-1
%High 
%s = 10;
%Medium 
%s = exp(2);
s = extra_params(1);
%Low 
%s = 50;
%time to peak of curve/s
%sp = exp(-0.3);
sp = extra_params(2);
p = sp/s;
epoch1 = t<deltaTa;
epoch2 = (t>=deltaTa) & (t<=(deltaTa+taua));
t2 = t(epoch2);
epoch3 = t>deltaTa+taua;
t3 = t(epoch3);
    
%PASL
if XASL==1
    AIF(epoch1) = 0;
    AIF(epoch2) = exp(-t2/T_1b) .* (1 - gammainc(s*(t2-deltaTa),1+sp,'upper') ); %% 
    AIF(epoch3) = exp(-t3/T_1b) .* ( gammainc(s*(t3-deltaTa-taua),1+sp,'upper') - gammainc(s*(t3-deltaTa),1+sp,'upper') ); %  
%CASL
elseif XASL==2
    AIF(epoch1) = 0;
    AIF(epoch2) = exp(-t2/T_1b) .* (1 - gammainc((s+1/T1b)*(t2-delt),1+sp,'upper') ); %% 
    AIF(epoch3) = exp(-t3/T_1b) .* ( gammainc((s+1/T_1b)*(t3-deltaTa-taua),1+sp,'upper') - gammainc((s+1/T_1b)*(t3-deltaTa),1+sp,'upper') ); %       
end
end





        
    


