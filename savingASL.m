function savingASL(hObject,eventdata,handles)

addpath('../fsl_matlab')
%% Obtaining the data

%% Modality
 XASL = get(handles.popupmenu1,'Value');

%% Parameters
 At = str2double(get(handles.edit1,'String'));
 taua = str2double(get(handles.edit7,'String'));
 T1_t = str2double(get(handles.edit2,'String'));
 T1_b = str2double(get(handles.edit3,'String'));
 fcalib = str2double(get(handles.edit4,'String'))/6000;
 alfa = str2double(get(handles.edit6,'String'));
 lambda = str2double(get(handles.edit5,'String'));

%% Buxton analytical solution
 Buxton = get(handles.checkbox7,'Value');

%% Dispersion
 Disp = get(handles.popupmenu2,'Value');
 extra_paramsD(1) = str2double(get(handles.edit8,'String'));
 extra_paramsD(2) = eval(get(handles.edit9,'String'));
 extra_paramsD(3) = eval(get(handles.edit10,'String'));
 extra_paramsD(4) = str2double(get(handles.edit11,'String'));
 extra_paramsD(5) = str2double(get(handles.edit12,'String'));
 extra_paramsD(6) = str2double(get(handles.edit31,'String'));
 extra_paramsD(7) = str2double(get(handles.edit32,'String'));

%% Residue
 Art_comp = get(handles.checkbox6,'Value');
 if Art_comp
     ArtBlVol = str2double(get(handles.edit27,'String'));
 else
     ArtBlVol = 0;
 end

 Imp_comp1 = get(handles.checkbox2,'Value');
 if Imp_comp1
     tau_pc = str2double(get(handles.edit14,'String'));
 else
     tau_pc = 0;
 end

 Tis_comp = get(handles.checkbox3,'Value');
 if Tis_comp
     tau_c = str2double(get(handles.edit15,'String'));
     Res_tissue = get(handles.popupmenu3,'Value');
     extra_paramsR(1) = str2double(get(handles.edit17,'String'));
     extra_paramsR(2) = str2double(get(handles.edit18,'String'));
     extra_paramsR(3) = str2double(get(handles.edit19,'String'));
 else
     tau_c = 0;
 end

 Imp_comp2 = get(handles.checkbox4,'Value');
 if Imp_comp2
     tau_v = str2double(get(handles.edit16,'String'));
     Res_venous = get(handles.popupmenu4,'Value');
     extra_paramsV(1) = str2double(get(handles.edit22,'String'));
 else
     tau_v = 0;
 end

 Vein_comp = get(handles.checkbox5,'Value');
 if Vein_comp
     VeinBlVol = str2double(get(handles.edit28,'String'));
 else
     VeinBlVol = 0;
 end

 params=[At taua T1_t T1_b fcalib alfa lambda tau_pc tau_c tau_v];

%% t vector
set(handles.edit21,'String','');
t1 = eval(get(handles.edit26,'String'));
t = interp(t1,1000);

    %% Vectors
deltat = eval(get(handles.edit24,'String'));
f = eval(get(handles.edit25,'String'))/6000;
    

Signal = zeros(length(deltat),length(f),length(t)); 

for m=1:length(deltat)   
    for j=1:length(f)
        params(1) = deltat(m);
        %params(5) = f(j);
        deltat_matrix(m,j) = deltat(m);
        f_matrix(m,j) = f(j);
        
        %% Buxton analytical solution
        if Buxton
            Signal(m,j,:) = BuxtonAnalytical(t,params,XASL);
        else
        %% Compartment convolutions

            %% Arterial compartment
            if Art_comp
                SignalArt = 2*ArtBlVol*GetAIF(t,params,extra_paramsD,XASL,Disp);
            else
                SignalArt = zeros(1,length(t));
            end

            %% Impermeable arteriole compartment
            if Imp_comp1
                AIF_impcomp1 = GetAIF(t,params,extra_paramsD,XASL,Disp); 
                Imprm = GetrmImpPreCap(t,params);
                SignalPreCap = 2*f(j)*Convolution(AIF_impcomp1,Imprm,t);
                params(1) = params(1)+tau_pc;
            else 
                SignalPreCap = zeros(1,length(t));
            end

            %% Tissue compartment
            if Tis_comp
                AIF_tiscomp = GetAIF(t,params,extra_paramsD,XASL,Disp);
                Tissuerm = GetResidueTissue(t,params,extra_paramsR,Res_tissue);
                SignalTissue = 2*f(j)*Convolution(AIF_tiscomp,Tissuerm,t);
                params(1) = params(1)+tau_c;
            else
                SignalTissue = zeros(1,length(t));
            end

            %% Impermeable venous compartment
            if Imp_comp2
                AIF_impcomp2 =  GetAIF(t,params,extra_paramsD,XASL,Disp); 
                Venousrm = GetResidueVenous(t,params,extra_paramsV,Res_venous);
                SignalVenous = 2*f(j)*Convolution(AIF_impcomp2,Venousrm,t);
                params(1) = params(1)+tau_v;
            else
                SignalVenous = zeros(1,length(t));
            end

            %% Vein compartment
            if Vein_comp
                SignalVein = 2*VeinBlVol*GetAIF(t,params,extra_paramsD,XASL,Disp);
            else
                SignalVein = zeros(1,length(t));
            end

            %% Noise
            noisy = get(handles.checkbox1,'Value');
            if noisy
                stand = eval(get(handles.edit13,'String'));
                n = stand*randn(m,j,length(t));
            else
                n = zeros(1,length(t));
            end
            SignalArt = shiftdim(SignalArt,-1);
            SignalPreCap = shiftdim(SignalPreCap,-1);
            SignalTissue = shiftdim(SignalTissue,-1);
            SignalVenous = shiftdim(SignalVenous,-1);
            SignalVein = shiftdim(SignalVein,-1);
            n = shiftdim(n,-1);

            %% Signal
            Signal(m,j,:) = SignalArt+SignalPreCap+SignalTissue+SignalVenous+SignalVein+n;
        end
    end
end

Signal = Signal(:,:,1:1000:end);
Signal = Signal(:,:,2:end);
t = t(1:1000:end);
t = t(2:end)

Signal = shiftdim(Signal,-1);
deltat_M = shiftdim(deltat_matrix,-1);
f_M = shiftdim(f_matrix,-1);

delete('../fsl_basil/Signal.nii.gz')
save_avw(Signal,'../fsl_basil/Signal','f',[1 1 1 3]);
mask = ones(size(Signal(:,:,:,1)));
 
save_avw(mask,'../fsl_basil/mask','f',[1 1 1 3]);

save_avw(deltat_M,'../fsl_basil/deltatM','f',[1 1 1 3]);
save_avw(f_M,'../fsl_basil/fM','f',[1 1 1 3]);

'Signal was saved'




