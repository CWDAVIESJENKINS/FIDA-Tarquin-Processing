function[Spec_Av]=DWMRS_ReadDcm(IN , OutDir)

iterin=5; tmaxin=0.2; aaDomain='f';


List  = dir(sprintf('%s/*.dcm',IN));

Spec_cc=cell(length(List));
for J=1:length(List)
    Spec_cc{J}=io_loadspec_IMA(List(J).name , 2.8936 , 4000);
end


Ndir = 3; Ng = 6; RepSize = Ndir*Ng+1;Nav = length(List)/RepSize;


Counter = 1;
for J=RepSize:-1:1 % Loop over and seperate each diffusion codition
    Ind_Met = (RepSize:RepSize:CC.sz(2))-J+1;
    Spec_cc{Counter} = Subset(Spec_CC,Ind_Met);
    Counter = Counter+1;
end

% Would normally remove motion corrupted averages before summing, but not sure it
% will work for high b-values. However, I could try:
% 1) Model Frequency & phase drift using only the B0 acquisitions
% 2) Apply this to all intermediate conditions

for J=1:RepSize % Again loop over diffusion
    Rep = Spec_cc{J};
    
    fsPoly=100; phsPoly=1000; iter=1;
    fsCum=zeros(Rep.sz(Rep.dims.averages),1);
    phsCum=zeros(Rep.sz(Rep.dims.averages),1);
    
    while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
        close all
        tmax=tmaxin+0.04*randn(1);
        % FRange for SpecReg. Range between ~1.8-3.5
        fmin=1.8+0.04*randn(1); fmaxarray=[3.5,3.5,4.2,4.2,7.5,4.2]; fmax=fmaxarray(randi(6,1));
        
        % Perform spectral registration in the frequency domain to correct frequency and phase drifts.
        % Could use 'a' as final input, and supply the highest snr average
        [Rep,fs,phs]=op_alignAverages_fd(Rep,fmin,fmax,tmax,'n');
        
        % Freq & Phase shifts. Could run this for only B0 as a comparison.
        fsCum=fsCum+fs(:,1);
        phsCum=phsCum+phs(:,1);
        
        fsPoly=polyfit([1:Rep.sz(Rep.dims.averages)]',fs(:,1),1);
        phsPoly=polyfit([1:Rep.sz(Rep.dims.averages)]',phs(:,1),1);

        iter=iter+1;
    end

    % Average spectra, and leftshift if needed.
    Spec_Av{J} = op_leftshift(op_averaging(Rep),Rep.pointsToLeftshift);
    
 
    % Perform frequency and phase alignement using the water peak for the
    % water unsuppressed acquisition.
    Spec_Av{J} = FnPh_m(Spec_Av{J});
    
          
end

if exist('OutDir','var')
    cd(OutDir)
    mkdir('Met');mkdir('Ref');mkdir('Met_ECC')
    for J=1:RepSize
        cd([OutDir,'/Met'])
        io_writejmrui(Spec_Av{J} , sprintf('Spec_%i.txt',J));
        cd([OutDir,'/Ref'])
        io_writejmrui(Specw_Av{J} , sprintf('Spec_%i.txt',J));
        cd([OutDir,'/Met_ECC'])
        io_writejmrui(Spec_Av_ECC{J} , sprintf('Spec_%i.txt',J));
    end
end

end

function[Struct_Out] = Subset(Struct_In , Ind)
% Grab a subset of spectra specified by vector, Ind, and create a new
% structure.
    Struct_Out = Struct_In;
    Struct_Out.specs = Struct_Out.specs(:,Ind);
    Struct_Out.fids = Struct_Out.fids(:,Ind);
    Struct_Out.sz = size(Struct_Out.specs);
    Struct_Out.averages = length(Ind);
    Struct_Out.p.Ind = Ind; % Add this for now to debug Diffusion conditions
end

function[Out] = ECC_fids(data, water_data)
% Perform eddy current corection using water unsuppressed fids.

Kphase=unwrap(angle(water_data.fids));
fids=zeros(size(data.fids));

K_data=abs(data.fids);
Kphase_data=unwrap(angle(data.fids));
Kphase_corr=(Kphase_data)-Kphase;
fids=K_data.* exp(1i*Kphase_corr);

Out = data;
Out.fids=double(fids);
Out.specs = fftshift(ifft(Out.fids,[],Out.dims.t),Out.dims.t);

end

function[out] = FnPh_w(Spec)
    %Automatic phi0 and frequency shift for water unsuppressed data
    
    % Zero order phase shift using water unsuppressed peak
    Tempw_zp=op_zeropad(Spec,16);
    indexw=find(abs(Tempw_zp.specs)==max(abs(Tempw_zp.specs(Tempw_zp.ppm>4 & Tempw_zp.ppm<5.4))));%FindMax
    ph0w=-phase(Tempw_zp.specs(indexw))*180/pi;
    outw_ph=op_addphase(Spec,ph0w);
    Tempw_zp=op_addphase(Tempw_zp,ph0w);
    
    % 
    [~,frqShiftw]=op_ppmref(Tempw_zp,4,5.5,4.65);
    out=op_freqshift(outw_ph,frqShiftw); 
end

function[out] = FnPh_m(Spec)
    % Automatic phi0 correction for metabolite spectra using Creatine Peak.
    % May not be optimal in all cases.
    Temp_zp=op_zeropad(Spec,16);
    [Temp_ph,ph0]=op_autophase(Spec,2.9,3.1);
    Temp_zp=op_addphase(Temp_zp,ph0);
    
    % Frequency shift all spectra so that Creatine appears at 3.027 ppm:
    [~,frqShift]=op_ppmref(Temp_zp,2.9,3.1,3.027);
    out=op_freqshift(Temp_ph,frqShift);% Final metabolite spectra
end