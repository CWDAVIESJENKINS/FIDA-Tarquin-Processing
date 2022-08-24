function[OUT , List]=DWMRS_ReadDat_GFlip(file_Met , OutDir)


%% Initialise some parameters and load data

iterin=10; tmaxin=0.25; BadAv_nsd=2;

% Read in both datasets:
RawData=io_loadspec_twix(file_Met);

%% Steal some fields from twix data.
%Grab number of directions and diffusion gradients from metabolite acquisition
Ndir = RawData.p.alFree{4};
Ng = RawData.p.alFree{5};
RepSize = Ndir*Ng+1;
Nav = RawData.sz(3)./RepSize;

G = cell2mat(RawData.p.alFree(22:21+Ng));

RawData.p.RawDataLocation = file_Met; % Add data locations to structure


%% Get relative coil phases & amplitudes using B0 unsuppressed data
% Coils phased using 1st point in time domain, and weighting based on S/N^2.
% Hall et al. Neuroimage 2014
coilcombos = op_getcoilcombos(RawData,1,'h');
CC = op_addrcvrs(RawData,1,'h',coilcombos);

%CC.fids(:,497:end) = [];CC.specs(:,497:end) = [];CC.sz(2)=496;CC.averages = 496;% Cock up

%% Sepperate out by diffusion condition
Counter = 1;
for J=RepSize:-1:1 % Loop over and seperate each diffusion codition
    Ind_Met = (RepSize:RepSize:CC.sz(2))-J+1;
    Spec_cc{Counter} = Subset(CC,Ind_Met);
    Counter = Counter+1;
end
Spec2 = reshape(Spec_cc(2:end),2,Ng/2*Ndir)';
Spec_cc(2:end) = [];

for J=1:length(Spec2)
    Spec_cc{J+1} = Spec2{J,1};
    for K=1:Nav
        Mag = sqrt(abs(Spec2{J,1}.fids(:,K)) .* abs(Spec2{J,2}.fids(:,K)));
        Phase = mean([unwrap(angle(Spec2{J,1}.fids(:,K))),unwrap(angle(Spec2{J,2}.fids(:,K)))],2);
        Spec_cc{J+1}.fids(:,K) = Mag.*exp(Phase*1i);
        
    end
    Spec_cc{J+1}.specs = fftshift(ifft(Spec_cc{J+1}.fids,[],Spec_cc{J+1}.dims.t),Spec_cc{J+1}.dims.t);
end
RepSize = Ng/2*Ndir + 1;



%% Averaging, ECC, and freq/phase correction
Spec_Av=cell(1,RepSize);
OutMat=zeros(RepSize,CC.sz(1));

for J=1:RepSize % Again loop over diffusion conditions
    
    Rep = Spec_cc{J};
    %[Rep,Metric{J},badAverages{J}]=op_rmbadaverages(Spec_cc{J},3,'t');
    
    fsPoly=100; phsPoly=1000; iter=1;
    fsCum=zeros(Rep.sz(Rep.dims.averages),1);
    phsCum=zeros(Rep.sz(Rep.dims.averages),1);
    
    while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<=iterin
        close all
        tmax=tmaxin+0.03*randn(1);
        % FRange for Spectral registration
        fmin=1.7+0.1*randn(1); fmaxarray=[3.5+0.1*randn(1,3),4+0.1*randn(1,3),5.5+0.1*randn(1,1)]; fmax=fmaxarray(randi(6,1));
        
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
    
    OUT = Spec_Av;
    OutMat(J,:)=real(Spec_Av{J}.specs);
end



%% If Out directory specified, then save combined spectra in repsective directories & sumary plots
if exist('OutDir','var')
    Home=pwd;cd(OutDir)

    mkdir('Met');
    for J=1:RepSize
        cd([OutDir,'/Met'])
        io_writelcm(Spec_Av{J} , sprintf('Spec_%03i.txt',J) , Spec_Av{J}.te);
    end
    cd(Home)

    List.Met = FindFiles([OutDir,'/Met/*.txt']);List.Met=List.Met(:,1);
end


end %function

%% Get subset of spectra specified by Ind, and create a new structure.
function[Struct_Out] = Subset(Struct_In , Ind)
    Struct_Out = Struct_In;
    Struct_Out.specs = Struct_Out.specs(:,Ind);
    Struct_Out.fids = Struct_Out.fids(:,Ind);
    Struct_Out.sz = size(Struct_Out.specs);
    Struct_Out.averages = length(Ind);
    Struct_Out.p.Ind = Ind; % Add this for now to debug Diffusion conditions
end

%% Perform automatic phi0 and frequency shift for water unsuppressed data
function[out] = FnPh_w(Spec) 
    
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

%% Automatic phi0 & fshift correction for metabolite spectra using Creatine Peak.
%  (May not be optimal in all cases)
function[out] = FnPh_m(Spec)
    Temp_zp=op_zeropad(Spec,16);
    [Temp_ph,ph0]=op_autophase(Spec,2.9,3.1);
    Temp_zp=op_addphase(Temp_zp,ph0);
    
    % Frequency shift all spectra so that Creatine appears at 3.027 ppm:
    [~,frqShift]=op_ppmref(Temp_zp,2.8,3.1,3.027);
    out=op_freqshift(Temp_ph,frqShift);% Final metabolite spectra
end
