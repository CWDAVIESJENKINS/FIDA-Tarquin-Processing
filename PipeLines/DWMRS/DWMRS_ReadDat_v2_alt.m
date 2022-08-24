function[OUT , List]=DWMRS_ReadDat_v2_alt(IN , OutDir)
%function[Spec_Av, Specw_Av , Spec_Av_ECC]=DWMRS_ReadDat_v2(IN , OutDir)
% Read in raw dMRS data, perform coil coombination, seperate averages by
% diffussion condition, perform ECC, f & ph correction, and export spectra
% as LCModel .txt format for analysis.
% Input:    IN = cellarray containing {a,b}, where a is the location of the
%                metabolite acquisition, and b is the corresponding
%                referenceTest.sz = size(Test).specs
%           OutDir = Directory to save the averaged and correceted spectra
%                    for analysis.
% Output    OUT = structure containing processed metabolite data
%                 with/without ECC, and corresponding water references.

%% Initialise some parameters and load data

iterin=10; tmaxin=0.25;

file_Met = IN{1}; file_Ref= IN{2};

% Read in both datasets:
RawData=io_loadspec_twix(file_Met);
RawData_ref=io_loadspec_twix(file_Ref);

%Check that met and ref acquisitions use same sequence/special parameters
if ~(cell2mat(RawData.p.alFree) == cell2mat(RawData_ref.p.alFree))
    warning('Metabolite & reference scans may have different acquisition parameters!!!')
end

%% Steal some fields from twix data.
%Grab number of directions and diffusion gradients from metabolite acquisition
Ndir = RawData.p.alFree{4};
Ng = RawData.p.alFree{5};
BCond = Ndir*Ng+1;

RawData.p.RawDataLocation = file_Met; % Add data locations to structure
RawData_ref.p.RawDataLocation = file_Ref;

%% Get relative coil phases & amplitudes using B0 unsuppressed data
% Coils phased using 1st point in time domain, and weighting based on S/N^2.
% Hall et al. Neuroimage 2014
coilcombos = op_getcoilcombos(RawData_ref,1,'h');
CC = RemapSubSpec(op_addrcvrs(RawData,1,'h',coilcombos),BCond);
CCw = RemapSubSpec(op_addrcvrs(RawData_ref,1,'h',coilcombos),BCond);

Spec_av=CC;

fsPoly=100; phsPoly=1000; iter=1;
fsCum=zeros(Spec_av.sz(Spec_av.dims.averages),1);
phsCum=zeros(Spec_av.sz(Spec_av.dims.averages),1);

while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<=iterin
    close all
    tmax=tmaxin+0.03*randn(1);
    % FRange for Spectral registration
    fmin=1.7+0.1*randn(1); fmaxarray=[3.5+0.1*randn(1,3),4+0.1*randn(1,3),5.5+0.1*randn(1,1)]; fmax=fmaxarray(randi(6,1));

    % Perform spectral registration in the frequency domain to correct frequency and phase drifts.
    % Could use 'a' as final input, and supply the highest snr average
    [Spec_av,fs,phs]=op_alignAverages_fd(Spec_av,fmin,fmax,tmax,'n');

    % Freq & Phase shifts. Could run this for only B0 as a comparison.
    fsCum=fsCum+fs(:,1);
    phsCum=phsCum+phs(:,1);

    fsPoly=polyfit([1:Spec_av.sz(Spec_av.dims.averages)]',fs(:,1),1);
    phsPoly=polyfit([1:Spec_av.sz(Spec_av.dims.averages)]',phs(:,1),1);

    iter=iter+1;
end

% For water signal, perform spectral registration in the time
% domain. This is generally preffered for water, but need to check higher
% b-values to be sure.
Specw_av=op_alignAverages(CCw,5*tmax,'n');
    
% Average spectra, and leftshift if needed.
Spec_av = op_leftshift(op_averaging(Spec_av),Spec_av.pointsToLeftshift);
Specw_av = op_leftshift(op_averaging(Specw_av),Specw_av.pointsToLeftshift);


keyboard%#######
% Perform eddy current correction using the water acquisition. Then
% frequency and phase shift based on the creatine peak of the ECC
% spectrum
Spec_av_ECC = FnPh_m(ECC_fids(Spec_av, Specw_av));
    
% Perform frequency and phase alignment using the 3PPM Creatine peak,
% Without doing ECC.
Specw_av = FnPh_w(Specw_av);
    
% Perform frequency and phase alignement using the water peak for the
% water unsuppressed acquisition.
Spec_av = FnPh_m(Spec_av);

OUT.Spec_av_ECC = Spec_av_ECC;
OUT.Spec_av = Spec_av;
OUT.Specw_av = Specw_av;

    
keyboard%#######



end

%% Rearranges data to sepparate averages from diffusion conditions
function OUT = RemapSubSpec(IN,BCond)
    Nav = IN.sz(2)/BCond;
    OUT = IN;
    OUT.sz = [IN.sz(1),Nav,BCond];
    OUT.fids=single(zeros(OUT.sz));
    OUT.specs=single(zeros(OUT.sz));
    OUT.averages = Nav;OUT.subspecs = BCond;
    OUT.dims.subSpecs = 3;

    % This could probably be optimised...
    Map_Av = repelem(1:Nav,BCond);Map_B = repmat(1:BCond,1,Nav);
    for JJ = 1:IN.sz(2)
        OUT.specs(:,Map_Av(JJ),Map_B(JJ)) = IN.specs(:,JJ);
        OUT.fids(:,Map_Av(JJ),Map_B(JJ)) = IN.fids(:,JJ);
    end
end

%% Perform eddy current corection using water unsuppressed fids.
function[Out] = ECC_fids(data, water_data)
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

%% Perform automatic phi0 and frequency allignment for water unsuppressed data
function[out] = FnPh_w(Spec)
out = Spec;
    for JJ=1:Spec.sz(Spec.dims.subSpecs)%Loop over subspecs
        % Zero order phase shift using water unsuppressed peak
        Tempw_zp=op_zeropad(out,16);
        indexw=find(abs(Tempw_zp.specs(:,JJ))==max(abs(Tempw_zp.specs((Tempw_zp.ppm>4 & Tempw_zp.ppm<5.4), JJ))));%FindMax for sub spec
        ph0w=-phase(Tempw_zp.specs(indexw,JJ))*180/pi;%FindPhase at this point
        
        % Apply phase to subspecs, with and without zero-padding
        outw_ph=op_addphaseSubspec(out,ph0w,JJ); 
        Tempw_zp=op_addphaseSubspec(Tempw_zp,ph0w,JJ);

        [~,frqShiftw]=op_ppmref(Tempw_zp,4,4.5,4.65,JJ); %Find f-shift using zero-padded data
        out=op_freqshiftSubspec(outw_ph,frqShiftw); %Apply shift to data
    end
end

%% Automatic phi0 correction for metabolite spectra using Creatine Peak.
%  (May not be optimal in all cases)
function[out] = FnPh_m(Spec)
    %Currently looping over all subspecs to apply independent ph0 and f
    %correction. Might want to spend a bit more effort correlating between
    %subspecs at some point.
    
    for JJ=1:Spec.sz(Spec.dims.subSpecs)%Loop over subspecs
        %Temp = op_takesubspec(Spec,JJ);
        Temp_zp=op_zeropad(Spec,16);
        
        [Temp_ph,ph0]=op_autophase(Spec,2.8,3.2,0,JJ);
        Temp_zp=op_addphaseSubspec(Temp_zp,ph0,JJ);
    
        % Frequency shift all spectra so that Creatine appears at 3.027 ppm:
        [~,frqShift]=op_ppmref(Temp_zp,2.8,3.2,3.027,JJ);
        out=op_freqshift(Temp_ph,frqShift);% Final metabolite spectra
    end
end



