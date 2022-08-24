function[OUT , List]=DWMRS_ReadDat_v2b(IN , OutDir)
%function[Spec_Av, Specw_Av , Spec_Av_ECC]=DWMRS_ReadDat_v2(IN , OutDir)
% Read in raw dMRS data, perform coil coombination, seperate averages by
% diffussion condition, perform ECC, f & ph correction, and export spectra
% as LCModel .txt format for analysis.
% Input:    IN = cellarray containing {a,b}, where a is the location of the
%                metabolite acquisition, and b is the corresponding
%                reference
%           OutDir = Directory to save the averaged and correceted spectra
%                    for analysis.
% Output    OUT = structure containing processed metabolite data
%                 with/without ECC, and corresponding water references.

%% Initialise some parameters and load data

iterin=10; tmaxin=0.25; BadAv_nsd=2;

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
Ndir = RawData.p.alFree{4};Ng = RawData.p.alFree{5};RepSize = Ndir*Ng+1;
Nav = RawData.sz(3)/RepSize;Nav_ref = RawData_ref.sz(3)/RepSize;

RawData.p.RawDataLocation = file_Met; % Add data locations to structure
RawData_ref.p.RawDataLocation = file_Ref;

%% Get relative coil phases & amplitudes using B0 unsuppressed data
% Coils phased using 1st point in time domain, and weighting based on S/N^2.
% Hall et al. Neuroimage 2014
coilcombos = op_getcoilcombos(RawData_ref,1,'h');
CC = op_addrcvrs(RawData,1,'h',coilcombos);
CCw = op_addrcvrs(RawData_ref,1,'h',coilcombos);

%% Sepperate out by diffusion condition
Counter = 1;
for J=RepSize:-1:1 % Loop over and seperate each diffusion codition
    Ind_Met = (RepSize:RepSize:CC.sz(2))-J+1;
    Ind_Ref = (RepSize:RepSize:CCw.sz(2))-J+1;
    Spec_cc{Counter} = Subset(CC,Ind_Met);
    Specw_cc{Counter} = Subset(CCw,Ind_Ref);
    Counter = Counter+1;
end

%% Bad average removal

% May need to add lb for this at some point

% This is done by subtracting each average from the median, and then 
% calculating the root mean squared of this difference spectrum

for J=1:RepSize
    Rep = Spec_cc{J};
    iter=1;nbadAverages=1;nBadAvgTotal=0;
    while nbadAverages>0;
        [out_rm,metric,badAverages]=op_rmbadaverages(Rep,BadAv_nsd,'t',1);
        nbadAverages=length(badAverages);
        nBadAvgTotal=nBadAvgTotal+nbadAverages;
        Rep=out_rm;
        iter=iter+1;
    end
    Spec_cc2{J}=Rep;
    Spec_cc2{J}.p.NBadAv=nBadAvgTotal;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Repw = Specw_cc{J};
    
    iter=1;nbadAverages=1;nBadAvgTotal=0;
    if Nav_ref>=4
        while nbadAverages>0;
            [out_rm,metricw,badAverages]=op_rmbadaverages(Repw,BadAv_nsd,'t',1);
            nbadAverages=length(badAverages);
            nBadAvgTotal=nBadAvgTotal+nbadAverages;
            Repw=out_rm;
            iter=iter+1;
        end
        Specw_cc2{J}=Repw;
        Specw_cc2{J}.p.NBadAv=nBadAvgTotal;
    else
        Specw_cc2{J}=Repw;
    end
end

%% Averaging, ECC, and freq/phase correction
Spec_Av=cell(1,RepSize);Spec_Av_ECC=cell(1,RepSize);Specw_Av=cell(1,RepSize); %Pre-alocate cells
OutMat=zeros(RepSize,CC.sz(1));
for J=1:RepSize % Again loop over diffusion conditions
    
    Rep = Spec_cc2{J};
    %[Rep,Metric{J},badAverages{J}]=op_rmbadaverages(Spec_cc{J},3,'t');
    Repw = Specw_cc2{J};
    
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
    
    % For water signal, perform spectral registration in the time
    % domain. This is generally preffered for water, but need to check higher
    % b-values to be sure.
    if Nav_ref>1
        Repw=op_alignAverages(Repw,5*tmax,'n');
        Specw_Av{J} = op_leftshift(op_averaging(Repw),Repw.pointsToLeftshift);
    else
        Specw_Av{J} = op_leftshift(Repw,Repw.pointsToLeftshift);
        Specw_Av{J}.flags.averaged = 1;
        Specw_Av{J}.dims.averages = 0;
    end

    Cumf{J} = fsCum;
    Cumph{J} = phsCum;
    
    % Average spectra, and leftshift if needed.
    Spec_Av{J} = op_leftshift(op_averaging(Rep),Rep.pointsToLeftshift);
    
    % Perform eddy current correction using the water acquisition. Then
    % frequency and phase shift based on the creatine peak of the ECC
    % spectrum
    Spec_Av_ECC{J} = FnPh_m(ECC_fids(Spec_Av{J}, Specw_Av{J}));
    
    % Perform frequency and phase alignment using the 3PPM Creatine peak,
    % Without doing ECC.
    Specw_Av{J} = FnPh_w(Specw_Av{J});
    
    % Perform frequency and phase alignement using the water peak for the
    % water unsuppressed acquisition.
    Spec_Av{J} = FnPh_m(Spec_Av{J});
    
    OUT.Spec_Av = Spec_Av;
    OUT.Specw_Av = Specw_Av;
    OUT.Spec_Av_ECC = Spec_Av_ECC;
    OutMat(J,:)=real(Spec_Av_ECC{J}.specs);
end

%% Create a series of summary plots

% Top-down image across spectral range
Fig_TopDown = figure;
%Fig_TopDown.Position = [1,1,1133,920]; 
I=Spec_Av_ECC{1}.ppm>1.5 & Spec_Av_ECC{1}.ppm<4;
imagesc(Spec_Av_ECC{1}.ppm(I), 1:RepSize, OutMat(:,I));set(gca,'XDir','reverse');
colorbar;colormap hot;xlabel('Frequency (PPM)');ylabel('Diffusion condition');
%
Fig_Combined = figure;
%Fig_Combined.Position = [1,1,1133,920]; 
subplot2(3,1,1);plot(RawData.ppm,real(RawData.specs(:,:,1)));xlim([1.5 4]);set(gca, 'XTickLabel', []);title('Raw data');grid on;set(gca,'XDir','reverse');
subplot2(3,1,2);plot(RawData.ppm,real(CC.specs));xlim([1.5 4]);set(gca, 'XTickLabel', []);title('Coil combined data');grid on;set(gca,'XDir','reverse');
subplot2(3,1,3);plot(RawData.ppm,OutMat');xlim([1.5 4]);title('Final ECC data');grid on;set(gca,'XDir','reverse');
xlabel('Frequency (PPM)');

%% If Out directory specified, then save combined spectra in repsective directories & sumary plots
if exist('OutDir','var')
    Home=pwd;cd(OutDir)
    
    CJ_SavFig(Fig_TopDown,'ECC_TopDown',OutDir);
    CJ_SavFig(Fig_Combined,'CoilCombination',OutDir);
    
    mkdir('Met');mkdir('Ref');mkdir('Met_ECC')
    for J=1:RepSize
        cd([OutDir,'/Met'])
        io_writelcm(Spec_Av{J} , sprintf('Spec_%03i.txt',J) , Spec_Av{J}.te);
        cd([OutDir,'/Ref'])
        io_writelcm(Specw_Av{J} , sprintf('Spec_%03i.txt',J) , Spec_Av{J}.te);
        cd([OutDir,'/Met_ECC'])
        io_writelcm(Spec_Av_ECC{J} , sprintf('Spec_%03i.txt',J) , Spec_Av{J}.te);
    end
    cd(Home)

    List.Met = FindFiles([OutDir,'/Met/*.txt']);List.Met=List.Met(:,1);
    List.Met_ECC = FindFiles([OutDir,'/Met_ECC/*.txt']);List.Met=List.Met(:,1);
    List.Ref = FindFiles([OutDir,'/Ref/*.txt']);List.Ref=List.Ref{1,1};
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

function[Out] = ECC_fids2(data, water_data)
% Uses Itamar's svd code with hard coded TdFrac for ECC of low SNR data.

TdFrac = 4;
[td_synth, ~, ~]=svdfid(water_data.fids, TdFrac, water_data.spectralwidth , -water_data.txfrq/2e6, water_data.txfrq/2e6, 0, 200, 15, '' );

Kphase=unwrap(angle(td_synth));
fids=zeros(size(data.fids));

K_data=abs(data.fids);
Kphase_data=unwrap(angle(data.fids));
Kphase_corr=(Kphase_data)-Kphase;
fids=K_data.* exp(1i*Kphase_corr);

Out = data;
Out.fids=double(fids);
Out.specs = fftshift(ifft(Out.fids,[],Out.dims.t),Out.dims.t);

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









%%
% figure('position',[0 0 560 420]);
% subplot(1,2,1);
% plot(out_cc.ppm,out_cc.specs);xlim([1 5]);xlabel('Frequency (ppm)');ylabel('Amplitude (a.u.)');title('Before removal of bad averages:');
% subplot(1,2,2);
% plot(out_rm.ppm,out_rm.specs);xlim([1 5]);xlabel('Frequency (ppm)');ylabel('Amplitude (a.u.)');title('After removal of bad averages:');
% 
% figure('position',[0 500 560 420]);
% plot([1:length(metric{1})],metric{1},'.r',[1:length(metric{iter-1})],metric{iter-1},'x');xlabel('Scan Number');ylabel('Unlikeness metric (a.u.)');title('Unlikeness metrics before and after bad averages removal');
% legend('before','after');
% legend boxoff;


