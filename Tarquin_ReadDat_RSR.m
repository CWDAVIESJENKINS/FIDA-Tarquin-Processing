function[Spec_Out, Spec_Outw] = Tarquin_ReadDat_RSR(InName, N, OutName)
% function[Spec_Out, Spec_Outw] = Tarquin_ReadDat_RSR(InName, N, OutName)
% V0.7 (out.p.Version.ReadDat)
%
% Coil combination and averaging for Siemens Twix data. Uses adapted Gannet
% code for robust spectral registration and weighted averaging of
% particular subspectra.
%
% Input:     InName = Either, the location of the twix file, or cell array
%                of the form {Metab , Ref}.
% Optional   N = Index of sub-spectrum of index. Pulls subspectra from
%                multi-dimensional twix data e.g. MEGA edit-off
%            OutName = If OutName specified, then JMRUI.txt file created
%                for the processed data. Must be cell if Met&Ref used.
%
% V0.7 ADDED  Robust spectral registration and weighted averaging.
%      CHANGE Brought in line with dMRS pre-processing


%% Input handling and read data
if ~exist('N','var')
    N=1; % Define specific sub sepctrum
end

if iscell(InName)
    WaterFile = true;
    InName_ref = InName{2};InName = InName{1};
    RawDataw=io_loadspec_twix(InName_ref);
    RawDataw.p.RawDataLocation = InName_ref; % Add data locations to structure
    RawDataw = ReduceDim(RawDataw,N);
else
    WaterFile = false;
end

RawData=io_loadspec_twix(InName);
RawData.p.RawDataLocation = InName; % Add data locations to structure
RawData = ReduceDim(RawData,N);


%% Coil combination
if WaterFile
    coilcombos = op_getcoilcombos(RawDataw,1,'h');
    CCw = op_addrcvrs(RawDataw,1,'h',coilcombos);
else
    coilcombos = op_getcoilcombos(RawData,1,'h');
end
CC = op_addrcvrs(RawData,1,'h',coilcombos);

%% phase/freq correction & averaging

% Robust spectral registration, weighted averaging, and leftshift if needed
Spec_align = op_Gannet_RSR(CC);
Spec_Av = op_leftshift(op_Gannet_WeightedAverage(Spec_align), CC.pointsToLeftshift);
Spec_Out = FnPh_m(Spec_Av);
if WaterFile
    % time-domain SR for water if there
    tmax=1.25;
    Spec_avw=op_leftshift(op_averaging(op_alignAverages(CCw,5*tmax,'n')),CC.pointsToLeftshift);
    Spec_Outw = FnPh_w(Spec_avw);
else
    Spec_Outw = [];
end

%% Tidy up outputs and generate processed files
figure;imagesc(Spec_align.ppm,[1:Spec_align.averages],abs(Spec_align.specs)');
xlim([1.5 4.5]);xlabel('Frequency(PPM)');ylabel('Nav');title('Alligned averages');

if exist('OutName','var')
    %Alternative format:
    %io_writelcm(OutStructure , OutName , OutStructure.te); % For LCM format
    if WaterFile
        OutNamew = OutName{2};
        OutName = OutName{1};
        io_writejmrui(Spec_Outw , OutNamew); % For JMRUI format
    end
    io_writejmrui(Spec_Out , OutName); % For JMRUI format
end

end
%% Reduce subspec dims
function[IN] = ReduceDim(IN,N)
if IN.dims.extras>0
    warning('Extras Field detected, but not handled')
end

IN.averages = IN.averages/IN.subspecs;
IN.fids = IN.fids(:,:,:,N);
IN.specs = IN.specs(:,:,:,N);
IN.subspecs = 1;
IN.sz = size(IN.fids);
IN.dims.subSpecs = 0;
%IN.averages = 

end

%% Automatic phi0 & fshift correction for metabolite spectra using Creatine Peak.
%  (May not be optimal in all cases)
function[out] = FnPh_m(Spec)
    Temp_zp=op_zeropad(Spec,16);
    [Temp_ph,ph0]=op_autophase(Spec,2.9,3.1);
    Temp_zp=op_addphase(Temp_zp,ph0,0,4.65,1);
    
    % Frequency shift all spectra so that Creatine appears at 3.027 ppm:
    [~,frqShift]=op_ppmref(Temp_zp,2.8,3.1,3.027);
    out=op_freqshift(Temp_ph,frqShift);% Final metabolite spectra
end

%% Perform automatic phi0 and frequency shift for water unsuppressed data
function[out] = FnPh_w(Spec) 
    
    % Zero order phase shift using water unsuppressed peak
    Tempw_zp=op_zeropad(Spec,16);
    indexw=find(abs(Tempw_zp.specs)==max(abs(Tempw_zp.specs(Tempw_zp.ppm>4 & Tempw_zp.ppm<5.4))));%FindMax
    ph0w=-phase(Tempw_zp.specs(indexw))*180/pi;
    outw_ph=op_addphase(Spec,ph0w,0,4.65,1);
    Tempw_zp=op_addphase(Tempw_zp,ph0w,0,4.65,1);
    
    % 
    [~,frqShiftw]=op_ppmref(Tempw_zp,4,5.5,4.65);
    out=op_freqshift(outw_ph,frqShiftw); 
end