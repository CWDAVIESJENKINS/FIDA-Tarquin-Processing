function[out , out_noproc] = Tarquin_ReadDat(InName , OutName , N , ph0 , ph1 )
% V 0.6 (out.p.Version.ReadDat)
% Function for conversion of Siemens.dat to JMRUI.txt format for Tarquin
% analysis. This code is a modified version of an FID-A io example script.
% Input:    InName = Path to input .dat file.
% Optional: OutName = path to store .txt file. (If excluded, then data stored with .dat)
%           N = dimension of spectrum of interest. (For extracting sub-acquisitions of MPRESS, for example)
%           ph0 = Zero order phase correction
%           ph1 = First order phase correction
% Output:   OutName = Path to output .txt file.

%% Initalisation & Hacking
iterin=20;
tmaxin=0.2;
RawData=io_loadspec_twix(InName);

if ~exist('OutName','var')
    OutName = strrep(InName,'.dat','.txt');
end
RawData.p.RawDataLocation = InName; % Add data locations to structure
RawData.p.PreprocessedLocation = OutName;

%keyboard

if ~exist('N','var')
    N=1; % N = Which subspectrum to process. (N=1 for MEGA-sLASER 'off')
end
if ~exist('ph0','var')
    ph0 = 0;
end
if ~exist('ph1','var')
    ph1 = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RawData.fids = RawData.fids(:,:,:,N);   %List{9}
RawData.specs = RawData.specs(:,:,:,N); %
RawData.subspecs = 1;                   %
RawData.rawSubspecs= 1;                 %  
RawData.sz(4) = 1;                      %
RawData.dims.subSpecs = 0;              %
%RawData.averages = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read and combine data (automated)

CoilCombos=op_getcoilcombos(op_averaging(RawData),RawData.pointsToLeftshift+1,'h'); % Finds the relative coil phases and amplitudes, and weights the coils based on max signal.

[out_cc]=op_addrcvrs(RawData , RawData.pointsToLeftshift+1 , 'h' , CoilCombos); % Perform weighted coil recombination
out_noproc=op_leftshift(op_averaging(out_cc),out_cc.pointsToLeftshift); % Remove leading datapoints from the fid to get rid of 1st order phase errors.
out_rm = out_cc;

%% Bad average removal

% Check if user wants bad average removal
out_cc2=out_cc;
nBadAvgTotal=0;
nbadAverages=1;

nsd=3;
iter=1;
nbadAverages=1;
nBadAvgTotal=0;
out_cc2=out_cc;
while nbadAverages>0;
    [out_rm,metric{iter},badAverages]=op_rmbadaverages(out_cc2,nsd,'t'); % Bad averages are identified by subtracting from the median, and then calculating the RMS. If greater than 'nsd' specified, then it is discarded.
    nbadAverages=length(badAverages);

    nBadAvgTotal=nBadAvgTotal+nbadAverages

    out_cc2=out_rm;
    iter=iter+1;
end
out_cc2.nBadAvg = nBadAvgTotal;
%% Frequency correction

out_rm2= out_rm ; 
fsPoly=100; 
phsPoly=1000;
fsCum=zeros(out_rm2.sz(out_rm2.dims.averages),1);
phsCum=zeros(out_rm2.sz(out_rm2.dims.averages),1);
iter=1;
while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
    close all
    tmax=tmaxin+0.04*randn(1);
    fmin=1.8+0.04*randn(1);
    fmaxarray=[3.5,3.5,4.2,4.2,7.5,4.2];
    fmax=fmaxarray(randi(6,1));

    %op_RobustSpecReg(out_rm2);
    
    [out_aa,fs,phs]=op_alignAverages_fd(out_rm2,fmin,fmax,tmax,'n'); % Perform spectral registration in the time domain to correct frequency and phase drifts.
    fsCum=fsCum+fs(:,1);
    phsCum=phsCum+phs(:,1);

    fsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs(:,1),1);
    phsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs(:,1),1);

    out_rm2=out_aa;

    iter=iter+1;
end

out=op_leftshift(op_averaging(out_aa),out_aa.pointsToLeftshift);
%% Phase correction

% Add phases
out=op_addphase(out,ph0,ph1);
out_noproc=op_addphase(out_noproc,ph0,ph1);

%% Write data

% Write resulting spectrum to filename specified
if exist('OutName','var')
    io_writejmrui(out , OutName); % For JMRUI format
    %io_writelcm(out , OutName , out.te); % For LCM format
    
end

out.p.Version.ReadDat = '0.6';

end