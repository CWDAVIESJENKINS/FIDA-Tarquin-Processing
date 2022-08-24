Plotting%% Options:

WaterScan = false; % Is this a water-unsuppressed MRS acquisition

Exclude_pts = false; % Whether to exclude points based on SNR or fit quality, as below
SNR_Threshold = 3;      Q_Threshold = 1.5;

Locations.RawData =      {'/cubric/scratch/sapcj3/dMRS_2/Data/20200213_80/meas_MID842_dMRS_gmax80_5grad_3dir_met_FID7908.dat'
'/cubric/scratch/sapcj3/dMRS_2/Data/20200213_80/meas_MID841_dMRS_gmax80_5grad_3dir_ref_FID7907.dat'}; % Raw data location(s) - cell array if reference acquired, string otherwise
Locations.T1 =            '';
Locations.OutDIR =        '/cubric/scratch/sapcj3/dMRS_2/Analysis/20200213_80'; % Where to save outputs
%
Locations.ParameterFile = '/cubric/scratch/sapcj3/dMRS/ParaFile.txt'; %MRS parameter file

%% Initialise some parameters

MatFile = [Locations.OutDIR,filesep,'Data.mat'];
Progress= zeros(1,5);

%% MRS Pre-processing

Locations.SpecOut = [Locations.OutDIR,filesep,'SpecOut'];
mkdir(Locations.SpecOut);
[Spec, QA] = DWMRS_ReadDat_v5(Locations.RawData, Locations.SpecOut, WaterScan);
Progress(1)=1;
save(MatFile,'Spec','QA','Locations','Progress');

%% MRS Fitting
close all
if isfield(Spec.List, 'ECC')
    In = Spec.List.Met_ECC;
else
    In = Spec.List.Met;
end

Locations.SpecFit = [Locations.OutDIR,filesep,'SpecFit'];
mkdir(Locations.SpecFit);

Spec_Fit = cell(1,length(In));
Flag = zeros(1,length(In));SNR = zeros(1,length(In));Q = zeros(1,length(In));
parfor K = 1:length(In)
    Spec_Fit{K} = Tarquin_Run(In{K}, '' , Locations.ParameterFile, '1h_brain', sprintf('--echo %f --te1 %f',Spec.te/1000,Spec.te/2000));
    Flag(K) = Spec_Fit{K}.Flag;
    Q(K) = Spec_Fit{K}.Diag.Q;
    SNR(K) = Spec_Fit{K}.Diag.SNRmetab;
end
QA.Fit.Flag = Flag;
QA.Fit.Q_Error = Q;
QA.Fit.SNR = SNR;

Tarquin_Plot2(Spec_Fit{1});CJ_SavFig(gcf,'b0',Locations.SpecFit);
figure;plot(Spec.Diff_G,QA.Fit.SNR,'kx'); hold on;plot(Spec.Diff_G, ones(1,length(Spec.Diff_G)).*SNR_Threshold,'r')
xlabel('Nominal G');ylabel('SNR');legend('SNR','Threshold');CJ_SavFig(gcf,'SNR',Locations.SpecFit);

clear Flag Q SNR

Progress(2)=1;
save(MatFile,'Spec','QA','Locations','Progress','Spec_Fit');

%% Initial Diffusion Fit
close all
if Exclude_pts
    Excl = (QA.Fit.SNR < SNR_Threshold | QA.Fit.Q_Error < Q_Threshold);
else
    Excl = zeros(length(QA.Fit.SNR));
end

Locations.Diff1 = [Locations.OutDIR,filesep,'Diff1'];
mkdir(Locations.Diff1);
[b,Locations.NonLin] = DWMRS_Parameters(Spec, Locations.Diff1);
DiffFit1 = Tarquin_PipeLine_dMRS(Spec_Fit,b,Locations.Diff1);

Progress(3)=1;
save(MatFile,'Spec','QA','Locations','Progress','Spec_Fit','DiffFit1','b');

%% Non-linearity corrections
close all
if ~strcmp(Locations.T1,'')
    Locations.Diff2 = [Locations.OutDIR,filesep,'Diff2'];
    mkdir(Locations.Diff2);
    Locations.Nii = Gen_nii(Locations.T1, Locations.Diff2);
    Spec = Struct_CoRegister(Spec, Locations.Nii, Locations.Diff2);
    Locations.NonLin.Mask = Spec.T1Reg.mask.outfile;
    
    bc = DWMRS_NonLinearity(Locations,fileparts(Locations.NonLin.Mask));
    
    Progress(4) = 1;
    save(MatFile,'Spec','QA','Locations','Progress','Spec_Fit','DiffFit1','b','bc');
end

%% Final diffusion Fit
close all
if ~strcmp(Locations.T1,'')
    DiffFit2 = Tarquin_PipeLine_dMRS_Dist2(Spec_Fit, bc, Locations.Diff2);
    D1_NAA = [DiffFit2.TNAA.x1(2), DiffFit2.TNAA.x2(2), DiffFit2.TNAA.x3(2)];
    D1_Cho = [DiffFit2.TCho.x1(2), DiffFit2.TCho.x2(2), DiffFit2.TCho.x3(2)];
    Progress(5) = 1;
    save(MatFile,'Spec','QA','Locations','Progress','Spec_Fit','DiffFit1','b','bc','DiffFit2','D1_NAA','D1_Cho');
end

