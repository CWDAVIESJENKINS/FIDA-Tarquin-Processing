function[OUT]=run_simSlaserShaped_WAND2(Data)
% n=Data.sz(1); %= number of points in fid/spectrum
% sw=Data.spectralwidth; %= desired spectral width in [Hz]
% Bfield=Data.Bo; %= main magnetic field strength in [T]

n=Data.sz(1); %= number of points in fid/spectrum
sw=Data.spectralwidth; %= desired spectral width in [Hz]
Bfield=Data.Bo; %= main magnetic field strength in [T]
te = Data.te; %sLASER total echo time [ms]
%n=2048; %= number of points in fid/spectrum
%sw=5998.800240; %= desired spectral width in [Hz]
%Bfield=6.98; %= main magnetic field strength in [T]

lw=2; %= linewidth in [Hz]
SpinSystems = load('/home/sapcj3/Tools/FID-A-master/simulationTools/metabolites/spinSystems.mat'); %= spin system definition structure
rfPulse=io_loadRFwaveform('sampleAFPpulse_HS2_R15.RF','inv'); % adiabatic RF pulse shaped waveform 
refTp=10; %= RF pulse duration in [ms]
flipAngle=180; %= flip angle of refocusing pulses [degrees] (Optional.  Default = 180 deg)
centreFreq=3; %= centre frequency of the spectrum in [ppm] (Optional.  Default = 2.3)
thkX=2; %slice thickness of x refocusing pulse [cm]
thkY=2; %slice thickness of y refocusing pulse [cm]
fovX=2; %size of the full simulation Field of View in the x-direction [cm]
fovY=2; %size of the full simulation Field of View in the y-direction [cm]
nX=8; %Number of grid points to simulate in the x-direction
nY=8; %Number of grid points to simulate in the y-direction
 
ph1=[0 0 0 0];  %phase cycling scheme of first refocusing pulse
ph2=[0 0 90 90]; %phase cycling scheme of second refocusing pulse
ph3=[0 0 0 0]; %phase cycling scheme of third refocusing pulse
ph4=[0 90 0 90]; %phase cycling scheme of fourth refocusing pulse
%
AllSys = {'Cr','GABA','Gln','Glu','H2O','NAA','Scyllo'};AllSysL = length(AllSys); % List of spin systems to simulate

% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using PRESS 
%             sequence.

%set up spatial grid
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]

gamma=42577000; %gyromagnetic ratio

%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
rfPulse=rf_resample(rfPulse,100);

%sys=sysRef0ppm

Gx=(rfPulse.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
Gy=(rfPulse.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]

%Initialize structures:
out_posxy_rpc=cell(length(x),length(y),length(ph1));
out_posxy=cell(length(x),length(y));
out=struct([]);

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%if you do not have the parallel computing toolbox, uncomment the first
%for loop and delete "parfor X=1:length(x)"
%for X=1:length(x)


for METAB = 1:AllSysL
    %eval(['sys=sys' AllSys{METAB} ';']);
    sys=SpinSystems.(['sys',AllSys{METAB}]);
    parfor X=1:length(x)
        for Y=1:length(y)
            for m=1:length(ph1)
                disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) '; Y-position ' num2str(Y) ...
                    ' of ' num2str(length(y)) '; Phase cycle position ' num2str(m) ' of ' num2str(length(ph1)) '!!' AllSys{METAB}]);
                out_posxy_rpc{X}{Y}{m}=sim_sLASER_shaped(n,sw,Bfield,lw,sys,te,...
                  rfPulse,refTp,x(X),y(Y),Gx,Gy,ph1(m),ph2(m),ph3(m),ph4(m),centreFreq);
                if m==1 
                     out_posxy{X}{Y}=out_posxy_rpc{X}{Y}{m};
                 elseif m==2 || m==3
                    out_posxy{X}{Y}=op_addScans(out_posxy{X}{Y},out_posxy_rpc{X}{Y}{m},1);
                else
                    out_posxy{X}{Y}=op_addScans(out_posxy{X}{Y},out_posxy_rpc{X}{Y}{m});

                end
            end
                 out=op_addScans(out,out_posxy{X}{Y});
        end
    end
    %For consistent scaling across different shaped simulations, we need to :
    %1.  Scale down by the total number of simulations run (since these were
    %    all added together.
    numSims=(nX*nY*length(ph1));
    out=op_ampScale(out,1/numSims);

    %2.  Scale by the total size of the simulated region, relative to the size
    %    of the voxel.
    voxRatio=(thkX*thkY)/(fovX*fovY);
    out=op_ampScale(out,1/voxRatio);
    
    OUT{METAB} = out;save('Data','OUT')
    RF{METAB} = io_writelcmraw(OUT{METAB},[AllSys{METAB} '.RAW'],AllSys{METAB});
    
    %waitbar(METAB/AllSysL)
end
