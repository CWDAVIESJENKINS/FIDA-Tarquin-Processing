function [out,outw,out_noproc,outw_noproc]=DWMRS_ReadDat(IN)

% Need to re-do this.
% Write Eddy processing loop as a seperate function
% Water drift as an overall function
%1) Coil combos using water phase
%2) Water drift done prior to Averaging
%3) Seperate acquisitions before averaging

% This seems to be the location of the seq/special acquision parameters:
% Test.hdr.MeasYaps.sWipMemBlock.alFree

% water = true;


% DWMRS meeting
% Seperate out DCs here, and grab coil phases using B0 for all acquisitions
% Test this Vs other method
% Probably most efficient to interface this into Itimar's code then, rather
% than applying my own


% Leave out spec reg for now


iterin=10; tmaxin=0.2; aaDomain='f';

file_Met = IN{1}; file_Ref= IN{2};

% Read in both datasets:
RawData=io_loadspec_twix(file_Met);
RawData_ref=io_loadspec_twix(file_Ref);

%Check that met and ref acquisitions use same parameters
if ~(cell2mat(RawData.p.alFree) == cell2mat(RawData_ref.p.alFree))
    warning('Metabolite & reference scans may have different acquisition parameters!!!')
end

%Grab number of directions and diffusion gradients from met
Ndir = RawData.p.alFree{4};Ng = RawData.p.alFree{5};RepSize = Ndir*Ng+1;
Nav = RawData.sz(3)/RepSize;Nav_ref = RawData_ref.sz(3)/RepSize;
RawData.p.RawDataLocation = file_Met; % Add data locations to structure

%first step should be to combine coil channels.  To do this find the coil
%phases from the water unsuppressed data.
keyboard
% Finds the relative coil phases and amplitudes using 1st point in time
% domain, and weighting based on max signal amplitude
coilcombos=op_getcoilcombos(RawData_ref,1,'w');
% Performs coil combination using weightings in previous line
[outw_cc,fidw_pre,specw_pre,phw,sigw]=op_addrcvrs(RawData_ref,1,'w',coilcombos);
% Uses coil amps and phases from reference to combine met acquisition
[out_cc,fid_pre,spec_pre,ph,sig]=op_addrcvrs(RawData,1,'w',coilcombos);

%Simple averaging schemes for plotting routines
[out_av_cc,fid_av_pre,spec_av_pre]=op_addrcvrs(op_averaging(RawData),1,'w',coilcombos);
raw_av=op_averaging(RawData);





% Would normally remove motion corrupted averages here, but not sure it
% will work for high b-value


%out_rm2=out_rm;
%out_rm2 = out_cc;

keyboard
Cnt=1;
for J=RepSize:-1:1 % Loop over each set of b-values and directions
    fsPoly=100; phsPoly=1000;iter=1;
    Ind = (RepSize:RepSize:out_cc.sz(2))-J+1;
    
    out_rm2 = out_cc;
    out_rm2.specs = out_rm2.specs(:,Ind);
    out_rm2.fids = out_rm2.fids(:,Ind);
    out_rm2.sz = size(out_rm2.specs);
    out_rm2.averages = Nav;
    
    fscum=zeros(out_rm2.sz(out_cc.dims.averages),1);
    phscum=zeros(out_rm2.sz(out_cc.dims.averages),1);
    
    while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
        iter=iter+1
        close all
        tmax=0.25+0.03*randn(1);
        ppmmin=1.6+0.1*randn(1);
        ppmmaxarray=[3.5+0.1*randn(1,2),4+0.1*randn(1,3),5.5+0.1*randn(1,1)];
        ppmmax=ppmmaxarray(randi(6,1));

        [out_aa,fs,phs]=op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'y');
        [outw_aa,fs_w,phs_w]=op_alignAverages(outw_cc,5*tmax,'n');
        
        %[out_aa,fs,phs]=op_alignAverages_fd(out_rm2,fmin,fmax,tmax,'n'); % Perform spectral registration in the time domain to correct frequency and phase drifts.

        fsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs,1)
        phsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs,1)
        iter

        fscum=fscum+fs;
        phscum=phscum+phs;

        out_rm2=out_aa;

    end
    out_DC{Cnt} = out_aa;Cnt=Cnt+1;
    outw_DC{Cnt} = outw_aa;Cnt=Cnt+1;
        
end




%% CHECKPOINT
keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%



subplot(1,2,1);
plot(out_rm.ppm,real(out_rm.specs(:,:)));xlim([1 5]);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
title('Before','FontSize',12);
subplot(1,2,2);
plot(out_aa.ppm,real(out_aa.specs(:,:)));xlim([1 5]);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
title('After','FontSize',12);

plot([1:out_aa.sz(out_aa.dims.averages)],fscum,'.-','LineWidth',2);
xlabel('Scan Number','FontSize',10);
ylabel('Frequency Drift [Hz]','FontSize',10);
legend('Frequency Drift','Location','SouthEast');
title('Estimated Freqeuncy Drift','FontSize',12);

plot([1:out_aa.sz(out_aa.dims.averages)],phscum,'.-','LineWidth',2);
xlabel('Scan Number','FontSize',10);
ylabel('Phase Drift [Deg.]','FontSize',10);
legend('Phase Drift','Location','SouthEast');
title('Estimated Phase Drift','FontSize',12);


totalFreqDrift=mean(max(fscum)-min(fscum));
totalPhaseDrift=mean(max(phscum)-min(phscum));
close all
    














    
%now combine the averages averages
out_av=op_averaging(out_aa);
outw_av=op_averaging(outw_aa);



%now leftshift
out_ls=op_leftshift(out_av,out_av.pointsToLeftshift);
outw_ls=op_leftshift(outw_av,outw_av.pointsToLeftshift);


%now do automatic zero-order phase correction (Use Creatine Peak):
out_ls_zp=op_zeropad(out_ls,16);
[out_ph,ph0]=op_autophase(out_ls,2.9,3.1);
out_ls_zp=op_addphase(out_ls_zp,ph0);

%And now for water unsuppressed data (use water peak):
outw_ls_zp=op_zeropad(outw_ls,16);
indexw=find(abs(outw_ls_zp.specs)==max(abs(outw_ls_zp.specs(outw_ls_zp.ppm>4 & outw_ls_zp.ppm<5.5))));
ph0w=-phase(outw_ls_zp.specs(indexw))*180/pi;
outw_ph=op_addphase(outw_ls,ph0w);
outw_ls_zp=op_addphase(outw_ls_zp,ph0w);



%Frequency shift all spectra so that Creatine appears at 3.027 ppm:
[~,frqShift]=op_ppmref(out_ls_zp,2.9,3.1,3.027);
out=op_freqshift(out_ph,frqShift);
out_noproc=op_freqshift(out_noproc,frqShift);
%And now for water unsuppressed data (user water peak and set to 4.65 ppm):
if water
    [~,frqShiftw]=op_ppmref(outw_ls_zp,4,5.5,4.65);
    outw=op_freqshift(outw_ph,frqShiftw);
    outw_noproc=op_freqshift(outw_noproc,frqShiftw);
end

%Make figure to show the final spectrum:
h=figure('visible','off');
plot(out.ppm,out.specs,'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
legend('diff');
legend boxoff;
box off;
title('Result: Final Spectrum','FontSize',12);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',[0 0 20 10]);
saveas(h,[filestring '/report/figs/finalSpecFig'],'jpg');
saveas(h,[filestring '/report/figs/finalSpecFig'],'fig');


close all;

%RF=io_writelcm(out,[filestring '/' filestring '_lcm'],out.te);
%RF=io_writelcm(outw,[filestring '_w/' filestring '_w_lcm'],outw.te);

%write an html report: 
fid=fopen([filestring '/report/report.html'],'w+');
fprintf(fid,'<!DOCTYPE html>');
fprintf(fid,'\n<html>');
logoPath=which('FID-A_LOGO.jpg');
fprintf(fid,'\n<img src= " %s " width="120" height="120"></body>',logoPath);
fprintf(fid,'\n<h1>FID-A Processing Report</h1>');
fprintf(fid,'\n<h2>Processing pipeline applied to PRESS data using run_pressproc_auto.m</h2>');
fprintf(fid,'\n<p>FILENAME: %s/%s/%s </p>',pwd,filestring,file_Met);
fprintf(fid,'\n<p>DATE: %s </p>',date);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of multi-coil combination:</h2>');
fprintf(fid,'\n<img src= " %s/%s/report/figs/coilReconFig.jpg " width="800" height="400"></body>',pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of removal of bad averages:</h2>');
fprintf(fid,'\n<p>Original number of averages: \t%5.6f </p>',RawData.sz(RawData.dims.averages));
fprintf(fid,'\n<p>Number of bad Averages removed:  \t%5.6f </p>',nBadAvgTotal);
fprintf(fid,'\n<p>Number of remaining averages in processed dataset:  \t%5.6f </p>',out_rm.sz(out_rm.dims.averages));
fprintf(fid,'\n<p>Bad Averages Removal Threshold was:  \t%2.2f </p>',nsd);
fprintf(fid,'\n<img src= " %s/%s/report/figs/rmBadAvg_prePostFig.jpg " width="800" height="600"><img src= " %s/%s/report/figs/rmBadAvg_scatterFig.jpg " width="800" height="400">',pwd,filestring,pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of spectral registration:</h2>');
fprintf(fid,'\n<p>Total frequency drift was: \t%5.6f </p>',max(totalFreqDrift));
fprintf(fid,'\n<p>Total phase drift was: \t%5.6f </p>',max(totalPhaseDrift));
fprintf(fid,'\n<img src= " %s/%s/report/figs/alignAvgs_prePostFig.jpg " width="800" height="600">',pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n<img src= " %s/%s/report/figs/freqDriftFig.jpg " width="400" height="400"><img src="%s/%s/report/figs/phaseDriftFig.jpg " width="400" height="400">',pwd,filestring,pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Final Result:</h2>');
fprintf(fid,'\n<img src= " %s/%s/report/figs/finalSpecFig.jpg " width="800" height="400">',pwd,filestring);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



