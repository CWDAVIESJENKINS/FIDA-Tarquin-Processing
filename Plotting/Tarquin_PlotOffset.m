function [] = Tarquin_PlotOffset(FitCell,Xlim,Offset)
% function [] = Tarquin_PlotOffset(FitCell,Xlim,Offset)
% Function to plot series of Tarquin results
% Input:     FitCell=Cell array of fit results
%            Xlim=Plot bounds: [Lower, Upper];
%            Offset=Y-offset between spectra

if ~exist('Xlim','var')
    Xlim = [0.5 4];
end

for JJ=1:length(FitCell)
    hold on
    %Ind = find(FitCell{JJ}.Plot.ppm>Xlim(1) & FitCell{JJ}.Plot.ppm<Xlim(2));
    plot(FitCell{JJ}.Plot.ppm, FitCell{JJ}.Plot.Data + JJ*Offset, 'k');
    plot(FitCell{JJ}.Plot.ppm, FitCell{JJ}.Plot.Fit+FitCell{JJ}.Plot.Baseline + JJ*Offset, 'r');
end
xlim(Xlim)
grid on;xlabel('Frequency (PPM)');ylabel('Spectrum (A.u.)')
set(gca,'XDir','reverse')
end