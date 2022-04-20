function[] = Tarquin_Checker(Struct)
% function[] = Tarquin_Checker(Struct)
%
% function for checking Tarquin fit. Plots the fit, baseline, residue, and
% major neural metabolites, as well as some useful metabolite
% concentrations.
% V1.0


% Make a bif figure window
F = figure;
F.Position = [1,1,1133,920]; 


% Plot Data after baseline subtraction, and fit
subplot(2,2,1)
plot(Struct.Plot.ppm , Struct.Plot.Data-Struct.Plot.Baseline , 'k') 
hold on
plot(Struct.Plot.ppm , Struct.Plot.Fit , 'r') % Plot Tarquin fit
grid('minor');xlim([1 4.7]);set(gca, 'xdir' , 'reverse');
xlabel('Frequency (ppm)');ylabel('Spectrum (a.u.)')
title('Fit following baseline subtraction')


NormFactor = max(Struct.Plot.Data(Struct.Plot.ppm>1 & Struct.Plot.ppm<4.7));
subplot(2,2,2)
plot(Struct.Plot.ppm , Struct.Plot.Baseline./NormFactor , 'r') 
hold on
plot(Struct.Plot.ppm , (Struct.Plot.Data - (Struct.Plot.Fit + Struct.Plot.Baseline))./NormFactor , 'g') % Plot Residue
grid('minor');xlim([1 4.7]);set(gca, 'xdir' , 'reverse');
xlabel('Frequency (ppm)');ylabel('Normalised spectrum (a.u.)');
title('Baseline and residue normalised to max of spectrum')

subplot(2,2,3)
plot(Struct.Plot.ppm , Struct.Plot.Data-Struct.Plot.Baseline , 'k')
hold on
Leg{1} = 'Data';

Opt = {'NAA','NAAG','Cr','PCr','GPC'};
Col = {'m','r','b','c','g'};
for J=1:length(Opt) % Loop over metabolites to add them to plot
    plot(Struct.Plot.ppm , Struct.Plot.Metab.(Opt{J}),Col{J});
    Leg{J+1} = Opt{J};
end
grid('minor');xlim([1.8 4.2]);set(gca, 'xdir' , 'reverse');
L = legend(Leg);L.Location = 'northwest';
xlabel('Frequency (ppm)');ylabel('Spectrum (a.u.)')

subplot(2,2,4)
axis off
T1 = text(0,0.9,sprintf('Water FWHM = %.2fHz',Struct.Diag.waterFWHM_Hz));
T1.FontSize = 14;
T2 = text(0,0.7,sprintf('SNR = %.2f',Struct.Diag.SNR));
T2.FontSize = 14;
T3 = text(0,0.5,sprintf('Q, std(res)/std(noise) = %.2fHz',Struct.Diag.Q));
T3.FontSize = 14;

end