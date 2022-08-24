function[] = Tarquin_Plot(Struct , Opt)
% V0.6
% Function for plotting of a Tarquin fit. Plots data against Fit+Baseline.
% Also note: sum(MetaboliteProfiles)+Baseline = Fit.
% Input:    Struct = Tarquin fit in my matlab format.
% Optional  Opt = cell array of metabolite names to include in plot. e.g. 
%                 {'NAA','GABA','Cr'} will plot NAA, creatine, and GABA 
%                 along with fit, data, and residual. Can also add 'res' or
%                 'Res' to plot residual
% V0.6 - Rework with new structures, and made residual plot optional
if ~exist('Opt','var') % empty cell if no metabolites specified
    Opt = {};
elseif isstr(Opt)
    Opt = {Opt};
end

% Make a bif figure window
F = figure;
F.Position = [1,1,1133,920]; 
xlabel('Freq(PPM)');ylabel('Spectrum(a.u.)')

plot(Struct.Plot.ppm , Struct.Plot.Data , 'k') % Plot Data
hold on
plot(Struct.Plot.ppm , Struct.Plot.Fit + Struct.Plot.Baseline , 'r') % Plot Tarquin fit
%
Leg = {'Data','Fit'};

for J=1:length(Opt) % Loop over metabolites to add them to plot
    if (strcmp(Opt{J},'Res') || strcmp(Opt{J},'res'))
        plot(Struct.Plot.ppm , Struct.Plot.Data - (Struct.Plot.Fit + Struct.Plot.Baseline) , 'g') % Plot Residue
        Leg{J+2} = Opt{J};
    else
        plot(Struct.Plot.ppm , Struct.Plot.Metab.(Opt{J}))
        Leg{J+2} = Opt{J};
    end
end

grid('minor');xlim([1 4.5]);legend(Leg);set(gca, 'xdir' , 'reverse');

end