function[OUT] = Tarquin_Read(Loc)
%function[OUT] = Tarquin_Read(Loc)
% V1.0
% Function for the reading and reformatting Tarquin .csv data to matlab my
% own format. Then delete the files if required.
% Input:     Loc = Cell array containing locations of Tarquin fit result
%                  .csv, and if available, fitted profiles .csv
% Optional   Del = true/false. If ture, delete the output files. The
%                  default is to not delete.
% Output:    OUTPUT = Matlab data structure containing the output and
%                     individual metabolite fits

Tarquin_Config

if ischar(Loc) % If single file, assume it's fit results
    Loc={Loc,''};
end

if exist(Loc{1},'file')
%% Reading output 1
    disp('Reading Output Data ...')
    Temp = delimread(Loc{1},',','mixed');Temp = Temp.mixed;   %Read Out.csv into a cell array using delimread

    AmpsAndErr = Temp([3,6],4:end)';AmpsAndErr=reshape(AmpsAndErr(~isnan(cell2mat(AmpsAndErr))),[],2);%Amps and CRLBs
    LBasis = size(AmpsAndErr,1);%length of basis set
    Metabolites = strrep(Temp(2,4:LBasis+3),'-Cr','Cr');%MetaboliteNames
    ShiftAndDamp = Temp(12:end,[6,7]);ShiftAndDamp = reshape(ShiftAndDamp(~isnan(cell2mat(ShiftAndDamp))),[],2);%CS & dampings for single models
    
    %Create output table for metabolites
    OUT.Metab = AmpsAndErr;OUT.Metab(1:size(ShiftAndDamp,1),[3,4]) = ShiftAndDamp;
    OUT.Metab = cell2table(OUT.Metab);
    OUT.Metab.Properties.VariableNames = {'Amplitude','CRLB','Shit_PPM','Damping_Hz'};
    OUT.Metab.Properties.RowNames = Metabolites;
    
    %Populate fit diagnostics
    FitDiag = Temp([8,9],4:end)';
    for J=1:size(FitDiag,1)
        if ischar(FitDiag{J,1})
            OUT.Diag.(CleanUpName(FitDiag{J,1})) = FitDiag{J,2};
        end
    end
    
    %populate fit parameters
    FitParams = Temp(boolean([zeros(11,1);~cellfun(@(x) isnumeric(x), Temp(12:end,1))]),[1,2]);
    for J=1:size(FitParams,1)
        OUT.Param.(CleanUpName(FitParams{J,1})) = FitParams{J,2};
    end
    
    disp('... Done!')
%% Reading output 2    
    if exist(Loc{2},'file')
        disp('Reading Fit Data ...')
        Temp2 = delimread(Loc{2},',','mixed');Temp2 = Temp2.mixed; % Read in individual metabolite fits

        OUT.Plot.ppm = cell2mat(Temp2(3:end,1));
        OUT.Plot.Data = cell2mat(Temp2(3:end,2));
        OUT.Plot.Fit = cell2mat(Temp2(3:end,3));
        OUT.Plot.Baseline = cell2mat(Temp2(3:end,4));

        Metabolites2 = strrep(Temp2(2,5:end),'-Cr','Cr');
        for J=1:length(Metabolites2)
            OUT.Plot.Metab.(Metabolites2{J}) = cell2mat(Temp2(3:end,J+4));
        end
        disp('... Done!')
    else
        warning('No fit profiles found')
    end
 %% Tidy up   
    OUT.Flag=false;
    if Del
        delete(Loc{1},Loc{2}) % Remove .csv after it has been read in (Good pactice to clean up in case Tarquin fails)
    end
else
    OUT.Flag = true; % If output files not there, assume Tarquin failed, and flag. Tarquin_Run will add command line output
    warning('Warning! Outputs not found!')
end

end

function[OutStr] = CleanUpName(InStr) %Remove special characters from Tarquin fields
OutStr = strrep(strrep(strrep(strrep(InStr,' ',''),'(','_'),')',''),'/','p');
end