function[Struct , RUN] = Tarquin_Run(File,File_Ref,OutName,Para,Basis,Args)
%% function[Struct , RUN] = Tarquin_Run(File,File_Ref,OutName,Para,Basis,Args)
%
% V1.0
% Matlab wrapper for Tarquin. Parses commandline arguments using system,
% and loads data into a Matlab friendly format.
%
% Cite:
%   Wilson, Martin, et al. "A constrained least‐squares approach to the 
%   automated quantitation of in vivo 1H magnetic resonance spectroscopy 
%   data." Magnetic resonance in medicine 65.1 (2011): 1-12.
% Download:
%   http://tarquin.sourceforge.net/download.php
% Tarquin user manual:
%   http://tarquin.sourceforge.net/user_guide/tarquin_user_guide.html
%   Matlab wrapper for Tarquin. Calls Tarquin via commandline, and reads in Tarquin fit outputs
%
% 
% Arguments ignored if input as empty string: ''
% Input:     File = Path to spectrum
% Optional:  File_Ref = Water reference spectrum
%            OutName = Path/Name to save Tarquin outputs under
%            Para = Path to parameter file (src/common/console.cpp for 
%                   available options)
%            Basis = internal basis set to use : 1h_brain, 1h_brain_exmm, 
%                    1h_brain_glth, 1h_brain_gly_glth, 
%                    1h_brain_gly_cit_glth, 1h_brain_full, 1h_brain_le, 
%                    1h_brain_no_pcr, 1h_brain_metab_only, megapress_gaba, 
%                    braino, 31p_brain
%            Args = Additional arguments to use of the form: 
%                   "--Opt1 Arg1 --Opt2 Arg2"
% Output:    Struct = Matlab structure including fit statistics, and 
%                     individual metabolite plots
%            RUN = String used to call Tarquin (Mainly for debuging)

Tarquin_Config;

if ~exist(File,'file')
    error('Metabolite File: %s not found!',File)
end

if strcmp(File(end-3:end) , '.dcm') || strcmp(File(end-3:end) , '.IMA') %Siemens
    RUN = sprintf('%s --input %s --format siemens' ,TPath,File);
elseif strcmp(File(end-3:end) , '.txt') %JMRUI .txt format
    RUN = sprintf('%s --input %s --format jmrui_txt' ,TPath,File);
elseif strcmp(File(end-3:end) , '.dat') %siemens TWIX data
    NEW = strrep(File,'.dat','.txt');
    Temp = Tarquin_ReadDat(File,NEW);
    File = NEW;
    if exist('File_Ref','var') && ~strcmp('',File_Ref)
        NEW = strrep(File_Ref,'.dat','.txt');
        Temp.Ref = Tarquin_ReadDat(File_Ref,NEW);
        File_Ref = NEW;
    end
    RUN = sprintf('%s --input %s --format jmrui_txt' ,TPath,File);  % JMRUI.txt format
else
    error('Unexpected file extension')
end
   
if exist('File_Ref','var') && ~strcmp('',File_Ref)
    if ~exist(File_Ref)
        disp('Reference file not found. Proceeding anyway.')
    end
    RUN = sprintf('%s --input_w %s',RUN,File_Ref);              % If reference scan specified, include it
end
if exist('OutName','var') && ~strcmp('',Para)
    Out=[OutName,'_Out.csv'];Fit=[OutName,'_Fit.csv'];
else
    % If not outfile, then generate temp file and delete once in Matlab
    Out=[tempname,'_Out.csv'];Fit=[tempname,'_Fit.csv'];
    Del=1;
end
if exist('Para','var') && ~strcmp('',Para)
    if ~exist(Para,'file')
        warning('Parameter file not found. Proceeding anyway.')
    end
    RUN = sprintf('%s --para_file %s',RUN,Para);                % If parameter file specified, include it
end
if exist('Basis','var') && ~strcmp('',Basis)
    RUN = sprintf('%s --int_basis %s',RUN,Basis);               % If basis specified, use it
end
if exist('Args','var') && ~strcmp('',Args)
    RUN = sprintf('%s %s',RUN,Args);                            % If additional arguments included, add them
end

RUN = sprintf('%s --output_csv %s --output_fit %s',RUN,Out,Fit); % Finally, output .csv in scratch for matlab load

fprintf('%s\tFitting: %s...',datetime,File)
[ ~ , CMDout] = system(RUN); % Run Tarquin for these options
disp('... Done!')

Struct = Tarquin_Read({Out,Fit},Del); % Run coversion from Tarquin .csv to Matlab format

if exist('Temp','var')
    Temp.p.Version.Run = '1.0';
    Temp.Fit = Struct;
    Struct = Temp;
end

Struct.CMD_out = CMDout;

end
