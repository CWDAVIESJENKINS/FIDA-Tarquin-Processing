function[] = Tarquin_GenPara(DatFile , OutFile)
%V2 WIP
%function[] = Tarquin_GenPara(DatFile , OutFile)
% Function that reads twix header to generate a Tarquin parameter file for
% analysis. 

%% Setup
if ~exist('OutFile','var')
    OutFile = 'UntitledPara.txt';
end

Twix = mapVBVD(DatFile);
OUT='format jmrui_text';

%% Hard-coded fields
%water_eddy {true | false}";
%w_conc NMR visable water concentration (35880)";
%w_att Water attenuation (0.7)";
%lipid_filter      {true | false} remove signals upfield of 1.8ppm in hsvd";
%lipid_filter_freq Set ppm of lipid filter";

%auto_phase        {true | false}";
%max_phi0          the value of phi0 in rads";
%max_phi1          the value of phi1_max/fs/2";

%auto_ref          {true | false}";
%max_dref          the max deviation from ref allowed by auto_ref";
%ref_freq          frequency of a single peak to be used for auto referencing (ppm)";
%ref_signals       {1h_naa_cr_cho_lip | 1h_naa_cho | 1h_naa_cr_cho | 1h_cr_cho | 1h_naa | 1h_cr | 1h_cho | 1h_h2o | 1h_lip | 31p_pcr | 31p_pcr_gammaapt}";

%dyn_freq_corr     {true | false}";
%dref_signals      reference signals for dynamic frequency correction, options as above";

%basis_csv         path to basis (CSV files)";
%basis_xml         path to basis (precompiled XML files)";
%basis_lcm         path to basis (LCModel .basis format)";
%int_basis         {1h_brain | 1h_brain_exmm | 1h_brain_glth | 1h_brain_gly_glth | 1h_brain_gly_cit_glth | 1h_brain_full | 1h_brain_le | 1h_brain_no_pcr | 1h_brain_metab_only | megapress_gaba | braino | 31p_brain}";


AddSTR('water_eddy true');
AddSTR('w_conc 1');
AddSTR('w_att 1');
AddSTR('lipid_filter true');
AddSTR('auto_phase true');
AddSTR('auto_ref true');
AddSTR('ref_signals 1h_cr');


%% Read twix fields

% Sampling frequency in Hz
if iscell(Twix.hdr.MeasYaps.sRXSPEC.alDwellTime)
    Rx = 1e9./Twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1};
elseif isstruct(Twix.hdr.MeasYaps.sRXSPEC.alDwellTime)
    Rx = 1e9./Twix.hdr.MeasYaps.sRXSPEC.alDwellTime(1);
end
AddSTR(sprintf('fs %f',Rx));

% Transmitter frequency
AddSTR(sprintf('ft %f', Twix.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1}.lFrequency));

% Sequence name
Seq = Twix.hdr.Meas.SequenceFileName;
if strfind(Seq, 'eja_svs_mslaser')
    AddSTR('pul_seq mega_press')
elseif strfind(Seq, 'eja_svs_mpress')
    AddSTR('pul_seq mega_press')
elseif strfind(Seq, 'eja_svs_slaser')
    AddSTR('pul_seq slaser')
elseif strfind(Seq, 'eja_svs_press')
    AddSTR('pul_seq press')
else
    warning('Sequence - %s not on list at the moment!',Seq);
end
% To add:
%te1               te1 time in seconds for PRESS sequence";
%tm                tm time in seconds for STEAM sequence";
%acq_delay         acquistion delay time for pulse acquire seq in seconds";
%cpmg_pulses       number of pulses for CPMG sequence";

% echo time in seconds
if iscell(Twix.hdr.MeasYaps.alTE)
    TE = Twix.hdr.MeasYaps.alTE{1}*(1e-6);
elseif isstruct(Twix.hdr.MeasYaps.alTE)
    TE = twix_obj.hdr.MeasYaps.alTE(1)*(1e-6);
end
AddSTR(sprintf('echo %f',TE));








%% Write to file
fid = fopen(OutFile,'wt');
fprintf(OUT);
fclose(fid);

%%
function AddSTR(STR)
    OUT = sprintf('%s\n%s',OUT,STR);
end

end


