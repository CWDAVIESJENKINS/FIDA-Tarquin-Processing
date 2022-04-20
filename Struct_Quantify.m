function TarquinFit = Struct_Quantify(MRS_struct, TarquinFit)
%Should probably add an output summary figure, similar to Gannet.

Mode  = 'Approximate';
%Mode  = 'Explicit';


% 3T Constants
% From Wansapura et al. 1999 (JMRI)
%        T1          T2
% WM   832 +/- 10  79.2 +/- 0.6
% GM  1331 +/- 13  110 +/- 2
%
% From Lu et al. 2005 (JMRI)
% CSF T1 = 3817 +/- 424msec - but state may underestimated and that 4300ms
% is likely more accurate - but the reference is to an ISMRM 2001 abstract
% MacKay (last author) 2006 ISMRM abstract has T1 CSF = 3300 ms
% CSF T2 = 503.0 +/- 64.3 Piechnik MRM 2009; 61: 579
% However, other values from Stanisz et al:
% CPMG for T2, IR for T1
% T2GM = 99 +/ 7, lit: 71+/- 10 (27)
% T1GM = 1820 +/- 114, lit 1470 +/- 50 (29)
% T2WM = 69 +/-3 lit 56 +/- 4 (27)
% T1WM = 1084 +/- 45 lit 1110 +/- 45 (29)

% CWJ Added parameters for 7T data:
% Wyss et al 2013 
% 'Relaxation parameter mapping adapted for 7T and 
% validation against optimized single voxel MRS'

%CWJ Added if statement to correctly choose relaxation parameters for 3T
%    and 7T GABA MRS.
if MRS_struct.Bo == 6.9800
    T1w_WM    = 1.285;
    T2w_WM    = 0.037;
    T1w_GM    = 1.954;
    T2w_GM    = 0.045;
    T1w_CSF   = 3.867;
    T2w_CSF   = 0.311;
elseif MRS_struct.Bo == 2.89362000000000
    T1w_WM    = 0.832;
    T2w_WM    = 0.0792;
    T1w_GM    = 1.331;
    T2w_GM    = 0.110;
    T1w_CSF   = 3.817;
    T2w_CSF   = 0.503;
end

% Determine concentration of water in GM, WM and CSF
% Gasparovic et al. 2006 (MRM) uses relative densities, ref to
% Ernst et al. 1993 (JMR)
% fGM = 0.78
% fWM = 0.65
% fCSH = 0.97
% such that
% concw_GM = 0.78 * 55.51 mol/kg = 43.30
% concw_WM = 0.65 * 55.51 mol/kg = 36.08
% concw_CSF = 0.97 * 55.51 mol/kg = 53.84

concw_GM  = 43.30*1e3;
concw_WM  = 36.08*1e3;
concw_CSF = 53.84*1e3;

TR = MRS_struct.tr/1000;
TE = MRS_struct.te/1000;
TR_water = TR;
TE_water = TE;


meanGMfra = mean(MRS_struct.T1Reg.tissue.GMfra); % average GM fraction across subjects
meanWMfra = mean(MRS_struct.T1Reg.tissue.WMfra); % average WM fraction across subjects

fracGM  = MRS_struct.T1Reg.tissue.GMfra;
fracWM  = MRS_struct.T1Reg.tissue.WMfra;
fracCSF = MRS_struct.T1Reg.tissue.CSFfra;


% calculate GM/WM portion of voxel
tissuefra = fracGM + fracWM;

Fn = fieldnames(MRS_struct.Fit.Metabolite);
for J=1:length(Fn)
    switch Mode
        case 'Approximate'
            if MRS_struct.Bo == 6.9800
                % Proton T1 relaxation times of metabolites in human occipital 
                % white and gray matter at 7 T, Lijing Xin et al. Mean value
                T1_Metab = 1.426; 

                % Localized 1H NMR spectroscopy in differentregions of human brain
                % in vivo at 7T: T2 relaxation times and concentrations of cerebral
                % metabolites, Małgorzata Marjańska et al. Mean value
                T2_Metab = 0.12175;
            elseif MRS_struct.Bo == 2.89362000000000;
                % Relaxation rates approximated as same rate as Glx
                T1_Metab = 1.23; % Posse et al. 2007 (MRM)
                T2_Metab = 0.18; % Ganji et al. 2012 (NMR Biomed)
            end
        case 'Explicit'
            error('Chris needs to sort this out...')
            %AIM:
            % Have 3T and 7T estimates for T1 and TE for major metabs and
            % use this to correct data.
            
    end

    %MRS_struct.Fit.Metabolite.(Fn{J}).SignalAmplitude
    
    TissCorr = ...
        (fracGM * concw_GM * (1-exp(-TR_water/T1w_GM)) * (exp(-TE_water/T2w_GM)) / ((1-exp(-TR/T1_Metab)) * (exp(-TE/T2_Metab))) + ...
        fracWM * concw_WM * (1-exp(-TR_water/T1w_WM)) * (exp(-TE_water/T2w_WM)) / ((1-exp(-TR/T1_Metab)) * (exp(-TE/T2_Metab))) + ...
        fracCSF * concw_CSF * (1-exp(-TR_water/T1w_CSF)) * (exp(-TE_water/T2w_CSF)) / ((1-exp(-TR/T1_Metab)) * (exp(-TE/T2_Metab))));
    % Apply tissue corection based on known concnetrations of water in GM,
    % WM, and CSF. Also adjust for relative differences in signals as a
    % function of T1 & T2
    TarquinFit.Metabolite.(Fn{J}).SignalAmplitude_CSF = TarquinFit.Metabolite.(Fn{J}).SignalAmplitude / tissuefra;
    TarquinFit.Metabolite.(Fn{J}).SignalAmplitude_TC = TarquinFit.Metabolite.(Fn{J}).SignalAmplitude*TissCorr;
    TarquinFit.Metabolite.(Fn{J}).SignalAmplitude_CSF_TC = TarquinFit.Metabolite.(Fn{J}).SignalAmplitude_CSF*TissCorr;
    TarquinFit.Metabolite.(Fn{J}).CRLB_TC = TarquinFit.Metabolite.(Fn{J}).CRLB*TissCorr; % Scale CRLB with Tissue Correction factor.
end

end