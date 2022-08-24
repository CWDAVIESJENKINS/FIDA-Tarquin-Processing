function [OUT , Stats , RES] = Tarquin_Stats(FIT , Info , Met , AdditionalFields)
% V0.2
% Function for the automated processing of an array of Tarquin fit results.
% ICCs and CVs are calcualated for all specified metabolites and ratios, and
% an output cell array created. 
% Input:    Fit = Cell vector of fit structure from Tarquin_Read.
%           Info = Cell array with three columns relating to the subject,
%                  session, and region of each Fit cell.
%           Met = Cell array of metabolites to be tabulated. e.g.
%                 {'Glu','GABA/Cr','GABA/NAA'}.
% Opt:      AdditionalFields = String of reference to sub structure of Fit
%                              to be included in output table.
% Output:   OUT = Output table containing the amplitude ratios, and for
%                 single metabolites, the CRLBs. Rows are sorted by
%                 participant, ancd columns by input order.
%           Stats = Structure containing the statistical analysis. ICCs are
%                   calculated for each region, along with CVs for both
%                   participant and region.
%           RES = 3D matrices of metabolite results.

% Read in subject, session, and Region info
Subj = cell2mat(CellGrab(Info , 1));
Sesh = cell2mat(CellGrab(Info , 2));
Region = CellGrab(Info , 3);

% Find the unique entries in each case
U_S = unique(Subj);
U_R = unique(Region);
U_Instance = unique(Sesh);

%Counter for the table rows
Counter1 = 2;

% Create a struct of NaNs to populate the results.
for J_M = 1:length(Met)
    for J = 1:length(U_R)
        RES.(strrep(Met{J_M} , '/' , '_')).(U_R{J}) = NaN( length(U_Instance) , length(U_S));
    end
end

%Loop over unique subjects, regions, and sessions to populate the OUT
%table and RES structure for each MET entry.
for J_S=U_S'
    F_S = find(Subj==U_S(J_S));
    for J_R = 1:length(U_R)
        F_R = find(~cellfun(@isempty,strfind(Region,U_R{J_R})));
        INT = intersect(F_R , F_S);
        for J_Instance=INT'
            OUT{Counter1,1} = Subj(J_Instance);
            OUT{Counter1,2} = Region{J_Instance};
            Counter2 = 3;
            for J_M = 1:length(Met)
                STR = strsplit(Met{J_M} , '/');
                if length(STR) == 1
                    %May need a try/catch?
                    OUT{1,Counter2} = sprintf('%s (Amp)',STR{1});
                    OUT{Counter1,Counter2} = FIT{J_Instance}.(STR{1}).SignalAmplitude;
                    RES.(strrep(Met{J_M} , '/' , '_')).(U_R{J_R})(Sesh(J_Instance) , Subj(J_Instance)) = FIT{J_Instance}.(STR{1}).SignalAmplitude;
                    Counter2 = Counter2 + 1;
                    OUT{1,Counter2} = sprintf('CRLB)');
                    OUT{Counter1,Counter2} = FIT{J_Instance}.(STR{1}).CRLB;
                    Counter2 = Counter2 + 1;
                elseif length(STR) == 2
                    OUT{1,Counter2} = sprintf('%s / %s',STR{1} , STR{2});
                    OUT{Counter1,Counter2} = FIT{J_Instance}.(STR{1}).SignalAmplitude / FIT{J_Instance}.(STR{2}).SignalAmplitude;
                    RES.(strrep(Met{J_M} , '/' , '_')).(U_R{J_R})(Sesh(J_Instance) , Subj(J_Instance)) = FIT{J_Instance}.(STR{1}).SignalAmplitude / FIT{J_Instance}.(STR{2}).SignalAmplitude;
                    Counter2 = Counter2 + 1;
                else
                    error('Too many "/"!! Need better code, or better typing skilz')
                end
                
            end
            if exist('AdditionalFields','var')
                for J_A=1:length(AdditionalFields)
                    STR = strsplit(AdditionalFields{J_A} , '.');
                    OUT{1,Counter2} = sprintf('%s',STR{length(AdditionalFields)});
                    OUT{Counter1,Counter2} = eval(sprintf('FIT{J_Instance}.%s',AdditionalFields{J_A}));
                    Counter2 = Counter2 + 1;
                end
            end
            Counter1 = Counter1 +1 ;
        end
    end
end

% Calculate ICCs for each region, and specified MET
for J_M =1:length(Met)
    S = 'cat(3';
    for J_R = 1:length(U_R)
        S = sprintf('%s , RES.(strrep(Met{J_M} , ''/'' , ''_'')).%s ',S,U_R{J_R});
        Stats.ICC.(strrep(Met{J_M} , '/' , '_')).(U_R{J_R}) = ICC(RES.(strrep(Met{J_M} , '/' , '_')).(U_R{J_R}) , '1-1');
        
    end
    S = sprintf('%s)',S);
    RES.(strrep(Met{J_M} , '/' , '_')).Full = eval(S);
end

% Calculate CVs for regions and participant, for each MET
C = Counter2+1;
for J_M =1:length(Met)
    S=size(RES.(strrep(Met{J_M} , '/' , '_')).Full);
    for J=1:S(2)
        Stats.CV_P.(strrep(Met{J_M} , '/' , '_'))(J) = mean(CV(permute(RES.(strrep(Met{J_M} , '/' , '_')).Full(:,J,:) , [1 3 2])));
    end
    for J=1:S(3)
        Stats.CV_R.(strrep(Met{J_M} , '/' , '_'))(J) = mean(CV(RES.(strrep(Met{J_M} , '/' , '_')).Full(:,:,J)));
    end
end


end