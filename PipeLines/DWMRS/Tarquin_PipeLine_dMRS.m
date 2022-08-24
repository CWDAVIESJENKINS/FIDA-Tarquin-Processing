function[DiffFit] = Tarquin_PipeLine_dMRS(Spec_Fit, b, OutLoc)

%% Extract data
for J=1:length(Spec_Fit)
    % Grab neuronal marker
    Ind = find(strcmp(Spec_Fit{J}.Metab.Properties.RowNames,'TNAA'));
    TNAA(J) = Spec_Fit{J}.Metab{Ind,1};
    CRLB_TNAA(J) = Spec_Fit{J}.Metab{Ind,2};
    
    % Grab Glial marker
    Ind = find(strcmp(Spec_Fit{J}.Metab.Properties.RowNames,'TCho'));
    TCho(J) = Spec_Fit{J}.Metab{Ind,1};
    CRLB_TCho(J) = Spec_Fit{J}.Metab{Ind,2};
    
    SNR(J) = Spec_Fit{J}.Diag.SNR;
       
end


%% detect directions using b-values
MaxbLoc = find(b==max(b));
NDir = length(MaxbLoc);

RepSize = (length(Spec_Fit)-1)./NDir;

Ind = [MaxbLoc'-RepSize+1 , MaxbLoc'];S = size(Ind);
% Rows of Ind are limits of b conditions, and columns are directions.


%% Run mono & bi-exponential models on data
Col = {'r','g','b'};
for J=1:S(1) %loop over directions
    DiffFit = MonoBiExp(b(Ind(J,1):Ind(J,2))',TNAA(Ind(J,1):Ind(J,2))', CRLB_TNAA(Ind(J,1):Ind(J,2))');
    errorbar(b(Ind(J,1):Ind(J,2))',TNAA(Ind(J,1):Ind(J,2))',CRLB_TNAA(Ind(J,1):Ind(J,2)),[Col{J},'x']);
    hold on
    plot(DiffFit.Mono.Fit,Col{J});
    D1(J) = DiffFit.Mono.Fit.D1;
    R2(J) = DiffFit.Mono.GOF.rsquare;
end

F = get(gca,'Children');
legend(F([1,3,5]),{'Dir1','Dir 2','Dir 3'});
grid minor;xlabel('Nominal b-value');ylabel('Signal (a.u.)');
set(gca,'YScale','log');title('Mono-exponential fit')
CJ_SavFig(gcf,'Mono_b_TNAA',OutLoc);


%% Make some pretty pictures



end