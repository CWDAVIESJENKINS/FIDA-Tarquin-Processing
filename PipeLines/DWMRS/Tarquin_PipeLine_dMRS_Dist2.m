function[Out] = Tarquin_PipeLine_dMRS_Dist2(Spec, FitCell, b, OutLoc)

%% Extract data from FitCell
for J=1:length(FitCell)
    % Grab neuronal marker
    Ind = find(strcmp(FitCell{J}.Metab.Properties.RowNames,'TNAA'));
    TNAA(J) = FitCell{J}.Metab{Ind,1};
    CRLB_TNAA(J) = FitCell{J}.Metab{Ind,2};
    
    % Grab Glial marker
    Ind = find(strcmp(FitCell{J}.Metab.Properties.RowNames,'TCho'));
    TCho(J) = FitCell{J}.Metab{Ind,1};
    CRLB_TCho(J) = FitCell{J}.Metab{Ind,2};
    
    SNR(J) = FitCell{J}.Diag.SNR;
end

%% Detect directions and fill vectors for fitting:
%MaxbLoc = find(b==max(b));
%NDir = length(MaxbLoc);
%RepSize = (length(FitCell)-1)./NDir;
%Ind = [MaxbLoc'-RepSize+1 , MaxbLoc'];S = size(Ind); % Rows of Ind are limits of b conditions, and columns are directions.

S = size(b);
[~,LoopD] = max(S);

has_b0=false;
if any(Spec.Diff_Dir==0)
    has_b0=true;
end

DIR = unique(Spec.Diff_Dir(Spec.Diff_Dir>0));

Func=@(A,xdata) ( xdata*[A(1); A(2)]);
Opt = optimoptions('lsqcurvefit','OptimalityTolerance',1e-100,'FunctionTolerance',1e-100);

bc1_m = median(b,2);
bc1_m_v = [ones(S(1),1),-bc1_m];

A = figure;
violinplot(b', 1000);
title('b-value distributions')
ylabel('b-value')
xlabel('Acquisition')

B = figure;
C = figure;

FuncP=@(A,xdata) exp(xdata*[A(1); A(2)]);
Col = {'r','g','b'};
Leg = {};
for JJ=1:length(DIR)
    Ind = find(Spec.Diff_Dir, DIR(JJ));
    if has_b0
        Ind = [1, Ind];
    end
    S1 = repmat(log(TNAA(Ind)),1,S(LoopD))';
    S2 = repmat(log(TCho(Ind)),1,S(LoopD))';

    b_p = [ones(1,numel(b(Ind,:))) ; -reshape( b(Ind,:) , 1 , numel(b(Ind,:)))]';
    
    [Out.TNAA.(sprintf('x%i',JJ)),Out.TNAA.(sprintf('resnorm%i',JJ))] = lsqcurvefit(Func,[-7;1e-5],double(b_p),S1,[-inf;-inf],[inf;inf],Opt);
    [Out.TCho.(sprintf('x%i',JJ)),Out.TCho.(sprintf('resnorm%i',JJ))] = lsqcurvefit(Func,[-7;1e-5],double(b_p),S2,[-inf;-inf],[inf;inf],Opt);
    
    figure(B)
    Y = FuncP(Out.TNAA.(sprintf('x%i',JJ)),bc1_m_v);
    errorbar( bc1_m(Ind), TNAA(Ind) ,CRLB_TNAA(Ind), CRLB_TNAA(Ind), std(b(Ind,:),0,2), std(b(Ind,:),0,2), [Col{JJ},'x']);
    hold on; plot(bc1_m(Ind) , Y(Ind), Col{JJ});

    figure(C);
    Y = FuncP(Out.TCho.(sprintf('x%i',JJ)),bc1_m_v);
    errorbar( bc1_m(Ind), TCho(Ind) ,CRLB_TCho(Ind), CRLB_TCho(Ind), std(b(Ind,:),0,2), std(b(Ind,:),0,2), [Col{JJ},'x']);
    hold on; plot(bc1_m(Ind) , Y(Ind), Col{JJ});

    Leg = [Leg,sprintf('Dir%i',JJ)];
end
figure(B);
title('NAA diffusivity')
set(gca,'YScale','log')
G = get(gca,'children');
legend(G(1:2:length(DIR)) , Leg);
xlabel('b-value');ylabel('Signal (A.u.)');
figure(C)
title('Choline diffusivity')
set(gca,'YScale','log')
legend(G(1:2:length(DIR)) , Leg);
xlabel('b-value');ylabel('Signal (A.u.)');

if exist('OutLoc','var')
    CJ_SavFig(A,'b_dist',OutLoc);
    CJ_SavFig(B,'TNAA',OutLoc);
    CJ_SavFig(C,'TCho',OutLoc);
end

end