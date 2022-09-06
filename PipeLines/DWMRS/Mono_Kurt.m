function[]=Mono_Kurt(Spec_Fit,bc,OutLoc,Excl)
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

%% Remove exclusions if specified
if exist('Excl','var')
    TNAA(Excl) = [];CRLB_TNAA(Excl) = [];
    TCho(Excl) = [];CRLB_TCho(Excl) = [];
    bc(Excl,:) = [];
    Title = 'KurtosisFits_excl';
else
    Title = 'KurtosisFits';
end

%% Plot and save
subplot(1,2,1)
FIT_TNAA = Fitting(bc,TNAA',CRLB_TNAA');
title('TNAA Kurtosis')
subplot(1,2,2)
FIT_TCho = Fitting(bc,TCho',CRLB_TCho');
title('TCho Kurtosis')

CJ_SavFig(gcf,Title,OutLoc)


end



function[OUT]=Fitting(bc,S,Wgt)
b=mean(bc,2);
b_std=std(bc');b_std=b_std(2:end);

Nb=3;
Ind = 1:Nb;


Y = log(S(2:end)./S(1)); % ln(S/S0)
X = double(b(2:end));
W = log(Wgt(2:end)./Wgt(1)); % ln(S/S0)

Ind2 = setdiff(1:length(X),Ind);

LBound_ADC = 0;    UBound_ADC = Inf;
LBound_Kurt = -Inf;   UBound_Kurt = Inf;

%% Mono exp

X1 = X(Ind);
Y1 = Y(Ind);
W1 = W(Ind);

MonoExp_Func = fittype('-x.*ADC','dependent',{'y'},'independent',{'x'},'coefficients',{'ADC'} );
Opt1 = fitoptions(MonoExp_Func);
Opt1.StartPoint = min(Y1)-max(Y1)/max(X1)-min(X1);
Opt1.Lower=LBound_ADC;Opt1.Upper=UBound_ADC;
Opt1.TolFun =1e-10;Opt1.TolX =1e-10;
if exist('W1','var')
    Opt1.Weights=abs(1./W1);
end

[Mono.Fit , Mono.GOF] = fit(X1, Y1, MonoExp_Func, Opt1);


%% Kurtosis

Kurt_Func = fittype('-x.*ADC + K/6.*(x.*ADC)^2','dependent',{'y'},'independent',{'x'},'coefficients',{'ADC','K'} );
Opt2 = fitoptions(Kurt_Func);
Opt2.StartPoint = [min(Y1)-max(Y1)/max(X1)-min(X1),Mono.Fit.ADC];
Opt2.Lower=[LBound_ADC,LBound_Kurt];Opt1.Upper=[UBound_ADC,UBound_Kurt];
Opt2.TolFun =1e-10;Opt1.TolX =1e-10;
if exist('W','var')
    Opt2.Weights=abs(1./W);
end

[Kurt.Fit , Kurt.GOF] = fit(X, Y, Kurt_Func, Opt2);




%% Plot
% A = errorbar(X1,Y1,abs(W1),'rx');hold on
% B = errorbar(X,Y,abs(W),'ko');
A = errorbar(X1,Y1,abs(W1),abs(W1),b_std(Ind),b_std(Ind),'rx');hold on
B = errorbar(X,Y,abs(W),abs(W),b_std,b_std,'ko');
C = plot(Mono.Fit,'r--');
D = plot(Kurt.Fit,'k');

legend([A,B,C,D],'All data','Included in mono-exponential fit','Mono-exponential fit','Kurtosis fit')

xlabel('effective b-value (s/mm^2)')
ylabel('ln(S/S_0)')
grid on

OUT.Mono = Mono;
OUT.Kurt = Kurt;

end