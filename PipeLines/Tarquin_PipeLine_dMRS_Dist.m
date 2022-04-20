function[Out] = Tarquin_PipeLine_dMRS_Dist(FitCell, b)

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

%% detect directions using b-values
%MaxbLoc = find(b==max(b));
%NDir = length(MaxbLoc);

%RepSize = (length(FitCell)-1)./NDir;

%Ind = [MaxbLoc'-RepSize+1 , MaxbLoc'];S = size(Ind);
% Rows of Ind are limits of b conditions, and columns are directions.

S = size(b);

[~,LoopD] = max(S);

% temporary solution
Ind1 = [1,2:5];
Ind2 = [1,8:11];
Ind3 = [1,14:17];

% Mono a bi exponential fit to data
MonoExp_Func = fittype('f1.*exp(-D1.*x)','dependent',{'y'},'independent',{'x'},'coefficients',{'f1','D1'} );
LBound_f = 0;    UBound_f = Inf;
LBound_D = 0;    UBound_D = Inf;

Opt1 = fitoptions(MonoExp_Func);
Opt1.StartPoint = [max(TNAA),1e-5];
Opt1.Lower=[LBound_f,LBound_D];Opt1.Upper=[UBound_f,UBound_D];
Opt1.TolFun =1e-10;Opt1.TolX =1e-10;
Opt2 = Opt1;
Opt2.StartPoint = [max(TCho),1e-5];

Opt2.Weights = 1./CRLB_TCho;

D1_NAA = zeros(3,S(2));
D1_Cho = zeros(3,S(2));
r2_NAA = zeros(3,S(2));
r2_Cho = zeros(3,S(2));

warning('off')
X = waitbar(0,'1','Name','Calculating D1 distribution');
for JJ = 1:S(2)
    waitbar(JJ./S(2),X,sprintf('Completed %i/%i',JJ,S(2)));
    
    Opt1.Weights = 1./CRLB_TNAA(Ind1);
    [MonoFit1 , GOF1] = fit(b(Ind1,JJ), TNAA(Ind1)', MonoExp_Func, Opt1);
    D1_NAA(1,JJ) = MonoFit1.D1;
    r2_NAA(1,JJ) = GOF1.rsquare;    
    Opt1.Weights = 1./CRLB_TNAA(Ind2);
    [MonoFit2 , GOF2] = fit(b(Ind2,JJ), TNAA(Ind2)', MonoExp_Func, Opt1);
    D1_NAA(2,JJ) = MonoFit2.D1;
    r2_NAA(2,JJ) = GOF2.rsquare;    
    Opt1.Weights = 1./CRLB_TNAA(Ind3);
    [MonoFit3 , GOF3] = fit(b(Ind3,JJ), TNAA(Ind3)', MonoExp_Func, Opt1);
    D1_NAA(3,JJ) = MonoFit3.D1;
    r2_NAA(3,JJ) = GOF3.rsquare;    
    %
    Opt2.Weights = 1./CRLB_TCho(Ind1);
    [MonoFit1 , GOF1] = fit(b(Ind1,JJ), TCho(Ind1)', MonoExp_Func, Opt2);
    D1_Cho(1,JJ) = MonoFit1.D1;
    r2_Cho(1,JJ) = GOF1.rsquare;    
    Opt2.Weights = 1./CRLB_TCho(Ind2);
    [MonoFit2 , GOF2] = fit(b(Ind2,JJ), TCho(Ind2)', MonoExp_Func, Opt2);
    D1_Cho(2,JJ) = MonoFit2.D1;
    r2_Cho(2,JJ) = GOF2.rsquare;    
    Opt3.Weights = 1./CRLB_TCho(Ind3);
    [MonoFit3 , GOF3] = fit(b(Ind3,JJ), TCho(Ind3)', MonoExp_Func, Opt2);
    D1_Cho(3,JJ) = MonoFit3.D1;
    r2_Cho(3,JJ) = GOF3.rsquare;    
end
warning('on')
%Opt1.Exclude=[8]

% Col = {'r','g','b'};
% for J=1:S(1) %loop over directions
%     DiffFit = MonoBiExp(b(Ind(J,1):Ind(J,2))',TNAA(Ind(J,1):Ind(J,2))', CRLB_TNAA(Ind(J,1):Ind(J,2))');
%     errorbar(b(Ind(J,1):Ind(J,2))',TNAA(Ind(J,1):Ind(J,2))',CRLB_TNAA(Ind(J,1):Ind(J,2)),[Col{J},'x']);
%     hold on
%     plot(DiffFit.Mono.Fit,Col{J});
%     D1(J) = DiffFit.Mono.Fit.D1;
%     R2(J) = DiffFit.Mono.GOF.rsquare;
%     
%     Out1{J} = TNAA(Ind(J,1):Ind(J,2))
% end
%F = get(gca,'Children');
%legend(F([1,3,5]),{'Dir1','Dir 2','Dir 3'});
%grid minor;xlabel('Nominal b-value');ylabel('Signal (a.u.)');
%set(gca,'YScale','log');
Out.D1_NAA = D1_NAA;
Out.D1_Cho = D1_Cho;
Out.r2_NAA = r2_NAA;
Out.r2_Cho = r2_Cho;
end