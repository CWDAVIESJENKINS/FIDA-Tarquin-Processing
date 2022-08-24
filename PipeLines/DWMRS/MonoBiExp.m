function[OUT]=MonoBiExp(X,Y,W)
% Mono a bi exponential fit to data
MonoExp_Func = fittype('f1.*exp(-D1.*x)','dependent',{'y'},'independent',{'x'},'coefficients',{'f1','D1'} );
Mono2Exp_Func = fittype('f1.*exp(-D1.*x)+f2','dependent',{'y'},'independent',{'x'},'coefficients',{'f1','f2','D1'} );
BiExp_Func = fittype('f1.*exp(-D1.*x)+f2.*exp(-D2.*x)','dependent',{'y'},'independent',{'x'},'coefficients',{'f1','D1','f2','D2'} );

%% Define upper and lower bounds for fitting
LBound_f = 0;    UBound_f = Inf;
LBound_D = 0;    UBound_D = Inf;

%% Fit Mono-exp
Opt1 = fitoptions(MonoExp_Func);
Opt1.StartPoint = [max(Y),1e-5];
Opt1.Lower=[LBound_f,LBound_D];Opt1.Upper=[UBound_f,UBound_D];
Opt1.TolFun =1e-10;Opt1.TolX =1e-10;
if exist('W','var')
    Opt1.Weights=1./W;
end
%Opt1.Exclude=[8]
[Mono.Fit , Mono.GOF] = fit(X, Y, MonoExp_Func, Opt1);

%% Fit Mono-exp with additional constant term
Opt2 = fitoptions(Mono2Exp_Func);
Opt2.StartPoint = [max(Y),0,1e-5];
Opt2.Lower=[LBound_f, LBound_f, LBound_D];Opt2.Upper=[UBound_f, UBound_f, UBound_D];
Opt2.TolFun =1e-10;Opt2.TolX =1e-10;
if exist('W','var')
    Opt2.Weights=1./W;
end
%Opt1.Exclude=[8]
[Mono2.Fit , Mono2.GOF] = fit(X, Y, Mono2Exp_Func, Opt2);


%% Fit Bi-exp
Opt3 = fitoptions(BiExp_Func);
Opt3.StartPoint = [Mono.Fit.f1,Mono.Fit.D1,1,1e-5]; % use mono exponential as startpoint for bi-exp fitting
Opt3.Lower=[LBound_f,LBound_D,LBound_f,LBound_D];Opt3.Upper=[UBound_f,UBound_D,UBound_f,UBound_D];
Opt3.TolFun =1e-10;Opt3.TolX =1e-10;
if exist('W','var')
    Opt3.Weights=1./W;
end
%Opt2.Exclude=[8]
[Bi.Fit , Bi.GOF] = fit(X, Y, BiExp_Func, Opt3);

OUT.Mono = Mono;
OUT.Mono2 = Mono2;
OUT.Bi = Bi;
end
%%

% Mono a bi exponential fit to data
%MonoExp_Func = fittype('exp(a.*x+b)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'} );
%BiExp_Func = fittype('exp(a.*x+b)+exp(c.*x+d)','dependent',{'y'},...
%'independent',{'x'},'coefficients',{'a','b','c','d'} );

% Define upper and lower bounds for fitting
%LBound_a = -Inf;    UBound_a = 0;
%LBound_b = -Inf;    UBound_b = Inf;

% Fit Mono-exp
% Opt1 = fitoptions(MonoExp_Func);
% Opt1.StartPoint = [-1e-5,1];
% Opt1.Lower=[LBound_a,LBound_b];
% Opt1.Upper=[UBound_a,UBound_b];
% [Mono.Fit , Mono.GOF] = fit(X, Y, MonoExp_Func, Opt1);
% 
% %% Fit Bi-exp
% Opt2 = fitoptions(BiExp_Func);
% use mono exponential as startpoint for bi-exp fitting
% Opt2.StartPoint = [Mono.Fit.a,Mono.Fit.b,-1e-5,0]; 
% Opt2.Lower=[LBound_a,LBound_b,LBound_a,LBound_b];
% Opt2.Upper=[UBound_a,UBound_b,UBound_a,UBound_b];
% [Bi.Fit , Bi.GOF] = fit(X, Y, BiExp_Func, Opt2);