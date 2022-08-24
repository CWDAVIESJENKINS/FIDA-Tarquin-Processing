function[Out] = Tarquin_PipeLine_dMRS_Dist2(FitCell, b, OutLoc)

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

% temporary solution, ignoring top two b-values

L = (length(FitCell)-1)/3;

Ind1 = [1,2:L-2];
Ind2 = [1,L+2:(2*L-2)];
Ind3 = [1,2*L+2:(3*L-2)];
Ind1b = [1,2:L];
Ind2b = [1,L+2:(2*L)];
Ind3b = [1,2*L+2:(3*L)];

%MonoExp_Func = fittype('x*[s0_p;D1]','dependent',{'y'},'independent',{'x'},'coefficients',{'s0_p','D1'} );

%keyboard
S1_p1 = repmat(log(TNAA(Ind1)),1,S(LoopD))';
S1_p2 = repmat(log(TNAA(Ind2)),1,S(LoopD))';
S1_p3 = repmat(log(TNAA(Ind3)),1,S(LoopD))';

S2_p1 = repmat(log(TCho(Ind1)),1,S(LoopD))';
S2_p2 = repmat(log(TCho(Ind2)),1,S(LoopD))';
S2_p3 = repmat(log(TCho(Ind3)),1,S(LoopD))';

b1 = [ones(1,numel(b(Ind1,:))) ; -reshape( b(Ind1,:) , 1 , numel(b(Ind1,:)))]';
b2 = [ones(1,numel(b(Ind2,:))) ; -reshape( b(Ind2,:) , 1 , numel(b(Ind2,:)))]';
b3 = [ones(1,numel(b(Ind3,:))) ; -reshape( b(Ind3,:) , 1 , numel(b(Ind3,:)))]';


Func=@(A,xdata) ( xdata*[A(1); A(2)]);

Opt = optimoptions('lsqcurvefit','OptimalityTolerance',1e-100,'FunctionTolerance',1e-100);

[Out.TNAA.x1,Out.TNAA.resnorm1] = lsqcurvefit(Func,[-7;1e-5],double(b1),S1_p1,[-inf;-inf],[inf;inf],Opt);
[Out.TNAA.x2,Out.TNAA.resnorm2] = lsqcurvefit(Func,[-7;1e-5],double(b2),S1_p2,[-inf;-inf],[inf;inf],Opt);
[Out.TNAA.x3,Out.TNAA.resnorm3] = lsqcurvefit(Func,[-7;1e-5],double(b3),S1_p3,[-inf;-inf],[inf;inf],Opt);

[Out.TCho.x1,Out.TCho.resnorm1] = lsqcurvefit(Func,[-7;1e-5],double(b1),S2_p1,[-inf;-inf],[inf;inf],Opt);
[Out.TCho.x2,Out.TCho.resnorm2] = lsqcurvefit(Func,[-7;1e-5],double(b2),S2_p2,[-inf;-inf],[inf;inf],Opt);
[Out.TCho.x3,Out.TCho.resnorm3] = lsqcurvefit(Func,[-7;1e-5],double(b3),S2_p3,[-inf;-inf],[inf;inf],Opt);

bc1_m = median(b,2);
bc1_m_v = [ones(S(1),1),-bc1_m];

%ERRORBAR(X,Y,YNEG,YPOS,XNEG,XPOS)
%boxplot(b(Ind1,:)', log(TNAA(Ind1)), 'Orientation','horizontal', 'Colors','r'); hold on


A = figure;
violinplot(b', 1000);
title('b-value distributions')
ylabel('b-value')
xlabel('Acquisition')
%
B = figure;
Func=@(A,xdata) exp(xdata*[A(1); A(2)]);
Y = Func(Out.TNAA.x1,bc1_m_v);
errorbar( bc1_m(Ind1b), TNAA(Ind1b) ,CRLB_TNAA(Ind1b), CRLB_TNAA(Ind1b), std(b(Ind1b,:),0,2), std(b(Ind1b,:),0,2), 'rx');
hold on; P1= plot(bc1_m(Ind1b) , Y(Ind1b), 'r');
Y = Func(Out.TNAA.x2,bc1_m_v);
errorbar( bc1_m(Ind2b), TNAA(Ind2b) ,CRLB_TNAA(Ind2b), CRLB_TNAA(Ind2b), std(b(Ind2b,:),0,2), std(b(Ind2b,:),0,2), 'bx');
hold on; P2= plot(bc1_m(Ind2b) , Y(Ind2b), 'b');
Y = Func(Out.TNAA.x3,bc1_m_v);
errorbar( bc1_m(Ind3b), TNAA(Ind3b) ,CRLB_TNAA(Ind3b), CRLB_TNAA(Ind3b), std(b(Ind3b,:),0,2), std(b(Ind3b,:),0,2), 'gx');
hold on; P3= plot(bc1_m(Ind3b) , Y(Ind3b), 'g');
title('NAA diffusivity')
set(gca,'YScale','log')
legend([P1,P2,P3], {'Dir1','Dir2','Dir3'});
xlabel('b-value');ylabel('Signal (A.u.)');
%
C = figure;
Func=@(A,xdata) exp(xdata*[A(1); A(2)]);
Y = Func(Out.TCho.x1,bc1_m_v);
errorbar( bc1_m(Ind1b), TCho(Ind1b) ,CRLB_TCho(Ind1b), CRLB_TCho(Ind1b), std(b(Ind1b,:),0,2), std(b(Ind1b,:),0,2), 'rx');
hold on; P1= plot(bc1_m(Ind1b) , Y(Ind1b), 'r');
Y = Func(Out.TCho.x2,bc1_m_v);
errorbar( bc1_m(Ind2b), TCho(Ind2b) ,CRLB_TCho(Ind2b), CRLB_TCho(Ind2b), std(b(Ind2b,:),0,2), std(b(Ind2b,:),0,2), 'bx');
hold on; P2= plot(bc1_m(Ind2b) , Y(Ind2b), 'b');
Y = Func(Out.TCho.x3,bc1_m_v);
errorbar( bc1_m(Ind3b), TCho(Ind3b) ,CRLB_TCho(Ind3b), CRLB_TCho(Ind3b), std(b(Ind3b,:),0,2), std(b(Ind3b,:),0,2), 'gx');
hold on; P3= plot(bc1_m(Ind3b) , Y(Ind3b), 'g');
title('Choline diffusivity')
set(gca,'YScale','log')
legend([P1,P2,P3], {'Dir1','Dir2','Dir3'});
xlabel('b-value');ylabel('Signal (A.u.)');

if exist('OutLoc','var')
    CJ_SavFig(A,'b_dist',OutLoc);
    CJ_SavFig(B,'TNAA',OutLoc);
    CJ_SavFig(C,'TCho',OutLoc);
end

end