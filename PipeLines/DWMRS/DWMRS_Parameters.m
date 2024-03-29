function[BVAL, Locations] = DWMRS_Parameters(Spec, Loc)

% if ~exist('Loc','var')
%     Loc = pwd;
%     fname = 'data';
% else 
%     [Loc,fname] = fileparts(Loc);
%     if isempty(fname)
%         fname = 'data';
%     end
% end
fname = 'data';
Locations.BVEC = [Loc,filesep,fname,'.bvec'];
Locations.BVAL = [Loc,filesep,fname,'.bval'];
Locations.EDDY = [Loc,filesep,fname,'.eddy_parameters'];

%% BVec
for J=1:length(Spec.Diff_Dir)
    BVEC(:,J) = LookupDIR(Spec.Diff_Dir(J));
end
 dlmwrite(Locations.BVEC,BVEC,' ');

%% Bval
b_in = [Spec.Diff_G
ones(1,length(Spec.Diff_G)).*Spec.p.RiseTime
ones(1,length(Spec.Diff_G)).*Spec.te]';

for J=1:length(Spec.Diff_Dir)
    BVAL(J) = LookupBVAL(b_in(J,:));
end
dlmwrite(Locations.BVAL,BVAL,' ');
%% Eddy params

EDDY = zeros(length(Spec.Diff_G),16);
dlmwrite(Locations.EDDY,EDDY,' ');

end
%%
function[Out] = LookupDIR(DIR)
    Vec = [
 0 0 0%0        
 0 0 0%1   
 0 0 0%2
 0 0 0%3
 0 0 0%4
 0 0 0%5
 0 0 0%6
 0 0 1%7
-1 0 0%8
 0 1 0%9
 0 0 0%10
 0 0 0%11
 0 0 0%12
];

Out = Vec(DIR+1,:);
end

%%
function[Out] = LookupBVAL(IN)
%Lookup table
%    G    Rise   te       b-val(pulled from Elena)
LookupTable = [
    59    12    74        1257
    89    12    74        2838
   118    12    74        4989
   177    12    74        11226
   236    12    74        19957
   295    12    74        31183
   50    12    70        557.76
   100    12    70        2260.6
   150    12    70        5107.4
   200    12    70        9099.3
   250    12    70        14237
   %295    12    70        19837
    75    12    70        1266.2
    59    12    70        871
   118    12    70        3482
   177    12    70        7835
   236    12    70        13929
   295    12    70        21763
    16     5    106       869
    32     5    106       3477
    48     5    106       7823
    64     5    106       13907
    80     5    106       21730
    16     6    116       1245
    24     6    116       2801
    32     6    116       4980
    48     6    116       11204
    64     6    116       19918
    80     6    116       31122
];
if IN(1)==0
    Out = 0;
    return
end

Search = find(sum(ismember(LookupTable(:,1:3), IN),2)==3);

if isempty(Search)
    error('Need more Elena values')
else
    Out = LookupTable(Search,4);
end
end