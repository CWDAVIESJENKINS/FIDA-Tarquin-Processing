function [] = CJ_SavFig( FigHandle , Title , SavLoc , Export , Overwrite )
%Save a figure as .fig and .png at specified location
% Input:    FigHandle - Figure number of target figure
%           Title - name to save figure as 
%           SavLoc - Directory to save the figure in
% Opt:      Export = export format (PAPER)
%           Overwrite = 1 to disable overwrite protection
% Output:   none

%Handle optional variable
if ~exist('SavLoc','var')
    SavLoc = pwd;
end
if ~exist('Overwrite','var')
    Overwrite = 1;
end

%We wish to return to current directory
Home = pwd;
%change to save location
cd(SavLoc)

% unless override is enabled, check tarrget directory for file of same name
% and ask whether to overwrite
if Overwrite == 0 && (~isempty(dir(sprintf('%s.fig',Title))) || ~isempty(dir(sprintf('%s.png',Title)))>0)
    Proceed = GUI_YesNo('UserData','File of that name already exists! Overwrite it?');
else
    Proceed = true;
end

if Proceed
    %create string with extension for filename
    PNG = sprintf('%s.png' , Title);
    FIG = sprintf('%s.fig' , Title);

    if exist('Export' , 'var')
        if ~strcmp(Export,'')
            sdf(FigHandle , Export)
        end
    end

    %save as fig
    savefig(FigHandle , FIG);

    %Old matlab versions seem to struggle with print
    if verLessThan('Matlab','R2015a')
        saveas(FigHandle , PNG);
    else
        %save as png
        print(FigHandle , PNG , '-dpng');
    end
end

%return to home directory
cd(Home)
end
