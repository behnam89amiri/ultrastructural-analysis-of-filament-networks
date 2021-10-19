%--------------------------------------------------------------------------
%
%  Main driver code for filopodia analysis.
%
%  Written by Behnam Amiri
%
%--------------------------------------------------------------------------

clc
close all
clear all
warning('off','all');
addpath('files')
prompt={'What do you want to analyse?'};
AnalysisType=questdlg(prompt,'Type of Analysis','Single filopodium','Multiple filopodia','Single filopodium'); % single or multiple filopodia
prompt2={['Pixel size (nm)'],'Min length of filaments (nm)','Max length of filaments (nm)','Min angle of filaments to the tip direction (deg)','Max angle of filaments to the tip direction (deg)','Min bendiness of filaments','Max bendiness of filaments','Min angle of filaments to the Zaxis (deg)','Max angle of filaments to the Zaxis (deg)','Diameter of filaments (nm)','Number of cross-sections along the structure'}; % filtering out filaments
definputs={'1','0','inf','0','inf','0','inf','30','inf','7','50'};
filterings=inputdlg(prompt2,'Include filaments in these ranges of properties',[1 70],definputs);
spatres=str2num(filterings{1}); %micrometere per pixel
FilamentDiameter=str2num(filterings{10}); %diameter of filament (to calculate volume fraction of filaments in the structure)
NumOfCrossSections=str2num(filterings{11}); 
currentfolder = pwd;
cd ..
parentFolder = pwd;
tmpind=exist('Data');
if tmpind~=7
    mkdir('Data');
end
cd('Data')
tmpind=exist('Filopodia_data');
if tmpind~=7
    mkdir('Filopodia_data');
end
cd('Filopodia_data')
defaultdir=pwd;
cd(currentfolder);
switch AnalysisType
    case 'Single filopodium'
        cd(defaultdir)
        [file,path]=uigetfile('*.txt','Select a text file to open');
        cd(currentfolder);
        filename = file;
        wb = waitbar(0,'Analysing the filopodium...');
        cd(path)
        data= importdata(filename); %[Fil# x y z]
        cd(currentfolder)
        MetricsDefFil
        cd(parentFolder)
        if ~exist('Analysis_results','dir')
            mkdir('Analysis_results');
        end
        cd('Analysis_results')
        if ~exist('filopodia_results','dir')
            mkdir('filopodia_results');
        end
        cd('filopodia_results')
        if ~exist('Results_of_single_filopodium_Analysis','dir')
            mkdir('Results_of_single_filopodium_Analysis');
        end
        cd('Results_of_single_filopodium_Analysis');
        save(strrep(filename,'.txt',''),'Data','data','Lengths','Bendiness','BendingEnergy','Angle','AngletoZ','numbFil','FilAxisLength',...
            'AvgZ','FilAxisX','FilAxisY','FilAxisZ','FilAxisdelX','FilAxisdelY','FilAxisdelZ','FilAxisdelLen','FilAxislen','filopmetrics');
        cd(currentfolder);
        close(wb);
    case 'Multiple filopodia'
        cd(parentFolder)
        if ~exist('Analysis_results','dir')
            mkdir('Analysis_results');
        end
        cd('Analysis_results')
        if ~exist('filopodia_results','dir')
            mkdir('filopodia_results');
        end
        cd(parentFolder);
        if ~exist('Data','dir')
            mkdir('Data');
        end
        cd(currentfolder)
        folder = uigetdir(defaultdir,'Select data directory');
        [~,folderName]=fileparts(folder(1:end-(folder(end)==filesep)));
        ResultsPath=fullfile(fullfile(fullfile(parentFolder,'Analysis_results'),'filopodia_results'),folderName) ;
        dirfiles = dir(folder);
        dirFlags = [dirfiles.isdir] & ~strcmp({dirfiles.name},'.') & ~strcmp({dirfiles.name},'..');
        subfolders = dirfiles(dirFlags);
        kfilopsAll=0;
        for cellcult = 1:length(subfolders )+1
            kfilops=0;
            if cellcult<=length(subfolders )
                CellCultures = subfolders(cellcult).name;
                files = dir(fullfile(strcat(folder, '/', CellCultures), '*.txt'));
                mkdir(ResultsPath,CellCultures)
                GroupFolderRes=fullfile(ResultsPath,CellCultures);
                GroupFolder=fullfile(folder,CellCultures);
            else
                if ~exist(ResultsPath,'dir')
                    mkdir(ResultsPath)
                end
                files = dir(fullfile(strcat(folder), '*.txt'));
                GroupFolderRes=ResultsPath;
                GroupFolder=folder;
                CellCultures=folderName;
            end
            wb = waitbar(0,['Analysing the filopodia of ' CellCultures ' group...']);
            for cellnum = 1:length(files)
                try
                    clearvars -except cellnum cellcult CellCultures files GroupFolder GroupFolderRes GroupFolderResAll  ResultsPath currentfolder defaultdir AnalysisType filterings spatres FilamentDiameter NumOfCrossSections folder subfolders kfilops kfilopsAll folderName wb
                    filename = files(cellnum).name;
                    cd(GroupFolder)
                    data= importdata(filename); %[Fil# x y z]
                    cd(currentfolder);
                    MetricsDefFil
                    cd(GroupFolderRes)
                    save(strrep(filename,'.txt',''),'Data','data','Lengths','Bendiness','BendingEnergy','Angle','AngletoZ','numbFil','FilAxisLength',...
                        'AvgZ','FilAxisX','FilAxisY','FilAxisZ','FilAxisdelX','FilAxisdelY','FilAxisdelZ','FilAxisdelLen','FilAxislen','filopmetrics');
                    cd(currentfolder);
                    kfilops=kfilops+1;
                    kfilopsAll=kfilopsAll+1;
                    display([num2str(kfilopsAll) ' filopodia has been analyzed']);
                catch    
                end
                waitbar((cellnum)/length(files),wb);
            end
            pause(0.1)
            close(wb);close all;
        end
end
cd(currentfolder)