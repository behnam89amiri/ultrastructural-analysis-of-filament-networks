%--------------------------------------------------------------------------
%
%  Code for filtering filaments in the lamellipodia.
%
%  Written by Behnam Amiri
%
%--------------------------------------------------------------------------

clc
close all
clear all
warning('off','all');
currentfolder = pwd;
addpath('files')
prompt={'What do you want to analyse?'};
AnalysisType=questdlg(prompt,'Type of Analysis','Single lamellipodium','Multiple lamellipodia','Single lamellipodium'); % single or multiple lamellipodia
prompt2={'Pixel size (nm)','Min length of filaments (nm)','Max length of filaments (nm)','Min angle of filaments to the leading edge direction (deg)','Max angle of filaments to the leading edge direction (deg)','Min bendiness of filaments','Max bendiness of filaments','Min angle of filaments to the Zaxis (deg)','Max angle of filaments to the Zaxis (deg)'}; % filtering out filaments
definputs={'1','0','inf','0','inf','0','inf','30','inf'};
filterings=inputdlg(prompt2,'Include filaments in these ranges of properties',[1 100],definputs);
prompt3={'Min length of filaments (nm)','Max length of filaments (nm)','Min angle of filaments to the leading edge direction (deg)','Max angle of filaments to the leading edge direction (deg)','Min bendiness of filaments','Max bendiness of filaments','Min angle of filaments to the Zaxis (deg)','Max angle of filaments to the Zaxis (deg)'}; % filtering out filaments
definputs2={'0','0','0','0','0','0','0','0'};
filterings2=inputdlg(prompt3,'Exclude filaments in these ranges of properties',[1 100],definputs2);
spatres=str2num(filterings{1});
cd ..
parentFolder = pwd;
tmpind=exist('Data');
if tmpind~=7
    mkdir('Data');
end
cd('Data')
tmpind=exist('Lamellipodia_data');
if tmpind~=7
    mkdir('Lamellipodia_data');
end
cd('Lamellipodia_data')
defaultdir=pwd;
cd(currentfolder)
switch AnalysisType
    case 'Single lamellipodium'
        cd(defaultdir)
        [file,path]=uigetfile('*.txt','Select a text file to open');
        filename = file;
        cd(path)
        [data,dlm,~]= importdata(filename); %[Fil# x y z]
        cd(currentfolder)
        SubFltrLam;
        cd(path)
        [~,MainFileName,NameExtension]=fileparts(filename);
        if strcmp(dlm,' ')==1; dlm='\t'; end
        dlmwrite([MainFileName,'_filtered',NameExtension],dataFiltered,dlm)
        cd(currentfolder);
    case 'Multiple lamellipodia'
        folder = uigetdir(defaultdir,'Original data directory');
        ResultsPath=[folder '_Filtered'];
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
                CellCultures='';
            end
            wb = waitbar(0,['filtering the lamellipodia of ' CellCultures ' group...']);
            for cellnum = 1:length(files)
                try
                    clearvars -except spatres width height cellnum cellcult CellCultures files GroupFolder GroupFolderRes GroupFolderResAll  ResultsPath currentfolder  AnalysisType filterings filterings2 spatres folder subfolders kfilops kfilopsAll wb
                    filename = files(cellnum).name;
                    cd(GroupFolder)
                    [data,dlm,~]= importdata(filename); %[Fil# x y z]
                    cd(currentfolder);
                    close all
                    SubFltrLam;
                    cd(GroupFolderRes)
                    if strcmp(dlm,' ')==1; dlm='\t'; end
                    dlmwrite(filename,dataFiltered,dlm)
                    cd(currentfolder);
                    kfilops=kfilops+1;
                    kfilopsAll=kfilopsAll+1;
                    display([num2str(kfilopsAll) ' lamellipodia has been filtered']);
                catch
                end
                waitbar((cellnum)/length(files),wb);
            end
            pause(0.1)
            close(wb);close all;
        end
end