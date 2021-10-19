%--------------------------------------------------------------------------
%
%  Code for cropping the region of interest.
%
%  Written by Behnam Amiri
%
%--------------------------------------------------------------------------

clc
close all
clear all
warning('off','all');
addpath('files')
currentfolder = pwd;
prompt={'What do you want to analyse?'};
AnalysisType=questdlg(prompt,'Type of Analysis','Single filopodium or lamellipodium','Multiple filopodia or lamellipodia','Single filopodium or lamellipodium'); % single or multiple lamellipodia
prompt2={['Pixel size (nm)'],'Width of the region of interest (nm)','Height of the region of interest (nm)'}; % filtering out filaments
definputs={'1','',''};
filterings=inputdlg(prompt2,'Inputs',[1 50],definputs);
spatres=str2num(filterings{1});
width=str2num(filterings{2});
height=str2num(filterings{3});
if numel(width)==0 || numel(height)==0
    error('Please enter values for height and width')
end
cd ..
parentFolder = pwd;
tmpind=exist('Data');
if tmpind~=7
    mkdir('Data');
end
cd('Data')
defaultdir=pwd;
cd(currentfolder)
switch AnalysisType
    case 'Single filopodium or lamellipodium'
        cd(defaultdir)
        [file,path]=uigetfile('*.txt','Select a text file to open');
        filename = file;
        cd(path)
        [data,dlm,~]= importdata(filename); %[Fil# x y z]
        cd(currentfolder)
        SubCr;
        cd(path)
        [~,MainFileName,NameExtension]=fileparts(filename);
        if strcmp(dlm,' ')==1; dlm='\t'; end
        dlmwrite([MainFileName,'_Cropped',NameExtension],dataCropped,dlm)
        cd(currentfolder);
    case 'Multiple filopodia or lamellipodia'
        folder = uigetdir(defaultdir,'Original data directory');
        ResultsPath=[folder '_Cropped'];
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
            wb = waitbar(0,['Analysing the lamellipodia of ' CellCultures ' group...']);
            for cellnum = 1:length(files)
                try
                    clearvars -except spatres width height cellnum cellcult CellCultures files GroupFolder GroupFolderRes GroupFolderResAll  ResultsPath currentfolder  AnalysisType filterings spatres folder subfolders kfilops kfilopsAll wb
                    filename = files(cellnum).name;
                    cd(GroupFolder)
                    [data,dlm,~]= importdata(filename); %[Fil# x y z]
                    cd(currentfolder);
                    close all
                    SubCr;
                    cd(GroupFolderRes)
                    if strcmp(dlm,' ')==1; dlm='\t'; end
                    dlmwrite(filename,dataCropped,dlm)
                    cd(currentfolder);
                    kfilops=kfilops+1;
                    kfilopsAll=kfilopsAll+1;
                    display([num2str(kfilopsAll) ' lamellipodia has been analyzed']);
                catch    
                end
                waitbar((cellnum)/length(files),wb);
            end
            pause(0.1)
            close(wb);close all;
        end
end


