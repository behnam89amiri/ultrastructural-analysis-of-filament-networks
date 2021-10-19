%--------------------------------------------------------------------------
%
%  Code for grouping analyzed lamellipodia files.
%
%  Written by Behnam Amiri
%
%--------------------------------------------------------------------------

clc
clear all
close all
warning off;
addpath('files')
LabelsDefLam;
LLabels=cell(1,size(Labels,2)+1);LLabels(1)={'Groups'};LLabels(2:size(Labels,2)+1)=Labels;
currentfolder=pwd;
cd ..
tmpind=exist('Analysis_results');
if tmpind~=7
    mkdir('Analysis_results');
end
cd('Analysis_results')
tmpind=exist('lamellipodia_results');
if tmpind~=7
    mkdir('lamellipodia_results');
end
cd('lamellipodia_results')
defaultfolder=pwd;
cd(currentfolder)
folder = uigetdir(defaultfolder,'Select the results directory that contain data of desired sub-groups');
wb = waitbar(0,'please wait...');
ResultsPath=folder;
dirfiles = dir(folder);
dirFlags = [dirfiles.isdir] & ~strcmp({dirfiles.name},'.') & ~strcmp({dirfiles.name},'..');
subfolders = dirfiles(dirFlags);
LamMetrics=[];
for cellcult = 1:length(subfolders )+1
    kfilops=0;CellNames={};   clear ddata
    try
        if cellcult<=length(subfolders )
            CellCultures = subfolders(cellcult).name;
            if strcmp(CellCultures,'DataSummary')==1
                continue
            end
            files = dir(fullfile(strcat(folder, '/', CellCultures), '*.mat'));
            GroupFolderResAll=fullfile(ResultsPath,'DataSummary');
            GroupFolder=fullfile(folder,CellCultures);
        else
            files = dir(fullfile(strcat(folder), '*.mat'));
            GroupFolderResAll=fullfile(ResultsPath,'DataSummary');
            GroupFolder=folder;
            CellCultures='';
        end
        clear LamMetrics;
        if length(files)>0
            for cellnum = 1:length(files)
                try
                    clearvars -except Labels LLabels CellNames ddata cellnum cellcult CellCultures files GroupFolder GroupFolderRes GroupFolderResAll  ResultsPath currentfolder  AnalysisType filterings spatres folder subfolders kfilops kfilopsAll LamMetrics wb
                    filename = files(cellnum).name;
                    cd(GroupFolder)
                    load(filename);
                    kfilops=kfilops+1;
                    LamMetrics(kfilops)=lammetrics;
                    CellNames{kfilops}=strrep(filename,'.mat','');
                    ddatat(1)=[lammetrics.FilamentNumber]';
                    ddatat(2)=[lammetrics.Density]';
                    ddatat(3)=[lammetrics.FilamentAreaAVG]';
                    ddatat(4)=[lammetrics.FilamentDensityAVG]';
                    ddatat(5)=[lammetrics.FilamentVolFrAVG]';
                    ddatat(6)=[lammetrics.FilamentAreaNumAVG]';
                    ddatat(7)=[lammetrics.FilamentHeightAVG]';
                    ddatat(8)=[lammetrics.FilamentAnisotropy]';
                    ddatat(9)=[lammetrics.FilamentLengthAVG]';
                    ddatat(10)=[lammetrics.FilamentAngleAVG]';
                    ddatat(11)=[lammetrics.FilamentBendinessAVG]';
                    ddatat(12)=[lammetrics.FilamentBendingEnergyAVG]';
                    ddatat(13)=[lammetrics.FilamentDensityBtT]';
                    ddatat(14)=[lammetrics.FilamentHeightBtT]';
                    ddatat(15)=[lammetrics.FilamentSeglengthBtT]';
                    ddatat(16)=[lammetrics.FilamentSegAngleBtT]';
                    ddatat(17)=[lammetrics.FilamentSegBendinessBtT]';
                    ddatat(18)=[lammetrics.FilamentSegBendingEnergyBtT]';
                    ddata(kfilops,:)=ddatat;
                catch
                end
            end
            celldata=cell(size(ddata,1)+1,size(ddata,2)+1);
            celldata(1,:)=LLabels;
            celldata(2:end,1)=CellNames;
            celldata(2:end,2:end)=num2cell(ddata);
            T = cell2table(celldata);
            
            if ~exist(GroupFolderResAll,'dir')
                mkdir(GroupFolderResAll)
            end
            if exist('LamMetrics')
                cd(GroupFolderResAll)
                save(['LamMetrics' CellCultures],'LamMetrics')
                writetable(T,['LamMetrics' CellCultures '.xls'])
                cd(currentfolder);
            end
        end
    catch
    end
    waitbar((cellcult)/(length(subfolders)+1),wb);
end
close(wb)
cd(currentfolder);
