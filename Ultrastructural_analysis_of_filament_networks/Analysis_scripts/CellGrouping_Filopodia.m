%--------------------------------------------------------------------------
%
%  Code for grouping analyzed filopodia files.
%
%  Written by Behnam Amiri
%
%--------------------------------------------------------------------------

clc
clear all
close all
warning off;
addpath('files')
LabelsDefFil;
LLabels=cell(1,size(Labels,2)+1);LLabels(1)={'Groups'};LLabels(2:size(Labels,2)+1)=Labels;
currentfolder=pwd;
cd ..
tmpind=exist('Analysis_results');
if tmpind~=7
    mkdir('Analysis_results');
end
cd('Analysis_results')
tmpind=exist('filopodia_results');
if tmpind~=7
    mkdir('filopodia_results');
end
cd('filopodia_results')
defaultfolder=pwd;
cd(currentfolder)
folder = uigetdir(defaultfolder,'Select the results directory that contain data of desired sub-groups');
wb = waitbar(0,'please wait...');
ResultsPath=folder;
dirfiles = dir(folder);
dirFlags = [dirfiles.isdir] & ~strcmp({dirfiles.name},'.') & ~strcmp({dirfiles.name},'..');
subfolders = dirfiles(dirFlags);
FilopMetrics=[];
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
        clear FilopMetrics;
        if length(files)>0
            for cellnum = 1:length(files)
                try
                    clearvars -except Labels LLabels CellNames ddata cellnum cellcult CellCultures files GroupFolder GroupFolderRes GroupFolderResAll  ResultsPath currentfolder  AnalysisType filterings spatres folder subfolders kfilops kfilopsAll FilopMetrics wb
                    filename = files(cellnum).name;
                    cd(GroupFolder)
                    load(filename);
                    kfilops=kfilops+1;
                    FilopMetrics(kfilops)=filopmetrics;
                    CellNames{kfilops}=strrep(filename,'.mat','');
                    ddatat(1)=[filopmetrics.FilamentNumber]';
                    ddatat(2)=[filopmetrics.Density]';
                    ddatat(3)=[filopmetrics.Length]';
                    ddatat(4)=[filopmetrics.Bendiness]';
                    ddatat(5)=[filopmetrics.BendingEnergy]';
                    ddatat(6)=[filopmetrics.FilamentAreaAVG]';
                    ddatat(7)=[filopmetrics.FilamentDensityAVG]';
                    ddatat(8)=[filopmetrics.FilamentVolFrAVG]';
                    ddatat(9)=[filopmetrics.FilamentAreaNumAVG]';
                    ddatat(10)=[filopmetrics.SegRoundnessAVG]';
                    ddatat(11)=[filopmetrics.SegIxAVG]';
                    ddatat(12)=[filopmetrics.SegIzAVG]';
                    ddatat(13)=[filopmetrics.FilamentAnisotropy]';
                    ddatat(14)=[filopmetrics.FilamentLengthAVG]';
                    ddatat(15)=[filopmetrics.FilamentAngleAVG]';
                    ddatat(16)=[filopmetrics.FilamentBendinessAVG]';
                    ddatat(17)=[filopmetrics.FilamentBendingEnergyAVG]';
                    ddatat(18)=[filopmetrics.FilamentAreaBtT]';
                    ddatat(19)=[filopmetrics.FilamentDensityBtT]';
                    ddatat(20)=[filopmetrics.FilamentAreaNumBtT]';
                    ddatat(21)=[filopmetrics.SegIxBtT]';
                    ddatat(22)=[filopmetrics.SegIzBtT]';
                    ddatat(23)=[filopmetrics.FilamentSeglengthBtT]';
                    ddatat(24)=[filopmetrics.FilamentSegAngleBtT]';
                    ddatat(25)=[filopmetrics.FilamentSegBendinessBtT]';
                    ddatat(26)=[filopmetrics.FilamentSegBendingEnergyBtT]';
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
            if exist('FilopMetrics')
                cd(GroupFolderResAll)
                save(['FilopMetrics' CellCultures],'FilopMetrics')
                writetable(T,['FilopMetrics' CellCultures '.xls'])
                cd(currentfolder);
            end
        end
    catch
    end
    waitbar((cellcult)/(length(subfolders)+1),wb);
end
close(wb)
cd(currentfolder)