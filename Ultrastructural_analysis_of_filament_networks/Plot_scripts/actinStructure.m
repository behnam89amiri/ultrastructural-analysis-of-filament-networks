%--------------------------------------------------------------------------
%
%  Code for visualization of the properties of analyzed structure groups.
%
%  Written by Behnam Amiri
%
%--------------------------------------------------------------------------

clc
close all
clear all
warning off
addpath('files')
currentfolder=pwd;
cd ..
tmpind=exist('Plot_results');
if tmpind~=7
    mkdir('Plot_results');
end
tmpind=exist('Analysis_results');
if tmpind~=7
    mkdir('Analysis_results');
end
cd('Analysis_results')
defaultfolder=pwd;
cd(currentfolder)
folder = uigetdir(defaultfolder,'Select the results directory that contain data of desired sub-groups');
cd(folder)
wb = waitbar(0,'please wait...');
ResultsPath=folder;
dirfiles=dir(fullfile(folder,'**\*.*'));
dirFlags = [dirfiles.isdir] & ~strcmp({dirfiles.name},'.') & ~strcmp({dirfiles.name},'..');
subfolders = dirfiles(dirFlags);
jj=1;        ii=1;   in=1;GroupNumber=[];
for cellcult = 1:length(subfolders )+1
    try
        CellGroup = subfolders(cellcult).name;
        if strcmp(CellGroup,'DataSummary')==1
            continue
        end
        cd(subfolders(cellcult).folder)
        cd(CellGroup)
        files=dir('*.mat');
        if numel(files)==0
            continue
        end 
        for cellnum = 1:length(files)
            try
                clear filopmetrics lammetrics
                filename = files(cellnum).name;
                load(filename);
                if exist('filopmetrics','var')==1
                    ddatat(1,jj)=filopmetrics.FilamentLengthAVG;
                    ddatat(2,jj)=filopmetrics.FilamentBendinessAVG;
                    ddatat(3,jj)=filopmetrics.Density;
                    ddatat(4,jj)=filopmetrics.FilamentAnisotropy;
                elseif exist('lammetrics','var')==1
                    ddatat(1,jj)=lammetrics.FilamentLengthAVG;
                    ddatat(2,jj)=lammetrics.FilamentBendinessAVG;
                    ddatat(3,jj)=lammetrics.Density;
                    ddatat(4,jj)=lammetrics.FilamentAnisotropy;
                end
                GroupNumber(jj) = ii;
                jj=jj+1;
            catch
            end     
        end
        GroupName{ii}=strrep(CellGroup,'_','-');
        ii=ii+1;
    catch
    end
end
variableNames={'length','bendiness','density','anisotropy'};
color=lines(numel(variableNames));
AllGroupNumber=unique(GroupNumber);
Groups=categorical(GroupNumber,AllGroupNumber,GroupName);
Z=zscore(ddatat');
if size(Z,1)<3
    error('There is not enough files for comparison')
end
[coefs,scores]=pca(Z);
figure;
bp=biplot(coefs(:,1:2),'Scores',scores(:,1:2),'MarkerSize',20,'VarLabels',variableNames);
tags=get(bp,'tag');
bpt=bp(strcmp(tags,'obsmarker'));
grp=findgroups(GroupNumber);
for i=1:max(grp)
    set(bpt(grp==i),'Color',color(i,:),'DisplayName',GroupName{i})
end
bpvar1=bp(strcmp(tags,'varline'));
bpvar2=bp(strcmp(tags,'varmarker'));
set(bpvar1,'Color','k')
set(bpvar2,'Color','k')
Inds=[];BP=[];tmp1=find((strcmp(tags,'obsmarker')));
for i=1:numel(GroupName)
    tmp2=tmp1(grp==i);
    tmp3=tmp2(1);
    Inds=[Inds tmp3];
    BP=[BP, bp(tmp3)];
end
legend(BP,GroupName)
variableNames={'length (nm)','bendiness','density (nm^{-2})','anisotropy'};
figure
gpl=gplotmatrix(ddatat',[],Groups,color,[],[],[],'grpbars',variableNames);
cd(currentfolder)
close(wb)