%--------------------------------------------------------------------------
%
%  Code for visualization of single structures together.
%
%  Written by Behnam Amiri
%
%--------------------------------------------------------------------------

clc
close all
clear all
warning('off','all');
currentfolder=pwd;
cd ..
tmpind=exist('Analysis_results');
if tmpind~=7
    mkdir('Analysis_results');
end
cd('Analysis_results')
defaultfolder=pwd;
cd(currentfolder)
folder = uigetdir(defaultfolder,'Select the results directory that contain data of desired sub-group');
cd(folder)
Cells=dir('*mat');
NumofCells=numel(Cells);
list={'length','relative orientation to axis','bendiness','relative orientation to Z-axis'};
[indx,td]=listdlg('ListString',list,'SelectionMode','single');
switch indx
    case 1
        defValue1=100;defValue2=1000;
    case 2
        defValue1=0;defValue2=90;
    case 3
        defValue1=1;defValue2=1.1;
    case 4
        defValue1=0;defValue2=90;
end
if floor(sqrt(NumofCells))*ceil(sqrt(NumofCells))>=NumofCells
    defrows=floor(sqrt(NumofCells));
    defcols=ceil(sqrt(NumofCells));
else
    defrows=ceil(sqrt(NumofCells));
    defcols=ceil(sqrt(NumofCells));
end
prompt2={'Number of subplot columns','Number of subplot rows','Thickness of filaments','Min value of the selected property for color coding','Max value of the selected property for color coding'};
definputs={num2str(defcols),num2str(defrows),'1',num2str(defValue1),num2str(defValue2)};
filterings=inputdlg(prompt2,'Inputs for the plots',[1 70],definputs);
plCol=str2num(filterings{1});
plRow=str2num(filterings{2});
filamentThickness=str2num(filterings{3});
minCustProp=str2num(filterings{4});
maxCustProp=str2num(filterings{5});
colresolution=50;
Cols=jet(colresolution);
minProp=minCustProp;maxProp=maxCustProp;
refProp=linspace(minProp,maxProp,colresolution);
ftmp=figure('visible','off');imagesc([minProp,maxProp]),ctmp=colorbar;ctmpTicks=ctmp.Ticks;ctmpTickLabels=ctmp.TickLabels;
close(ftmp);
pls=figure;
for i=1:NumofCells
    
    if i>plCol*plRow
        error('Number of files exceed number of subplots')
    end
    load(Cells(i).name)
    
    switch indx
        case 1
            ColCodePara=Lengths; ctitle='length (nm)';
        case 2
            ColCodePara=Angle; ctitle='relative orientation to axis (deg)';
        case 3
            ColCodePara=Bendiness; ctitle='bendiness';
        case 4
            ColCodePara=AngletoZ; ctitle='relative orientation to Z-axis (deg)';
    end
    subplot(plRow,plCol,i)
    c(i)=colorbar;
    hold on
    for j=1:length(Data)
        col=[interp1(refProp,Cols(:,1)',max(minProp,min(maxProp,ColCodePara(j))))    interp1(refProp,Cols(:,2)',max(minProp,min(maxProp,ColCodePara(j))))    interp1(refProp,Cols(:,3)',max(minProp,min(maxProp,ColCodePara(j))))];
        plot3(Data(j).X,Data(j).Y,Data(j).Z,'color',col,'LineWidth',filamentThickness)
    end
    view(2);
    axis tight
    axis equal
    title(strrep(strrep(Cells(i).name,'.mat',''),'_','-'))
    box on;
    ax=gca;ax.XTickLabel={};ax.YTickLabel={};ax.ZTickLabel={};
    c(i).Ticks=interp1([minProp,maxProp],[0,1],ctmpTicks);
    c(i).TickLabels=ctmpTickLabels;
    colormap(Cols)
    c(i).Label.String=ctitle;
end
cd(currentfolder)
