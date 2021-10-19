clc
clear all
close all
currentfolder=pwd;
cd ..
cd('Analysis_results')
cd('lamellipodia_results')
defaultfolder=pwd;
cd(currentfolder)
folder = uigetdir(defaultfolder,'Select the results directory that contain data of Z-angle investigation');
cd(folder);
if exist('DataSummary')~=7
    error('Please run the script CellGrouping_Lamellipodia')
else
cd('DataSummary');load('LamMetrics.mat')
end
cd(currentfolder)
for ii=1:numel(LamMetrics)
    Zangle(ii)=nanmean(LamMetrics(ii).FilamentAngletoZ);
    CrossSectionalDensity(ii)=LamMetrics(ii).FilamentAreaNumAVG;
end
figure
plot(Zangle,CrossSectionalDensity,'d','MarkerFaceColor','b')
ax=gca;ax.FontSize=14; ylabel( 'Cross-sectional number of filaments') ;
xlabel('Mean angle of filaments to Z-axis')
box on

