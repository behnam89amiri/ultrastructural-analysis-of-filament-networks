%create actin structures
clear all
clc
close all
currentFolde=pwd;
cd ..
cd('Data');cd('Lamellipodia_data')
resultsFolder='ArtificialStructureData';
if exist(resultsFolder)~=7
    mkdir(resultsFolder)
end
cd(resultsFolder);ResultsFolder=pwd;
cd(currentFolde) 
Inputs={'Total number of filaments in structure','Minimum Z-angle condition','Z-angle increment','Fraction of tilted filaments','SD of Z-angle distribution',...
    'Avg length of filaments','SD of length of filaments','Avg angle of filaments','SD of engle of filaments'};
DefaultInputs={'400','10','20','0.1','0','75','10','0','10'};
Ans=inputdlg(Inputs,'Inputs',[1 60],DefaultInputs);
numberofFilaments=str2double(Ans{1});
MinZangCond=str2double(Ans{2});
IncZangCond=str2double(Ans{3});
FractZangCond=str2double(Ans{4});
ZAngleSigma=str2double(Ans{5});
LengthMu=str2double(Ans{6});
LengthSigma=str2double(Ans{7});
AngleMu=str2double(Ans{8});
AngleSigma=str2double(Ans{9});

Conditions=[90:-IncZangCond:MinZangCond];
maxX=400;maxY=300;maxZ=50;
SegmentLength=5;
PercentageOfTiltedFils=round(FractZangCond*numberofFilaments);
Ans2=questdlg('Do you want to see the structures?');
Plotting=strcmp(Ans2,'Yes');
for Icondition=1:numel(Conditions)
    
    AngleDistribution=normrnd(AngleMu,AngleSigma,[numberofFilaments,1]);
    LengthDistribution=normrnd(LengthMu,LengthSigma,[numberofFilaments,1]);
    ZAngleDistribution1=[];
      for j=1:Icondition
    ZAngleMu=Conditions(j);
    ZAngleDistributiontmp=normrnd(ZAngleMu,ZAngleSigma,[PercentageOfTiltedFils,1]);
    ZAngleDistribution1=[ZAngleDistribution1;ZAngleDistributiontmp];
      end
    ZAngleDistribution1=sign(rand(numel(ZAngleDistribution1),1)-0.5*ones(numel(ZAngleDistribution1),1)).*ZAngleDistribution1;
    ZAngleDistribution2=normrnd(90,ZAngleSigma,[numberofFilaments-numel(ZAngleDistribution1),1]);
    ZAngleDistribution=[ZAngleDistribution1;ZAngleDistribution2];
    ZAngleDistribution=ZAngleDistribution(randperm(numberofFilaments));
    
    XnormalizedVect=sind(ZAngleDistribution).*cosd(AngleDistribution);
    YnormalizedVect=sind(ZAngleDistribution).*sind(AngleDistribution);
    ZnormalizedVect=cosd(ZAngleDistribution);
    
    XMidPoint=rand([numberofFilaments,1])*maxX;
    YMidPoint=rand([numberofFilaments,1])*maxY;
    ZMidPoint=rand([numberofFilaments,1])*maxZ;
    
    XStartPoint=XMidPoint-XnormalizedVect.*LengthDistribution/2;
    XEndPoint=XMidPoint+XnormalizedVect.*LengthDistribution/2;
    YStartPoint=YMidPoint-YnormalizedVect.*LengthDistribution/2;
    YEndPoint=YMidPoint+YnormalizedVect.*LengthDistribution/2;
    ZStartPoint=ZMidPoint-ZnormalizedVect.*LengthDistribution/2;
    ZEndPoint=ZMidPoint+ZnormalizedVect.*LengthDistribution/2;
    
    NumberOfPointsOnFilaments=round(LengthDistribution/SegmentLength);
    AllFilamentData=[];
    
    tmp=max(abs(ZAngleDistribution));
    numCol=100;Cols=jet(numCol);
    if Plotting==1;ff=figure;hold on; end
    for i=1:numberofFilaments
        x=linspace(XStartPoint(i),XEndPoint(i),NumberOfPointsOnFilaments(i));
        y=linspace(YStartPoint(i),YEndPoint(i),NumberOfPointsOnFilaments(i));
        z=linspace(ZStartPoint(i),ZEndPoint(i),NumberOfPointsOnFilaments(i));
        OutofBox=find(x<0 | x>maxX | y<0 | y>maxY | z<0 | z>maxZ);
        x(OutofBox)=[];y(OutofBox)=[];z(OutofBox)=[];
        thisFilamentData=[i*ones(size(x')) x' y' z'];
        AllFilamentData=[AllFilamentData;thisFilamentData];
        
        if Plotting==1
            j=1+round(numCol*(90-abs(ZAngleDistribution(i)))/(95));
            plot3(x,y,z,'linewidth',1.5,'color',Cols(j,:))
        end
        
    end
    if Plotting==1
        axis equal;ax=gca;ax.XTickLabel={};ax.YTickLabel={};ax.ZTickLabel={};
        view([-20 30])
        pause
    end
    
    cd(ResultsFolder)
    dlmwrite(['Zangle=' num2str(round(ZAngleMu)) '.txt'],AllFilamentData,' ')
    cd ..
end


