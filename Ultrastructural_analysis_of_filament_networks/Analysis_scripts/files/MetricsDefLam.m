data(:,[2 3 4])=data(:,[2 3 4])*spatres;
FilsInds=unique(data(:,1));
numbFil=length(FilsInds);
DeletIndx=[]; %remoive filaments smaller than 3 points.
for i=1:numbFil
    iindi=FilsInds(i);
    indx=find(data(:,1)==iindi);
    if numel(indx)<3
      DeletIndx=[DeletIndx;indx];  
    end
end
data(DeletIndx,:)=[];
for i=1:numbFil
    iindi=FilsInds(i);
    indx=find(data(:,1)==iindi);
    if length(indx)<=3
        data(indx,:)=[];
    end
end
FilsInds=unique(data(:,1));
numbFil=length(FilsInds);
for i=1:numbFil
    iindi=FilsInds(i);
    indx=find(data(:,1)==iindi);
    Data(i).ind=i;
    Data(i).X=data(indx,2);
    Data(i).Y=data(indx,3);
    Data(i).Z=data(indx,4);
end
figure(1); hold on
filChange=find(diff(data(:,1)));
dataplot=data;
dataplot(filChange,:)=nan;
plot(dataplot(:,2),dataplot(:,3),'color','k','LineWidth',1.5);
axis equal;
LEdirisOkay=0;
prompt4={'Do you confirm the direction of lamellipodium?'};
while LEdirisOkay==0
    title({['file: \color{red}'  strrep(strrep(filename,'_','-'),'.txt','')],'\color{black}Select two points as start and end points of a vector' , 'which shows the direction of lamellipodium (direction of leading egde)'})
    dirofLE=ginput(2);
    LEdir=quiver(dirofLE(1,1),dirofLE(1,2),(dirofLE(2,1)-dirofLE(1,1)),(dirofLE(2,2)-dirofLE(1,2)),'LineWidth',5,'color','r');
    ConfrimDir=questdlg(prompt4,'','No','Yes','Yes');
    if strcmp(ConfrimDir,'Yes')==1
        LEdirisOkay=1;
    else
        delete(LEdir)
    end
end
close(figure(1));

% Mean axis of the filament
FilAxisEnd2EndLen=sqrt((dirofLE(2,1)-dirofLE(1,1)).^2+(dirofLE(2,2)-dirofLE(1,2)).^2);
projectOfFirstpointonAxis=min((data(:,2:3)- repmat(dirofLE(1,:),[size(data,1),1]) )*(dirofLE(2,:)-dirofLE(1,:))'/FilAxisEnd2EndLen);
tmp1=dirofLE(1,:)+projectOfFirstpointonAxis*(dirofLE(2,:)-dirofLE(1,:))/FilAxisEnd2EndLen;
projectOfLastpointonAxis=max((data(:,2:3)- repmat(dirofLE(1,:),[size(data,1),1]) )*(dirofLE(2,:)-dirofLE(1,:))'/FilAxisEnd2EndLen);
tmp2=dirofLE(1,:)+projectOfLastpointonAxis*(dirofLE(2,:)-dirofLE(1,:))/FilAxisEnd2EndLen;
dirofLE=[tmp1;tmp2];

AvgZ=mean(data(:,4));
FilAxisX=linspace(dirofLE(1,1),dirofLE(2,1),500);
FilAxisY=linspace(dirofLE(1,2),dirofLE(2,2),500);
FilAxisZ=AvgZ*ones(size(FilAxisX));
FilAxisdelX=gradient(FilAxisX);
FilAxisdelY=gradient(FilAxisY);
FilAxisdelZ=gradient(FilAxisZ);
FilAxisdelLentmp=sqrt(diff(FilAxisX).^2+diff(FilAxisY).^2);
FilAxisdelLen=[FilAxisdelLentmp FilAxisdelLentmp(end)];
FilAxislen=[0 cumsum(FilAxisdelLentmp)];
FilAxisLength=FilAxislen(end);
FilAxisEnd2EndLen=sqrt((FilAxisX(end)-FilAxisX(1)).^2+(FilAxisY(end)-FilAxisY(1)).^2+(FilAxisZ(end)-FilAxisZ(1)).^2);
FilsInds=unique(data(:,1));

numbFil=length(FilsInds);
for i=1:numbFil
    iindi=FilsInds(i);
    indx=find(data(:,1)==iindi);
    Data(i).ind=i;
    Data(i).X=data(indx,2);
    Data(i).Y=data(indx,3);
    Data(i).Z=data(indx,4);
    Data(i).delX=diff(data(indx,2));
    Data(i).delY=diff(data(indx,3));
    Data(i).delZ=diff(data(indx,4));
    Data(i).dell=sqrt(Data(i).delX.^2 +Data(i).delY.^2 +Data(i).delZ.^2);
    DelX(i)=Data(i).X(end)-Data(i).X(1); % end to end vector
    DelY(i)=Data(i).Y(end)-Data(i).Y(1);
    DelZ(i)=Data(i).Z(end)-Data(i).Z(1);
    DelL(i)=sqrt(DelX(i)^2+DelY(i)^2+DelZ(i)^2); % end to end distance
    Lengths(i)=sum(Data(i).dell);
    Bendiness(i)=Lengths(i)/DelL((i));
    [L,R,k] = curvature(data(indx,2:4));
    K=1./R; K(isnan(K))=0; K(1)=[]; K(end)=[];
    LengthsPoints=(Data(i).dell(1:end-1)+Data(i).dell(2:end))/2;
    BendingEnergy(i)=sum(LengthsPoints.*K.^2)./(Lengths(i).^2);
    LengthofAxisSegment=sqrt((FilAxisX(end)-FilAxisX(1)).^2+(FilAxisY(end)-FilAxisY(1)).^2+(FilAxisZ(end)-FilAxisZ(1)).^2);
    dxofAxisSegment=(FilAxisX(end)-FilAxisX(1))/LengthofAxisSegment;
    dyofAxisSegment=(FilAxisY(end)-FilAxisY(1))/LengthofAxisSegment;
    dzofAxisSegment=(FilAxisZ(end)-FilAxisZ(1))/LengthofAxisSegment;
    % and then we find the angle between end2end vector of filament and corresponding segment on the axis of filopodium
    Angle(i)=acosd((DelX(i).*dxofAxisSegment+DelY(i).*dyofAxisSegment+DelZ(i).*dzofAxisSegment)./(sqrt(DelX(i).^2+DelY(i).^2+DelZ(i).^2)));
    AngletoZ(i)=acosd(DelZ(i)./(sqrt(DelX(i).^2+DelY(i).^2+DelZ(i).^2)));
    AngletoXY(i)=acosd((DelX(i).*dxofAxisSegment+DelY(i).*dyofAxisSegment)./(sqrt(DelX(i).^2+DelY(i).^2))*(sqrt(dxofAxisSegment.^2+dyofAxisSegment.^2)));
    % dot product of filament axes and structure axis
    FilDir=sign([(dirofLE(2,:)-dirofLE(1,:)) 0]*[DelX(i);DelY(i);DelZ(i)]);
    if FilDir>=0 %the filament alignment is simuilar to that of structure
        PEndPosition(i,:)=[Data(i).X(1) Data(i).Y(1) Data(i).Z(1)]; BEndPosition(i,:)=[Data(i).X(end) Data(i).Y(end) Data(i).Z(end)];
    else
        BEndPosition(i,:)=[Data(i).X(1) Data(i).Y(1) Data(i).Z(1)]; PEndPosition(i,:)=[Data(i).X(end) Data(i).Y(end) Data(i).Z(end)];
    end
end
Angle=90-abs(Angle-90);
AngletoZ=90-abs(AngletoZ-90);
AngletoXY=90-abs(AngletoXY-90);
AngleDist=histcounts(AngletoXY,[0:90]);
AngleDist=AngleDist/sum(AngleDist);
Anisotropy=sum((AngleDist-1/90).^2);
%Barbd and pointed ends
NumOfCrossSectionsEnds=100;
WhichLen=linspace(0,1,NumOfCrossSectionsEnds+1);
for jj=1:length(WhichLen)
    whichLen=WhichLen(jj);
    whichind=round(interp1(FilAxislen,[1:length(FilAxislen)],whichLen*FilAxisLength));
    whichX=FilAxisX(whichind);whichY=FilAxisY(whichind);whichZ=FilAxisZ(whichind);
    whichdelX=FilAxisdelX(whichind)/FilAxisdelLen(whichind);whichdelY=FilAxisdelY(whichind)/FilAxisdelLen(whichind);whichdelZ=FilAxisdelZ(whichind)/FilAxisdelLen(whichind);
    FilAxisPoint=[whichX whichY whichZ]';
    FilAxisVector=[whichdelX whichdelY whichdelZ]';
    tmpBarbed(:,jj)=sign(FilAxisVector(1)*(BEndPosition(:,1)-FilAxisPoint(1))+FilAxisVector(2)*(BEndPosition(:,2)-FilAxisPoint(2))+FilAxisVector(3)*(BEndPosition(:,3)-FilAxisPoint(3)));
    tmpPointed(:,jj)=sign(FilAxisVector(1)*(PEndPosition(:,1)-FilAxisPoint(1))+FilAxisVector(2)*(PEndPosition(:,2)-FilAxisPoint(2))+FilAxisVector(3)*(PEndPosition(:,3)-FilAxisPoint(3)));
end
[I,J]=find(diff(tmpBarbed')');
BarbedEndNorm(I)=J/NumOfCrossSectionsEnds;
[I,J]=find(diff(tmpPointed')');
PointedEndNorm(I)=J/NumOfCrossSectionsEnds;
BarbedEnd=BarbedEndNorm*FilAxisEnd2EndLen; PointedEnd=PointedEndNorm*FilAxisEnd2EndLen;
%Filtering filaments based on length, angle, bendiness
ind=[];delind=[];FilsInds=unique(data(:,1));
for i=1:numbFil
    if Lengths(i)<str2num(filterings{2}) || Lengths(i)>str2num(filterings{3}) || Angle(i)<str2num(filterings{4}) ...
            || Angle(i)>str2num(filterings{5}) || Bendiness(i)<str2num(filterings{6}) || Bendiness(i)>str2num(filterings{7}) ...
            || AngletoZ(i)<str2num(filterings{8}) || AngletoZ(i)>str2num(filterings{9})
        ind=[ind ; i];
        delind=[delind;find(data(:,1)==FilsInds(i))];
    end
end
Data(ind)=[]; data(delind,:)=[]; DelX(ind)=[];  DelY(ind)=[]; DelZ(ind)=[]; DelL(ind)=[];
Lengths(ind)=[]; Bendiness(ind)=[]; BendingEnergy(ind)=[]; Angle(ind)=[]; AngletoZ(ind)=[]; AngletoXY(ind)=[];
PointedEnd(ind)=[];BarbedEnd(ind)=[];PointedEndNorm(ind)=[];BarbedEndNorm(ind)=[]; BEndPosition(ind,:)=[]; PEndPosition(ind,:)=[];
numbFil=size(unique(data(:,1)),1);
%% properties of filopodium at the tip, base and beam
WhichLen=linspace(0.02,0.98,NumOfCrossSections);
d = delaunayTriangulation(data(:,2), data(:,3), data(:,4));
for jj=1:length(WhichLen)
    whichLen=WhichLen(jj);
    whichind=round(interp1(FilAxislen, [1:length(FilAxislen)] ,whichLen*FilAxisLength));
    whichX=FilAxisX(whichind);whichY=FilAxisY(whichind);whichZ=FilAxisZ(whichind);
    whichdelX=FilAxisdelX(whichind)/FilAxisdelLen(whichind);whichdelY=FilAxisdelY(whichind)/FilAxisdelLen(whichind);whichdelZ=FilAxisdelZ(whichind)/FilAxisdelLen(whichind);
    FilAxisPoint=[whichX whichY whichZ]';
    FilAxisVector=[whichdelX whichdelY whichdelZ]';
    j=1;CrossedFilaments=[];Xcrs=[];Ycrs=[];Zcrs=[];
    for i=1:numbFil
        pointsinFil=[1:length(Data(i).X)];
        tmp=FilAxisVector(1)*(Data(i).X-FilAxisPoint(1))+FilAxisVector(2)*(Data(i).Y-FilAxisPoint(2))+FilAxisVector(3)*(Data(i).Z-FilAxisPoint(3));
        [~,DupInd]=unique(tmp);
        duplicateValues=setdiff([1:length(tmp)],DupInd);
        tmp(duplicateValues)=[];pointsinFil(duplicateValues)=[];
        if tmp(1)*tmp(end)<0 % if the filament number i crosses the plain of interest
            SolInd=interp1(tmp,pointsinFil',0);
            Xcrs(j)=interp1(pointsinFil,Data(i).X(pointsinFil),SolInd);
            Ycrs(j)=interp1(pointsinFil,Data(i).Y(pointsinFil),SolInd);
            Zcrs(j)=interp1(pointsinFil,Data(i).Z(pointsinFil),SolInd);
            CrossedFilaments(j)= i;
            j=j+1;
        end
    end
    Seglength(jj)=nanmean(Lengths(CrossedFilaments));
    SegBendiness(jj)=nanmean(Bendiness(CrossedFilaments));
    SegBendingEnergy (jj)=nanmean(BendingEnergy(CrossedFilaments));
    SegAngle(jj)=nanmean(Angle(CrossedFilaments));
    SegFilNum(jj)=length(CrossedFilaments);
    R1=FilAxisVector; % new x is the normal vector of cross section plain
    R3=[0;0;1];
    R2=cross(R3,R1);
    R=[R1 R2 R3]';
    if numel(Xcrs)>0
        CoorRt=R*([Xcrs;Ycrs;Zcrs]-repmat(FilAxisPoint,[1,size(Xcrs,2)]));
    else
        CoorRt=[];
    end
    if size(CoorRt,2)>2
        kk=convhull(CoorRt(2,:),CoorRt(3,:));
        shp=alphaShape(CoorRt(2,:)',CoorRt(3,:)',inf);
        SegArea(jj)=area(shp);
        SegDensity(jj)=SegFilNum(jj)/SegArea(jj);
        SegHeight(jj)=max(CoorRt(3,:))-min(CoorRt(3,:));
    else % if there are less than 3 filament in this cross section
        SegArea(jj)=nan;
        SegDensity(jj)=nan;
        SegPerimeter(jj)=nan;
        SegRoundness(jj)=nan;
        SegLength(jj)=nan;
        SegWidth(jj)=nan;
        SegElongation(jj)=nan;
        Ix(jj)=nan;
        Iz(jj)=nan;
    end
end
%% Relative configuration of filaments
dataDirx=[];dataDiry=[];dataDirz=[];   LastPoinsX=[]; LastPoinsY=[];
LastPoinsZ=[]; FirstPoinsX=[];  FirstPoinsY=[]; FirstPoinsZ=[];
for  i=1:length(Data)
    dataDirx=[dataDirx ; gradient(Data(i).X)];
    dataDiry=[dataDiry ; gradient(Data(i).Y)];
    dataDirz=[dataDirz ; gradient(Data(i).Z)];
    LastPoinsX=[LastPoinsX Data(i).X(end)];
    LastPoinsY=[LastPoinsY Data(i).Y(end)];
    LastPoinsZ=[LastPoinsZ Data(i).Z(end)];
    FirstPoinsX=[FirstPoinsX Data(i).X(1)];
    FirstPoinsY=[FirstPoinsY Data(i).Y(1)];
    FirstPoinsZ=[FirstPoinsZ Data(i).Z(1)];
end
dataDirl=sqrt(dataDirx.^2+dataDiry.^2+dataDirz.^2);
dataDirx=dataDirx./dataDirl;
dataDiry=dataDiry./dataDirl;
dataDirz=dataDirz./dataDirl;
numberInvestigations=200;
everyPoint=round(numel(data(:,2))/numberInvestigations);
AllPoints=[1:everyPoint:numel(data(:,2))];
PairwiseFilsAngles=[];PairwiseFilsDistances=[];
for jj=1:length(AllPoints)
    filAxisPoint=[data(AllPoints(jj),2) data(AllPoints(jj),3) data(AllPoints(jj),4)]';
    filAxisVector=[dataDirx(AllPoints(jj)) dataDiry(AllPoints(jj)) dataDirz(AllPoints(jj))]';
    SideofStartPoints=filAxisVector(1)*(FirstPoinsX-filAxisPoint(1))+filAxisVector(2)*(FirstPoinsY-filAxisPoint(2))+filAxisVector(3)*(FirstPoinsZ-filAxisPoint(3));
    SideofEndPoints=filAxisVector(1)*(LastPoinsX-filAxisPoint(1))+filAxisVector(2)*(LastPoinsY-filAxisPoint(2))+filAxisVector(3)*(LastPoinsZ-filAxisPoint(3));
    FilsCrossing=find(SideofStartPoints.*SideofEndPoints<0);
    j=1;Xcrs=[];Ycrs=[];Zcrs=[]; Delxtmp=[];Delytmp=[]; Delztmp=[];
    for i=1:numel(FilsCrossing)
        tmp=filAxisVector(1)*(Data(FilsCrossing(i)).X-filAxisPoint(1))+filAxisVector(2)*(Data(FilsCrossing(i)).Y-filAxisPoint(2))+filAxisVector(3)*(Data(FilsCrossing(i)).Z-filAxisPoint(3));
        if numel(find(tmp(1:end-1).*tmp(2:end)==0))>0
            tmp=tmp+0.0001;
        end
        solInd=find(tmp(1:end-1).*tmp(2:end)<0);
        SolIndFloor=solInd(1);
        decimal=abs(tmp(SolIndFloor))/(abs(tmp(SolIndFloor))+abs(tmp(SolIndFloor+1)));
        Xcrs(j)= Data(FilsCrossing(i)).X(SolIndFloor)+ decimal*(Data(FilsCrossing(i)).X(SolIndFloor+1)-Data(FilsCrossing(i)).X(SolIndFloor));
        Ycrs(j)= Data(FilsCrossing(i)).Y(SolIndFloor)+ decimal*(Data(FilsCrossing(i)).Y(SolIndFloor+1)-Data(FilsCrossing(i)).Y(SolIndFloor));
        Zcrs(j)= Data(FilsCrossing(i)).Z(SolIndFloor)+ decimal*(Data(FilsCrossing(i)).Z(SolIndFloor+1)-Data(FilsCrossing(i)).Z(SolIndFloor));
        delxtmp=Data(FilsCrossing(i)).X(SolIndFloor+1)-Data(FilsCrossing(i)).X(SolIndFloor);
        delytmp=Data(FilsCrossing(i)).Y(SolIndFloor+1)-Data(FilsCrossing(i)).Y(SolIndFloor);
        delztmp=Data(FilsCrossing(i)).Z(SolIndFloor+1)-Data(FilsCrossing(i)).Z(SolIndFloor);
        delltmp=sqrt(delxtmp.^2+delytmp.^2+delztmp.^2);
        Delxtmp(j)=delxtmp./delltmp;Delytmp(j)=delytmp./delltmp;Delztmp(j)=delztmp./delltmp;
        j=j+1;
    end
    R1=filAxisVector; % new x is the normal vector of cross section plain
    R3=[0;0;1];
    R2=cross(R3,R1);
    R=[R1 R2 R3]';
    if numel(Xcrs)>0
        CoorRt=R*([Xcrs;Ycrs;Zcrs]-repmat(filAxisPoint,[1,size(Xcrs,2)]));
        tmpDistances=squareform(pdist(CoorRt(2:3,:)'));
        tmpAngles=acosd(min(1,abs(Delxtmp'*Delxtmp + Delytmp'*Delytmp + Delztmp'*Delztmp)));
        uppertri=triu(ones(size(tmpDistances)))-eye(size(tmpDistances));
        tmp2Distances=tmpDistances(uppertri==1);
        tmp2Angles=tmpAngles(uppertri==1);
    else
        CoorRt=[];
        tmp2Angles=[];tmp2Distances=[];
    end
    PairwiseFilsAngles=[PairwiseFilsAngles; tmp2Angles];
    PairwiseFilsDistances=[PairwiseFilsDistances; tmp2Distances];
end
%%
% The average lamellipodium properties

lammetrics.Name=filename;
%General filopodium properties
FilopShape=alphaShape(data(:,2),data(:,3),data(:,4),100);
lammetrics.Length=FilAxisLength;
lammetrics.Vol=volume(FilopShape);
lammetrics.Density=pi*(FilamentDiameter/2)^2 *sum(Lengths)/volume(FilopShape); %volume fraction of filaments
lammetrics.FilamentNumber=numbFil; %total number of filaments

%filament properties
lammetrics.FilamentLength=Lengths;
lammetrics.FilamentAngle=Angle;
lammetrics.FilamentBendiness=Bendiness;
lammetrics.FilamentBendingEnergy=BendingEnergy;
lammetrics.FilamentAngletoZ=AngletoZ;
lammetrics.FilamentAngletoXY=AngletoXY;
lammetrics.FilamentAnisotropy=Anisotropy;
lammetrics.FilamentBarbedEnds=BarbedEnd;
lammetrics.FilamentPointedEnds=PointedEnd;
lammetrics.FilamentBarbedEndsNorm=BarbedEndNorm;
lammetrics.FilamentPointedEndsNorm=PointedEndNorm;
lammetrics.FilamentBEndPosition= BEndPosition;
lammetrics.FilamentPEndPosition= PEndPosition;
%Avg and SD of filaments properties
lammetrics.FilamentLengthSD=nanstd(Lengths); %SD of length of filaments
lammetrics.FilamentLengthAVG=nanmean(Lengths); %Avg of length of filaments
lammetrics.FilamentAngleSD=nanstd(Angle); %SD of angle of filaments
lammetrics.FilamentAngleAVG=nanmean(Angle); %Avg of angle of filaments
lammetrics.FilamentBendinessSD=nanstd(Bendiness); %SD of bendiness of filaments
lammetrics.FilamentBendinessAVG=nanmean(Bendiness); %Avg of bendiness of filaments
lammetrics.FilamentBendingEnergySD=nanstd(BendingEnergy); %SD of bending energy of filaments
lammetrics.FilamentBendingEnergyAVG=nanmean(BendingEnergy); %Avg of bending energy of filaments

%Cross sectional properties
lammetrics.WhichLen=WhichLen;
lammetrics.Seglength=Seglength;
lammetrics.SegAngle=SegAngle;
lammetrics.SegBendiness=SegBendiness;
lammetrics.SegBendingEnergy=SegBendingEnergy;
lammetrics.SegArea=SegArea;
lammetrics.SegFilNum=SegFilNum;
lammetrics.SegDensity=SegDensity;
lammetrics.SegVolFraction=pi*(FilamentDiameter/2)^2 *SegDensity;%cross sectional volume fraction
lammetrics.SegHeight=SegHeight;
lammetrics.PairwiseFilsAngles=PairwiseFilsAngles;
lammetrics.PairwiseFilsDistances=PairwiseFilsDistances;

%Avg and SD of cross sectional properties
lammetrics.FilamentSeglengthAVG=nanmean(Seglength); %average of crossection area
lammetrics.FilamentSeglengthSD=nanstd(Seglength); %SD of crossection area
lammetrics.FilamentSegAngleAVG=nanmean(SegAngle); %average of crossection area
lammetrics.FilamentSegAngleSD=nanstd(SegAngle); %SD of crossection area
lammetrics.FilamentSegBendinessAVG=nanmean(SegBendiness); %average of crossection area
lammetrics.FilamentSegBendinessSD=nanstd(SegBendiness); %SD of crossection area
lammetrics.FilamentSegBendingEnergyAVG=nanmean(SegBendingEnergy); %average of crossection area
lammetrics.FilamentSegBendingEnergySD=nanstd(SegBendingEnergy); %SD of crossection area
lammetrics.FilamentHeightAVG=nanmean(SegHeight); %average of crossection area
lammetrics.FilamentHeightSD=nanstd(SegHeight); %SD of crossection area
lammetrics.FilamentDensityAVG=nanmean(SegDensity); %average of crossection area
lammetrics.FilamentDensitySD=nanstd(SegDensity); %SD of crossection area
lammetrics.FilamentVolFrAVG=nanmean(pi*(FilamentDiameter/2)^2 *SegDensity); %average of crossection area
lammetrics.FilamentVolFrSD=nanstd(pi*(FilamentDiameter/2)^2 *SegDensity); %SD of crossection area
lammetrics.FilamentAreaAVG=nanmean(SegArea); %average of crossection area
lammetrics.FilamentAreaSD=nanstd(SegArea); %SD of crossection area
lammetrics.FilamentAreaNumAVG=nanmean(SegFilNum); %AVG of crossection number of filaments in each cross section
lammetrics.FilamentAreaNumSD=nanstd(SegFilNum);

% Base to tip ratio of cross sectional parameters
lammetrics.FilamentSeglengthBtT=nanmean(Seglength(1:round(length(Seglength)/2)))/nanmean(Seglength(round(length(Seglength)/2):end)); %average of crossection area
lammetrics.FilamentSegAngleBtT=nanmean(SegAngle(1:round(length(SegAngle)/2)))/nanmean(SegAngle(round(length(SegAngle)/2):end)); %average of crossection area
lammetrics.FilamentSegBendinessBtT=nanmean(SegBendiness(1:round(length(SegBendiness)/2)))/nanmean(SegBendiness(round(length(SegBendiness)/2):end)); %average of crossection area
lammetrics.FilamentSegBendingEnergyBtT=nanmean(SegBendingEnergy(1:round(length(SegBendingEnergy)/2)))/nanmean(SegBendingEnergy(round(length(SegBendingEnergy)/2):end)); %average of crossection area
lammetrics.FilamentHeightBtT=nanmean(SegHeight(1:round(length(SegHeight)/2)))/nanmean(SegHeight(round(length(SegHeight)/2):end)); %average of crossection area
lammetrics.FilamentDensityBtT= nanmean(SegDensity(1:round(length(SegDensity)/2)))/nanmean(SegDensity(round(length(SegDensity)/2):end)); %average of crossection area


