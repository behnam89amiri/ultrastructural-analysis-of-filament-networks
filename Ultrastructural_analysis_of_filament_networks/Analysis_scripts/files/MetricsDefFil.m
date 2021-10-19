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
prompt4={'Do you confirm the selected point?'};
while LEdirisOkay==0
    title({['file: \color{red}'  strrep(strrep(filename,'_','-'),'.txt','')],'\color{black}Select a point which is closer to the base of filopodium'})
    CorofBase=ginput(1);
    LEdir=plot(CorofBase(1,1),CorofBase(1,2),'r*','LineWidth',5);
    ConfrimDir=questdlg(prompt4,'','No','Yes','Yes');
    if strcmp(ConfrimDir,'Yes')==1
        LEdirisOkay=1;
    else
        delete(LEdir)
    end
end
close(figure(1));
AvgZ=nanmean(data(:,4));
tmpp=find(diff(data(:,1))~=0);
EndPoints1=[tmpp; numel(data(:,1))];
EndPoints2=[1; tmpp+1 ];
EndPoints=union(EndPoints1,EndPoints2);
DelatEnds=data(EndPoints1,2:3)- data(EndPoints2,2:3);
MeanAngleofFilamentstoX=mean(abs(atand(DelatEnds(:,2)./DelatEnds(:,1))));
if MeanAngleofFilamentstoX<45 % if the filament axis is closer to current X axis, place it on new x
    [Xssort,indXssort]=sort(data(EndPoints,2));
    percOfEndsAtEdge=20;
    tmpPrc1=find(Xssort<prctile(Xssort,percOfEndsAtEdge));
    tmpPrc2=find(Xssort>prctile(Xssort,100-percOfEndsAtEdge));
    Xmin=mean(data(EndPoints(indXssort(tmpPrc1)),2));
    Xmax=mean(data(EndPoints(indXssort(tmpPrc2)),2));
    FilX=Xmax-Xmin; FilY=mean(data(EndPoints(indXssort(tmpPrc2)),3))-mean(data(EndPoints(indXssort(tmpPrc1)),3));
    FilEnd2EndLength=sqrt(FilX^2+FilY^2);
    cosFilEnd2EndTheta=FilX/FilEnd2EndLength;
    sinFilEnd2EndTheta=FilY/FilEnd2EndLength;
else
    [Yssort,indYssort]=sort(data(EndPoints,3));
    percOfEndsAtEdge=10;
    tmpPrc1=find(Yssort<prctile(Yssort,percOfEndsAtEdge));
    tmpPrc2=find(Yssort>prctile(Yssort,100-percOfEndsAtEdge));
    Ymin=mean(data(EndPoints(indYssort(tmpPrc1)),3));
    Ymax=mean(data(EndPoints(indYssort(tmpPrc2)),3));
    FilX=mean(data(EndPoints(indYssort(tmpPrc1)),2))-mean(data(EndPoints(indYssort(tmpPrc2)),2)); FilY=Ymax-Ymin;
    FilEnd2EndLength=sqrt(FilX^2+FilY^2);
    cosFilEnd2EndTheta=FilX/FilEnd2EndLength;
    sinFilEnd2EndTheta=FilY/FilEnd2EndLength;
end
X=data(:,2)*cosFilEnd2EndTheta+data(:,3)*sinFilEnd2EndTheta;
Y=data(:,3)*cosFilEnd2EndTheta-data(:,2)*sinFilEnd2EndTheta;
Z=data(:,4);
CoorofBase(1)=CorofBase(1)*cosFilEnd2EndTheta+CorofBase(2)*sinFilEnd2EndTheta;
CoorofBase(2)=CorofBase(2)*cosFilEnd2EndTheta+CorofBase(1)*sinFilEnd2EndTheta;
if  abs(CoorofBase(1)-min(X))>abs(CoorofBase(1)-max(X))  % if indicated base point is closer to the current end
    X=-X;Y=-Y;
end
data(:,2)=X;data(:,3)=Y;data(:,4)=Z;
[Xmin,indXmin]=min(data(:,2));[Xmax,indXmax]=max(data(:,2));
[Ymin,indYmin]=min(data(:,3));[Ymax,indYmax]=max(data(:,3));

% Mean axis of the filament
[Xsort,indX]=sort(data(:,2));
Ysort=data(indX,3);
p=polyfit(Xsort,Ysort,2);
FilAxisX=linspace(Xmin,Xmax,500);
FilAxisY=polyval(p,FilAxisX);
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
clear Data;
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
    % we find the segment on the filopodium which correspond to the filament (by finding closest point on axis)
    [~,endpointindx]=min((Data(i).X(end)-FilAxisX).^2 + (Data(i).Y(end)-FilAxisY).^2 + (Data(i).Z(end)-FilAxisZ).^2);
    [~,startpointindx]=min((Data(i).X(1)-FilAxisX).^2 + (Data(i).Y(1)-FilAxisY).^2 + (Data(i).Z(1)-FilAxisZ).^2);
    LengthofAxisSegment=sqrt((FilAxisX(endpointindx)-FilAxisX(startpointindx)).^2+(FilAxisY(endpointindx)-FilAxisY(startpointindx)).^2+(FilAxisZ(endpointindx)-FilAxisZ(startpointindx)).^2);
    dxofAxisSegment=(FilAxisX(endpointindx)-FilAxisX(startpointindx))/LengthofAxisSegment;
    dyofAxisSegment=(FilAxisY(endpointindx)-FilAxisY(startpointindx))/LengthofAxisSegment;
    dzofAxisSegment=(FilAxisZ(endpointindx)-FilAxisZ(startpointindx))/LengthofAxisSegment;
    % and then we find the angle between end2end vector of filament and corresponding segment on the axis of filopodium
    Angle(i)=acosd((DelX(i).*dxofAxisSegment+DelY(i).*dyofAxisSegment+DelZ(i).*dzofAxisSegment)./(sqrt(DelX(i).^2+DelY(i).^2+DelZ(i).^2)));
    AngletoZ(i)=acosd(DelZ(i)./(sqrt(DelX(i).^2+DelY(i).^2+DelZ(i).^2)));
    AngletoXY(i)=acosd((DelX(i).*dxofAxisSegment+DelY(i).*dyofAxisSegment)./(sqrt(DelX(i).^2+DelY(i).^2))*(sqrt(dxofAxisSegment.^2+dyofAxisSegment.^2)));
% dot product of filament axes and structure axis
    FilDir=sign([FilAxisX(end)-FilAxisX(1) FilAxisY(end)-FilAxisY(1) 0]*[DelX(i);DelY(i);DelZ(i)]);
    if FilDir>=0 %the filament alignment is similar to that of structure
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
    whichind=round(interp1(FilAxislen,[1:length(FilAxislen)],whichLen*FilAxisLength));
    whichX=FilAxisX(whichind);whichY=FilAxisY(whichind);whichZ=FilAxisZ(whichind);
    whichdelX=FilAxisdelX(whichind)/FilAxisdelLen(whichind);whichdelY=FilAxisdelY(whichind)/FilAxisdelLen(whichind);whichdelZ=FilAxisdelZ(whichind)/FilAxisdelLen(whichind);
    FilAxisPoint=[whichX whichY whichZ]';
    FilAxisVector=[whichdelX whichdelY whichdelZ]';
    j=1;CrossedFilaments=[];Xcrs=[];Ycrs=[];Zcrs=[];
    for i=1:numbFil
        pointsinFil=length(Data(i).X);
        tmp=FilAxisVector(1)*(Data(i).X-FilAxisPoint(1))+FilAxisVector(2)*(Data(i).Y-FilAxisPoint(2))+FilAxisVector(3)*(Data(i).Z-FilAxisPoint(3));
        if tmp(1)*tmp(end)<0 % if the filament number i crosses the plain of interest
            SolInd=interp1(tmp,[1:pointsinFil]',0);
            Xcrs(j)=interp1([1:pointsinFil],Data(i).X,SolInd);
            Ycrs(j)=interp1([1:pointsinFil],Data(i).Y,SolInd);
            Zcrs(j)=interp1([1:pointsinFil],Data(i).Z,SolInd);
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
        SegPerimeter(jj)=perimeter(shp);
        SegRoundness(jj)=4*pi*SegArea(jj)/(SegPerimeter(jj)^2);
        corners=minBoundingBox([CoorRt(2,kk);CoorRt(3,kk)]);
        dist1=sqrt((corners(1,1)-corners(1,2)).^2 +(corners(2,1)-corners(2,2)).^2);
        dist2=sqrt((corners(1,3)-corners(1,2)).^2 +(corners(2,3)-corners(2,2)).^2);
        SegLength(jj)=max(dist1,dist2);
        SegWidth(jj)=min(dist1,dist2);
        SegElongation(jj)=SegWidth(jj)/SegLength(jj);
        Ix(jj)=sum((CoorRt(2,:)-mean(CoorRt(2,:))).^2);
        Iz(jj)=sum((CoorRt(3,:)-mean(CoorRt(3,:))).^2);
        %     Distt=[10:10:100]*spatres;
        %     for ijk=1:length(Distt)
        %         distt=Distt(ijk);
        %     Kf(ijk)=kfunction2(CoorRt(2:3,:)',distt);
        %     end
        %     Lf=sqrt(Kf/pi);
        
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
%%  Relative configuration of filaments
dataDirx=[];dataDiry=[];dataDirz=[];    LastPoinsX=[]; LastPoinsY=[];
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
% The average filopodium properties

filopmetrics.Name=filename;
%General filopodium properties
FilopShape=alphaShape(data(:,2),data(:,3),data(:,4),100);
filopmetrics.Length=FilAxisLength;
filopmetrics.Thickness=2*sqrt(nanmean(SegArea))/pi;
filopmetrics.Vol=volume(FilopShape);
filopmetrics.Density=pi*(FilamentDiameter/2)^2 * sum(Lengths)/volume(FilopShape); %volume fraction of filaments
filopmetrics.Bendiness=FilAxisLength/FilAxisEnd2EndLen;
[L,R,k] = curvature([FilAxisX' FilAxisY' FilAxisZ']);
K=1./R; K(isnan(K))=0; K(1)=[]; K(end)=[];
LengthsPoints=(FilAxisdelLen(1:end-1)+FilAxisdelLen(2:end))/2;
filopmetrics.BendingEnergy=sum(LengthsPoints(1:end-1)'.*K.^2)./(FilAxisLength.^2);%average of bending energy
filopmetrics.FilamentNumber=numbFil; %total number of filaments

%filament properties
filopmetrics.FilamentLength=Lengths;
filopmetrics.FilamentAngle=Angle;
filopmetrics.FilamentBendiness=Bendiness;
filopmetrics.FilamentBendingEnergy=BendingEnergy;
filopmetrics.FilamentAngletoZ=AngletoZ;
filopmetrics.FilamentAngletoXY=AngletoXY;
filopmetrics.FilamentAnisotropy=Anisotropy;
filopmetrics.FilamentBarbedEnds=BarbedEnd;
filopmetrics.FilamentPointedEnds=PointedEnd;
filopmetrics.FilamentBarbedEndsNorm=BarbedEndNorm;
filopmetrics.FilamentPointedEndsNorm=PointedEndNorm;
filopmetrics.FilamentBEndPosition= BEndPosition;
filopmetrics.FilamentPEndPosition= PEndPosition;
%Avg and SD of filaments properties
filopmetrics.FilamentLengthSD=nanstd(Lengths); %SD of length of filaments
filopmetrics.FilamentLengthAVG=nanmean(Lengths); %Avg of length of filaments
filopmetrics.FilamentAngleSD=nanstd(Angle); %SD of angle of filaments
filopmetrics.FilamentAngleAVG=nanmean(Angle); %Avg of angle of filaments
filopmetrics.FilamentBendinessSD=nanstd(Bendiness); %SD of bendiness of filaments
filopmetrics.FilamentBendinessAVG=nanmean(Bendiness); %Avg of bendiness of filaments
filopmetrics.FilamentBendingEnergySD=nanstd(BendingEnergy); %SD of bending energy of filaments
filopmetrics.FilamentBendingEnergyAVG=nanmean(BendingEnergy); %Avg of bending energy of filaments

%Cross sectional properties
filopmetrics.WhichLen=WhichLen;
filopmetrics.Seglength=Seglength;
filopmetrics.SegAngle=SegAngle;
filopmetrics.SegBendiness=SegBendiness;
filopmetrics.SegBendingEnergy=SegBendingEnergy;
filopmetrics.SegFilNum=SegFilNum;
filopmetrics.SegDensity=SegDensity;
filopmetrics.SegVolFraction=pi*(FilamentDiameter/2)^2 *SegDensity; %cross sectional volume fraction
filopmetrics.SegArea=SegArea;
filopmetrics.SegPerimeter=SegPerimeter;
filopmetrics.SegRoundness=SegRoundness;
filopmetrics.SegElongation=SegElongation;
filopmetrics.SegLength=SegLength;
filopmetrics.SegWidth=SegWidth;
filopmetrics.Ix=Ix;
filopmetrics.Iz=Iz;

filopmetrics.PairwiseFilsAngles=PairwiseFilsAngles;
filopmetrics.PairwiseFilsDistances=PairwiseFilsDistances;

%Avg and SD of cross sectional properties
filopmetrics.FilamentSeglengthAVG=nanmean(Seglength); %average of crossection area
filopmetrics.FilamentSeglengthSD=nanstd(Seglength); %SD of crossection area
filopmetrics.FilamentSegAngleAVG=nanmean(SegAngle); %average of crossection area
filopmetrics.FilamentSegAngleSD=nanstd(SegAngle); %SD of crossection area
filopmetrics.FilamentSegBendinessAVG=nanmean(SegBendiness); %average of crossection area
filopmetrics.FilamentSegBendinessSD=nanstd(SegBendiness); %SD of crossection area
filopmetrics.FilamentSegBendingEnergyAVG=nanmean(SegBendingEnergy); %average of crossection area
filopmetrics.FilamentSegBendingEnergySD=nanstd(SegBendingEnergy); %SD of crossection area
filopmetrics.FilamentAreaAVG=nanmean(SegArea); %average of crossection area
filopmetrics.FilamentAreaSD=nanstd(SegArea); %SD of crossection area
filopmetrics.FilamentDensityAVG=nanmean(SegDensity); %average of crossection area
filopmetrics.FilamentDensitySD=nanstd(SegDensity); %SD of crossection area
filopmetrics.FilamentVolFrAVG=nanmean(pi*(FilamentDiameter/2)^2 *SegDensity); %average of crossection area
filopmetrics.FilamentVolFrSD=nanstd(pi*(FilamentDiameter/2)^2 *SegDensity); %SD of crossection area
filopmetrics.FilamentAreaNumAVG=nanmean(SegFilNum); %AVG of crossection number of filaments in each cross section
filopmetrics.FilamentAreaNumSD=nanstd(SegFilNum); %SD of crossection number of filaments in each cross section
filopmetrics.SegPerimeterAVG=nanmean(SegPerimeter);
filopmetrics.SegPerimeterSD=nanstd(SegPerimeter);
filopmetrics.SegRoundnessAVG=nanmean(SegRoundness);
filopmetrics.SegRoundnessSD=nanstd(SegRoundness);
filopmetrics.SegElongationAVG=nanmean(SegElongation);
filopmetrics.SegElongationSD=nanstd(SegElongation);
filopmetrics.SegLengthAVG=nanmean(SegLength);
filopmetrics.SegLengthSD=nanstd(SegLength);
filopmetrics.SegWidthAVG=nanmean(SegWidth);
filopmetrics.SegWidthSD=nanstd(SegWidth);
% Base to tip ratio of cross sectional parameters
filopmetrics.FilamentSeglengthBtT=nanmean(Seglength(1:round(length(Seglength)/2)))/nanmean(Seglength(round(length(Seglength)/2):end)); %average of crossection area
filopmetrics.FilamentSegAngleBtT=nanmean(SegAngle(1:round(length(SegAngle)/2)))/nanmean(SegAngle(round(length(SegAngle)/2):end)); %average of crossection area
filopmetrics.FilamentSegBendinessBtT=nanmean(SegBendiness(1:round(length(SegBendiness)/2)))/nanmean(SegBendiness(round(length(SegBendiness)/2):end)); %average of crossection area
filopmetrics.FilamentSegBendingEnergyBtT=nanmean(SegBendingEnergy(1:round(length(SegBendingEnergy)/2)))/nanmean(SegBendingEnergy(round(length(SegBendingEnergy)/2):end)); %average of crossection area
filopmetrics.FilamentAreaBtT=nanmean(SegArea(1:round(length(SegArea)/2)))/nanmean(SegArea(round(length(SegArea)/2):end)); %average of crossection area
filopmetrics.FilamentDensityBtT= nanmean(SegDensity(1:round(length(SegDensity)/2)))/nanmean(SegDensity(round(length(SegDensity)/2):end)); %average of crossection area
filopmetrics.FilamentAreaNumBtT=nanmean(SegFilNum(1:round(length(SegFilNum)/2)))/nanmean(SegFilNum(round(length(SegFilNum)/2):end)); %AVG of crossection number of filaments in each cross section
filopmetrics.SegPerimeterBtT= nanmean(SegPerimeter(1:round(length(SegPerimeter)/2)))/nanmean(SegPerimeter(round(length(SegPerimeter)/2):end));
filopmetrics.SegRoundnessBtT= nanmean(SegRoundness(1:round(length(SegRoundness)/2)))/nanmean(SegRoundness(round(length(SegRoundness)/2):end));
filopmetrics.SegElongationBtT= nanmean(SegElongation(1:round(length(SegElongation)/2)))/nanmean(SegElongation(round(length(SegElongation)/2):end));
filopmetrics.SegLengthBtT= nanmean(SegLength(1:round(length(SegLength)/2)))/nanmean(SegLength(round(length(SegLength)/2):end));
filopmetrics.SegWidthBtT= nanmean(SegWidth(1:round(length(SegWidth)/2)))/nanmean(SegWidth(round(length(SegWidth)/2):end));
%moment of inertia
filopmetrics.SegIxAVG=nanmean(Ix);
filopmetrics.SegIxSD=nanstd(Ix);
filopmetrics.SegIxBtT= nanmean(Ix(1:round(length(Ix)/2)))/nanmean(Ix(round(length(Ix)/2):end));
filopmetrics.SegIzAVG=nanmean(Iz);
filopmetrics.SegIzSD=nanstd(Iz);
filopmetrics.SegIzBtT= nanmean(Iz(1:round(length(Iz)/2)))/nanmean(Iz(round(length(Iz)/2):end));% filopmetrics.FilamentHomogeneity=; %SD of crossection number of filaments in each cross section
filopmetrics.SegIx2IzAVG=nanmean(Ix)/nanmean(Iz);


% filopmetrics.I=;%average of crossection second moment of area
% filopmetrics.Ix=;%average of crossection second moment of area
% filopmetrics.Iz=;%average of crossection second moment of area

