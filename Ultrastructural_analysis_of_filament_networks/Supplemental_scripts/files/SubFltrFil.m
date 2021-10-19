data(:,[2 3 4])=data(:,[2 3 4])*spatres;
FilsInds=unique(data(:,1));
numbFil=length(FilsInds);
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
prompt4={'Do you confirm the selected point?'};
while LEdirisOkay==0
    title({['file: \color{red}' strrep(strrep(filename,'_','-'),'.txt','')],'\color{black}Select a point which is closer to the base of filopodium'})
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

AvgZ=mean(data(:,4));
[Xmin,indXmin]=min(data(:,2));[Xmax,indXmax]=max(data(:,2));
[Ymin,indYmin]=min(data(:,3));[Ymax,indYmax]=max(data(:,3));
if range(data(:,2))>range(data(:,3)) % if the filament axis is closer to current X axis, place it on new x
    FilX=Xmax-Xmin; FilY=data(indXmax,3)-data(indXmin,3);
    FilEnd2EndLength=sqrt(FilX^2+FilY^2);
    cosFilEnd2EndTheta=FilX/FilEnd2EndLength;
    sinFilEnd2EndTheta=FilY/FilEnd2EndLength;
else
    FilX=data(indYmax,2)-data(indYmin,2); FilY=Ymax-Ymin;
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
FilAxisX=[Xmin:Xmax];
FilAxisY=polyval(p,FilAxisX);
FilAxisZ=AvgZ*ones(size(FilAxisX));
FilAxisdelX=diff(FilAxisX);
FilAxisdelY=diff(FilAxisY);
FilAxisdelZ=diff(FilAxisZ);
FilAxisdelLen=sqrt(FilAxisdelX.^2+FilAxisdelY.^2);
FilAxislen=cumsum(FilAxisdelLen);
FilAxisLength=sum(FilAxisdelLen);
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
    
    %
    Lengths(i)=sum(Data(i).dell);
    Bendiness(i)=Lengths(i)/DelL((i));
    %
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
end
Angle=90-abs(Angle-90);
AngletoZ=90-abs(AngletoZ-90);
%Filtering filaments based on length, angle, bendiness
ind=[];delind=[];FilsInds=unique(data(:,1));
for i=1:numbFil
    if Lengths(i)<str2num(filterings{2}) || Lengths(i)>str2num(filterings{3}) || Angle(i)<str2num(filterings{4}) ...
            || Angle(i)>str2num(filterings{5}) || Bendiness(i)<str2num(filterings{6}) || Bendiness(i)>str2num(filterings{7}) ...
            || AngletoZ(i)<str2num(filterings{8}) || AngletoZ(i)>str2num(filterings{9}) || ...
            (Lengths(i)>str2num(filterings2{1}) && Lengths(i)<str2num(filterings2{2})) || ...
            (Angle(i)>str2num(filterings2{3}) && Angle(i)<str2num(filterings2{4})) || ...
            (Bendiness(i)>str2num(filterings2{5}) && Bendiness(i)<str2num(filterings2{6})) || ...
            (AngletoZ(i)>str2num(filterings2{7}) && AngletoZ(i)<str2num(filterings2{8}))
        
        ind=[ind ; i];
        delind=[delind;find(data(:,1)==FilsInds(i))];
    end
end
Data(ind)=[]; data(delind,:)=[]; DelX(ind)=[];  DelY(ind)=[]; DelZ(ind)=[]; DelL(ind)=[];
Lengths(ind)=[]; Bendiness(ind)=[]; BendingEnergy(ind)=[]; Angle(ind)=[]; AngletoZ(ind)=[];
numbFil=size(unique(data(:,1)),1);

%%
dataFiltered=data;
figure(2);hold on
IDs=unique(dataFiltered(:,1));
for i=1:length(IDs)
    ThisFil=find(dataFiltered(:,1)==IDs(i));
    plot3(dataFiltered(ThisFil,2),dataFiltered(ThisFil,3),dataFiltered(ThisFil,4),'k','LineWidth',1.5);
end
axis equal
xlabel('X (nm)'); ylabel('Y (nm)');
title({'result after filtering','press any key to continue'})
pause
close all