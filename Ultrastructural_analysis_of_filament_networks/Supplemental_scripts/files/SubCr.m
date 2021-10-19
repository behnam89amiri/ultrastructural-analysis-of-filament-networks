method=1; %1 is cutting filaments    2: is deleting filaments
data(:,[2 3 4])=data(:,[2 3 4])*spatres;
figure(1);hold on
IDs=unique(data(:,1));

filChange=find(diff(data(:,1)));
dataplot=data;
dataplot(filChange,:)=nan;
plot(dataplot(:,2),dataplot(:,3),'color','k','LineWidth',1.5);

Xc=(min(data(:,2))+max(data(:,2)))/2;
Yc=(min(data(:,3))+max(data(:,3)))/2;
axis equal;
xlabel('X (nm)');ylabel('Y (nm)');zlabel('Z (nm)')
LEdirisOkay=0;
prompt4={'Do you confirm the selected region?'};
while LEdirisOkay==0
    title({['file: \color{red}' strrep(strrep(filename,'_','-'),'.txt','')],'\color{black}Select the rectangular region of interest'})
    pause(0.1);
    roi=drawrectangle('Position',[Xc-(width/2) Yc-(height/2) width height],'Rotatable',1);
    pause(0.1);
    wait(roi)
    % Pts=roi.Position;
    Pts=roi.Vertices;
    angle=-roi.RotationAngle;
    cosFilEnd2EndTheta=cosd(angle);
    sinFilEnd2EndTheta=sind(angle);
    X=data(:,2)*cosFilEnd2EndTheta+data(:,3)*sinFilEnd2EndTheta;
    Y=data(:,3)*cosFilEnd2EndTheta-data(:,2)*sinFilEnd2EndTheta;
    Z=data(:,4);
    Vertices(1:4,1)=Pts(:,1)*cosFilEnd2EndTheta+Pts(:,2)*sinFilEnd2EndTheta;
    Vertices(1:4,2)=Pts(:,2)*cosFilEnd2EndTheta-Pts(:,1)*sinFilEnd2EndTheta;
    
    ConfrimDir=questdlg(prompt4,'','No','Yes','Yes');
    if strcmp(ConfrimDir,'Yes')==1
        LEdirisOkay=1;
    else
        delete(roi)
    end
end
datanew=[data(:,1) X Y Z];
FilsInds=unique(data(:,1));
numbFil=numel(FilsInds);
j=numbFil+1;
OutsiderPoints=[];dataCropped=datanew;
for i=1:numbFil
    iindi=FilsInds(i);
    indx=find(data(:,1)==iindi);
    NewID=i*ones(size(datanew(indx,1)));
    X1=min(Vertices(:,1));Y1=min(Vertices(:,2)); X2=max(Vertices(:,1));Y2=max(Vertices(:,2));
    outsiderPoints= find(datanew(indx,2)>X2 | datanew(indx,2)<X1 | datanew(indx,3)<Y1 | datanew(indx,3)>Y2 );
    if method==1
        [insiderPoints,insiderPointsInd]= setdiff([1:length(indx)],outsiderPoints);
        cutPoints=[find(diff(insiderPoints)>1) length(insiderPoints)];
        if length(cutPoints)>1
            for ii=2:length(cutPoints)
                NewID(insiderPointsInd(cutPoints(ii-1)+1:cutPoints(ii)))=j;
                j=j+1;
            end
        end
        dataCropped(indx,1)=NewID;
        OutsiderPoints=[OutsiderPoints;indx(outsiderPoints)];
    elseif method==2
        if isempty(outsiderPoints) ~=1
            dataCropped(indx,:)=[] ;
        end
    end
end
dataCropped(OutsiderPoints,:)=[];
%%
figure(2);hold on
filChange=find(diff(dataCropped(:,1)));
dataplot=dataCropped;
dataplot(filChange,:)=nan;
plot(dataplot(:,2),dataplot(:,3),'color','k','LineWidth',1.5);
plot([Vertices(:,1);Vertices(1,1)],[Vertices(:,2);Vertices(1,2)],'r','Linewidth',2)
axis equal
xlabel('new X (nm)'); ylabel('new Y (nm)');
title({'result after cropping','press any key to continue'})
pause
close all