function varargout = GUIPlotlam(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUIPlotlam_OpeningFcn, ...
    'gui_OutputFcn',  @GUIPlotlam_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUIPlotlam is made visible.
function GUIPlotlam_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIPlotlam (see VARARGIN)

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
tmpind=exist('lamellipodia_results');
if tmpind~=7
    mkdir('lamellipodia_results');
end
cd('lamellipodia_results')
[file,path]=uigetfile('*.mat','Select a lamellipodium file from analysis results');
cd(path);
load(file);
cd(currentfolder)
Cols=parula(50);
handles.file=file;
filamentThickness=1;
relColCode=1;
minCustProp=[];
maxCustProp=[];
Prop=1;
handles.filamentThickness=filamentThickness;
handles.relColCode=relColCode;
handles.minCustProp=minCustProp;
handles.maxCustProp=maxCustProp;
handles.Lengths=Lengths;
handles.Angle=Angle;
handles.AngletoZ=AngletoZ;
handles.Bendiness=Bendiness;
% handles.data=data;
handles.Data=Data;
handles.Cols=Cols;
handles.Prop=Prop;
handles.colBar=0;
handles.c=[];
handles.ax=0;
handles.shape=0;
handles.FilAxisX=FilAxisX;
handles.FilAxisY=FilAxisY;
handles.FilAxisZ=FilAxisZ;
handles.FilAxislen=FilAxislen;
handles.FilAxisdelX=FilAxisdelX;
handles.FilAxisdelY=FilAxisdelY;
handles.FilAxisdelZ=FilAxisdelZ;
handles.FilAxisdelLen=FilAxisdelLen;
handles.FilAxisLength=FilAxisLength;
handles.Length=lammetrics.Length;
handles.Height=lammetrics.FilamentHeightAVG;
handles.Number=lammetrics.FilamentNumber;
handles.data=data;

CrsSec2=0;
handles.CrsSec2=CrsSec2;
CrsSec3=0;
handles.CrsSec3=CrsSec3;
CrsSecDist=0.5;
handles.CrsSecDist=CrsSecDist;



if relColCode==1
    minProp=min(Lengths);maxProp=max(Lengths);
else
    minProp=minCustProp;maxProp=maxCustProp;
end
colresolution=50;
Cols=parula(colresolution);
refProp=linspace(minProp,maxProp,colresolution);
%color code angle
axes(handles.axes1);
% subplot(1,2,1)
hold on
for j=1:length(Data)
    col=[interp1(refProp,Cols(:,1)',max(minProp,min(maxProp,Lengths(j))))    interp1(refProp,Cols(:,2)',max(minProp,min(maxProp,Lengths(j))))    interp1(refProp,Cols(:,3)',max(minProp,min(maxProp,Lengths(j))))];
    plot3(Data(j).X,Data(j).Y,Data(j).Z,'color',col,'LineWidth',filamentThickness)
end
box on;
axis equal
rotate3d on
ax=gca;ax.XTickLabel={};ax.YTickLabel={};ax.ZTickLabel={};
% Choose default command line output for GUIPlotlam
handles.output = hObject;
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);

axes(handles.axes4)
axis off

axes(handles.axes5)
title(strrep(handles.file,'.mat',''),'FontSize',12)
text(0.01,0.8,['Length (back-front axis): ' num2str(round(handles.Length)) 'nm'],'Units','normalized','FontUnits','normalized')
text(0.01,0.5,['Height: ' num2str(round(handles.Height))  'nm'],'Units','normalized','FontUnits','normalized')
text(0.01,0.2,['Number of filaments: ' num2str(round(handles.Number))],'Units','normalized','FontUnits','normalized')
axis off

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIPlotlam wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIPlotlam_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = get(hObject, 'String');
val = get(hObject,'Value');

switch str{val};
    case 'length'
        Prop=1;
    case 'angle to the leading edge direction'
        Prop=2;
    case 'bendiness'
        Prop=3;
    case 'angle to the Z-axis'
        Prop=4;
end
handles.Prop=Prop;
guidata(hObject,handles)


% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

relColCode = (get(hObject,'Value')) ;
handles.relColCode=relColCode;
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ax=get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.shape=get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of togglebutton2


% --- Executes on button press in togglebutton3.
function togglebutton3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.colBar=get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of togglebutton3



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filamentThickness=str2num(get(hObject,'String'));
handles.filamentThickness=filamentThickness;
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

minCustProp=str2double(get(hObject,'String'));
handles.minCustProp=minCustProp;
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
maxCustProp=str2double(get(hObject,'String'));
handles.maxCustProp=maxCustProp;
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colresolution=50;
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case 'parula'
        Cols=parula(colresolution);
    case 'jet'
        Cols=jet(colresolution);
    case 'autumn'
        Cols=autumn(colresolution);
    case 'cool'
        Cols=cool(colresolution);
    case 'copper'
        Cols=copper(colresolution);
    case 'winter'
        Cols=winter(colresolution);
    case 'summer'
        Cols=summer(colresolution);
end
handles.Cols=Cols;
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1

% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  set(hObject,'Units','normalized','outerposition',[0 0 1 1])



axes(handles.axes1)
cla;
axes(handles.axes1)
cla;
str = get(hObject, 'String');
val = get(hObject,'Value');
filamentThickness=handles.filamentThickness;
relColCode=handles.relColCode;
minCustProp=handles.minCustProp;
maxCustProp=handles.maxCustProp;
Lengths=handles.Lengths;
Angle=handles.Angle;
AngletoZ=handles.AngletoZ;
Bendiness=handles.Bendiness;
% data=handles.data;
Data=handles.Data;
Cols=handles.Cols;
Prop=handles.Prop;
FilAxisX=handles.FilAxisX;
FilAxisY=handles.FilAxisY;
FilAxisZ=handles.FilAxisZ;
FilAxislen=handles.FilAxislen;
FilAxisdelX=handles.FilAxisdelX;
FilAxisdelY=handles.FilAxisdelY;
FilAxisdelZ=handles.FilAxisdelZ;
FilAxisdelLen=handles.FilAxisdelLen;
FilAxisLength=handles.FilAxisLength;
data=handles.data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D plot
switch Prop;
    case 1
        if relColCode==1
            minProp=min(Lengths);maxProp=max(Lengths);
        else
            minProp=minCustProp;maxProp=maxCustProp;
        end
        
        refProp=linspace(minProp,maxProp,length(Cols));
        %color code angle
        axes(handles.axes1)
        %         subplot(1,2,1)
        hold on
        for j=1:length(Data)
            col=[interp1(refProp,Cols(:,1)',max(minProp,min(maxProp,Lengths(j))) )    interp1(refProp,Cols(:,2)',max(minProp,min(maxProp,Lengths(j))))    interp1(refProp,Cols(:,3)',max(minProp,min(maxProp,Lengths(j))))];
            plot3(Data(j).X,Data(j).Y,Data(j).Z,'color',col,'LineWidth',filamentThickness)
        end
        box on;
        ctitle='length (nm)';
        
    case 2
        
        if relColCode==1
            minProp=min(Angle);maxProp=max(Angle);
        else
            minProp=minCustProp;maxProp=maxCustProp;
        end
        refProp=linspace(minProp,maxProp,length(Cols));
        %color code angle
        axes(handles.axes1)
        %         subplot(1,2,1)
        hold on
        for j=1:length(Data)
            col=[interp1(refProp,Cols(:,1)',max(minProp,min(maxProp,Angle(j))) )    interp1(refProp,Cols(:,2)',max(minProp,min(maxProp,Angle(j))))    interp1(refProp,Cols(:,3)',max(minProp,min(maxProp,Angle(j))))];
            plot3(Data(j).X,Data(j).Y,Data(j).Z,'color',col,'LineWidth',filamentThickness)
        end
        box on;
        ctitle='angle to the L.E. direction (deg)';
        
    case 3
        if relColCode==1
            minProp=min(Bendiness);maxProp=max(Bendiness);
        else
            minProp=minCustProp;maxProp=maxCustProp;
        end
        refProp=linspace(minProp,maxProp,length(Cols));
        %color code angle
        axes(handles.axes1)
        %         subplot(1,2,1)
        hold on
        for j=1:length(Data)
            col=[interp1(refProp,Cols(:,1)',max(minProp,min(maxProp,Bendiness(j))) )    interp1(refProp,Cols(:,2)',max(minProp,min(maxProp,Bendiness(j))))    interp1(refProp,Cols(:,3)',max(minProp,min(maxProp,Bendiness(j))))];
            plot3(Data(j).X,Data(j).Y,Data(j).Z,'color',col,'LineWidth',filamentThickness)
        end
        box on;
        ctitle='bendiness';
        
    case 4
        if relColCode==1
            minProp=min(AngletoZ);maxProp=max(AngletoZ);
        else
            minProp=minCustProp;maxProp=maxCustProp;
        end
        
        refProp=linspace(minProp,maxProp,length(Cols));
        %color code angle
        axes(handles.axes1)
        %         subplot(1,2,1)
        hold on
        for j=1:length(Data)
            col=[interp1(refProp,Cols(:,1)',max(minProp,min(maxProp,AngletoZ(j))) )    interp1(refProp,Cols(:,2)',max(minProp,min(maxProp,AngletoZ(j))))    interp1(refProp,Cols(:,3)',max(minProp,min(maxProp,AngletoZ(j))))];
            plot3(Data(j).X,Data(j).Y,Data(j).Z,'color',col,'LineWidth',filamentThickness)
        end
        box on;
        ctitle='angle to the Z-axis (deg)';
end
axis equal
rotate3d on


if handles.shape==1 %shape
    axes(handles.axes1)
    Bound=boundary(data(:,2:4));trisurf(Bound,data(:,2),data(:,3),data(:,4),'FaceColor','black','FaceAlpha',0.3,'EdgeAlpha',0.01)
end

if handles.ax==1  %axis
    axes(handles.axes1)
    plot3(FilAxisX,FilAxisY,FilAxisZ,'Linewidth',5,'color','k');
    quiver3(FilAxisX(end-1),FilAxisY(end-1),FilAxisZ(end-1),FilAxisX(end)-FilAxisX(end-3),FilAxisY(end)-FilAxisY(end-3),FilAxisZ(end)-FilAxisZ(end-3),'Linewidth',15,'color','k')
end

if handles.colBar==1  %colorbar
    ftmp=figure('visible','off');imagesc([minProp,maxProp]),ctmp=colorbar;ctmpTicks=ctmp.Ticks;ctmpTickLabels=ctmp.TickLabels;
    close(ftmp);
    axes(handles.axes1)
    colorbar('off')
    c=colorbar;
    caxis([0 1]);
    c.Ticks=interp1([minProp,maxProp],[0,1],ctmpTicks);
    c.TickLabels=ctmpTickLabels;
    colormap(Cols)
    c.Label.String=ctitle;
    c.Position=[0.6 0.5 0.015 0.4];
else
    colorbar('off')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Properties
axes(handles.axes5)
cla
text(0.01,0.8,['Length (back-front axis): ' num2str(round(handles.Length)) 'nm'],'Units','normalized','FontUnits','normalized')
text(0.01,0.5,['Height: ' num2str(round(handles.Height))  'nm'],'Units','normalized','FontUnits','normalized')
text(0.01,0.2,['Number of filaments: ' num2str(round(handles.Number))],'Units','normalized','FontUnits','normalized')
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D plot


if handles.CrsSec2==1 ||  handles.CrsSec3==1
    
    whichLen=0.03+handles.CrsSecDist*0.94;
    
    whichind=round(interp1(FilAxislen,[1:length(FilAxislen)],whichLen*FilAxisLength));
    whichX=FilAxisX(whichind);whichY=FilAxisY(whichind);whichZ=FilAxisZ(whichind);
    whichdelX=FilAxisdelX(whichind)/FilAxisdelLen(whichind);whichdelY=FilAxisdelY(whichind)/FilAxisdelLen(whichind);whichdelZ=FilAxisdelZ(whichind)/FilAxisdelLen(whichind);
    
    d = delaunayTriangulation(data(:,2), data(:,3), data(:,4));
    % d.convexHull;
    FilAxisPoint=[whichX whichY whichZ]';
    FilAxisVector=[whichdelX whichdelY whichdelZ]';
    px = sliceDelaunay(d , FilAxisVector' , FilAxisPoint');
    
    
    % filamentThickness2d=filamentThickness*10;
    colresolution=50;
    
    switch Prop;
        case 1
            if relColCode==1
                minProp=min(Lengths);maxProp=max(Lengths);
            else
                minProp=minCustProp;maxProp=maxCustProp;
            end
            refProp=linspace(minProp,maxProp,length(Cols));
            j=1;tmpLL=[];
            for i=1:length(Data)
                pointsinFil=length(Data(i).X);
                tmp=FilAxisVector(1)*(Data(i).X-FilAxisPoint(1))+FilAxisVector(2)*(Data(i).Y-FilAxisPoint(2))+FilAxisVector(3)*(Data(i).Z-FilAxisPoint(3));
                if tmp(1)*tmp(end)<0 % if the filament number i crosses the plain of interest
                    SolInd=interp1(tmp,[1:pointsinFil]',0);
                    Xcrs(j)=interp1([1:pointsinFil],Data(i).X,SolInd);
                    Ycrs(j)=interp1([1:pointsinFil],Data(i).Y,SolInd);
                    Zcrs(j)=interp1([1:pointsinFil],Data(i).Z,SolInd);
                    CrossedFilaments(j)= i;
                    tmpLL=[tmpLL Lengths(j)];
                    col2d(j,:)=[interp1(refProp,Cols(:,1)',max(minProp,min(maxProp,Lengths(i))))    interp1(refProp,Cols(:,2)',max(minProp,min(maxProp,Lengths(i))))   interp1(refProp,Cols(:,3)',max(minProp,min(maxProp,Lengths(i))))];
                    j=j+1;
                end
            end
            
        case 2
            if relColCode==1
                minProp=min(Angle);maxProp=max(Angle);
            else
                minProp=minCustProp;maxProp=maxCustProp;
            end
            refProp=linspace(minProp,maxProp,length(Cols));
            j=1;
            for i=1:length(Data)
                pointsinFil=length(Data(i).X);
                tmp=FilAxisVector(1)*(Data(i).X-FilAxisPoint(1))+FilAxisVector(2)*(Data(i).Y-FilAxisPoint(2))+FilAxisVector(3)*(Data(i).Z-FilAxisPoint(3));
                if tmp(1)*tmp(end)<0 % if the filament number i crosses the plain of interest
                    SolInd=interp1(tmp,[1:pointsinFil]',0);
                    Xcrs(j)=interp1([1:pointsinFil],Data(i).X,SolInd);
                    Ycrs(j)=interp1([1:pointsinFil],Data(i).Y,SolInd);
                    Zcrs(j)=interp1([1:pointsinFil],Data(i).Z,SolInd);
                    CrossedFilaments(j)= i;
                    col2d(j,:)=[interp1(refProp,Cols(:,1)',max(minProp,min(maxProp,Angle(i))))   interp1(refProp,Cols(:,2)',max(minProp,min(maxProp,Angle(i))))     interp1(refProp,Cols(:,3)',max(minProp,min(maxProp,Angle(i)))) ];
                    j=j+1;
                end
            end
            
        case 3
            if relColCode==1
                minProp=min(Bendiness);maxProp=max(Bendiness);
            else
                minProp=minCustProp;maxProp=maxCustProp;
            end
            refProp=linspace(minProp,maxProp,length(Cols));
            j=1;
            for i=1:length(Data)
                pointsinFil=length(Data(i).X);
                tmp=FilAxisVector(1)*(Data(i).X-FilAxisPoint(1))+FilAxisVector(2)*(Data(i).Y-FilAxisPoint(2))+FilAxisVector(3)*(Data(i).Z-FilAxisPoint(3));
                if tmp(1)*tmp(end)<0 % if the filament number i crosses the plain of interest
                    SolInd=interp1(tmp,[1:pointsinFil]',0);
                    Xcrs(j)=interp1([1:pointsinFil],Data(i).X,SolInd);
                    Ycrs(j)=interp1([1:pointsinFil],Data(i).Y,SolInd);
                    Zcrs(j)=interp1([1:pointsinFil],Data(i).Z,SolInd);
                    CrossedFilaments(j)= i;
                    col2d(j,:)=[interp1(refProp,Cols(:,1)',max(minProp,min(maxProp,Bendiness(i))))   interp1(refProp,Cols(:,2)',max(minProp,min(maxProp,Bendiness(i))))   interp1(refProp,Cols(:,3)',max(minProp,min(maxProp,Bendiness(i))))];
                    j=j+1;
                end
            end
            
        case 4
            if relColCode==1
                minProp=min(AngletoZ);maxProp=max(AngletoZ);
            else
                minProp=minCustProp;maxProp=maxCustProp;
            end
            refProp=linspace(minProp,maxProp,length(Cols));
            j=1;
            for i=1:length(Data)
                pointsinFil=length(Data(i).X);
                tmp=FilAxisVector(1)*(Data(i).X-FilAxisPoint(1))+FilAxisVector(2)*(Data(i).Y-FilAxisPoint(2))+FilAxisVector(3)*(Data(i).Z-FilAxisPoint(3));
                if tmp(1)*tmp(end)<0 % if the filament number i crosses the plain of interest
                    SolInd=interp1(tmp,[1:pointsinFil]',0);
                    Xcrs(j)=interp1([1:pointsinFil],Data(i).X,SolInd);
                    Ycrs(j)=interp1([1:pointsinFil],Data(i).Y,SolInd);
                    Zcrs(j)=interp1([1:pointsinFil],Data(i).Z,SolInd);
                    CrossedFilaments(j)= i;
                    col2d(j,:)=[interp1(refProp,Cols(:,1)',max(minProp,min(maxProp,AngletoZ(i))))   interp1(refProp,Cols(:,2)',max(minProp,min(maxProp,AngletoZ(i))))   interp1(refProp,Cols(:,3)',max(minProp,min(maxProp,AngletoZ(i))))];
                    j=j+1;
                end
            end
            
    end
    
    %Rotation
    R1=FilAxisVector; % new x is the normal vector of cross section plain
    R3=[0;0;1];
    R2=cross(R3,R1);
    R=[R1 R2 R3]';
    pxRt=R*(px-repmat(FilAxisPoint,[1,size(px,2)]));
    CoorRt=R*([Xcrs;Ycrs;Zcrs]-repmat(FilAxisPoint,[1,size(Xcrs,2)]));
    kk=convhull(CoorRt(2,:),CoorRt(3,:));
end
if handles.CrsSec2==1
    axes(handles.axes4)
    scatter(CoorRt(2,:),CoorRt(3,:),filamentThickness*25,col2d,'filled')
    patch( CoorRt(2,kk),CoorRt(3,kk),'k','FaceColor','k','FaceAlpha',0.25,'EdgeColor','none')
    axis equal
    ax=gca;ax.XTickLabel={};ax.YTickLabel={};
    axis off
else
    axes(handles.axes4)
    cla
end
if  handles.CrsSec3==1
    axes(handles.axes1)
    tmppp=(R^(-1)*CoorRt(:,kk))+repmat(FilAxisPoint,[1,size(CoorRt(:,kk),2)]);
    patch( tmppp(1,:),tmppp(2,:),tmppp(3,:),'k','FaceColor','k','FaceAlpha',0.65,'EdgeColor','none')
end



% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CrsSecDist=get(hObject,'Value');
guidata(hObject,handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fig2D=figure;
hax2d=handles.axes4;
hax_new = copyobj(hax2d, Fig2D);
set(hax_new, 'Units','normal','Position', [0.05 0.05 0.9 0.9]);
name=strrep(handles.file,'.mat','');
curntPath=pwd;
cd ..
cd('Plot_results')
tmpind=exist('Plots_of_single_cells');
if tmpind~=7
    mkdir('Plots_of_single_cells');
end
cd('Plots_of_single_cells')
print(Fig2D,[name 'CS'],'-djpeg','-r400');
print(Fig2D,[name 'CS'],'-dsvg','-r400');
close(Fig2D)
cd(curntPath)
%
% F = getframe(handles.axes1);
% Image = frame2im(F);
% imwrite(Image, 'Image.jpg')


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fig3D=figure;
hax3d=handles.axes1;
hax3dCL=handles.axes1.Colorbar;
hax_new = copyobj(hax3d, Fig3D);
hax_newCL = copyobj([hax3dCL,hax3d], Fig3D);colormap(handles.Cols);
set(hax_new, 'Units','normal','Position', [0.02 0.05 0.82 0.9]);
pause(0.2)
set(hax_newCL, 'Units','normal','Position', [0.86 0.3 0.02 0.4],'YAxisLocation','right');
pause(0.2)
name=strrep(handles.file,'.mat','');
curntPath=pwd;
cd ..
cd('Plot_results')
tmpind=exist('Plots_of_single_cells');
if tmpind~=7
    mkdir('Plots_of_single_cells');
end
cd('Plots_of_single_cells')
pause(0.5);
print(Fig3D,[name '3D'],'-djpeg','-r400');
print(Fig3D,[name '3D'],'-dsvg','-r400');
close(Fig3D)
cd(curntPath)
%


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CrsSec3=get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CrsSec2=get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton2
