function varargout = Plots_Group_Cells_Filopodia(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Plots_Group_Cells_Filopodia_OpeningFcn, ...
    'gui_OutputFcn',  @Plots_Group_Cells_Filopodia_OutputFcn, ...
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


% --- Executes just before Plots_Group_Cells_Filopodia is made visible.
function Plots_Group_Cells_Filopodia_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Plots_Group_Cells_Filopodia (see VARARGIN)

% Choose default command line output for Plots_Group_Cells_Filopodia
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Plots_Group_Cells_Filopodia wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Plots_Group_Cells_Filopodia_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
warning off
addpath('files')
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
cd(folder)
tmpind=exist('DataSummary');
if tmpind~=7
    error('Please run the CellGrouping-Filopodia script');
end
ResultsPath=fullfile(folder,'DataSummary');
cd(ResultsPath);
groups=dir('*mat');
cd(currentfolder)

handles.groups=groups;
handles.ResultsPath=ResultsPath;
handles.folder=folder;
handles.currentfolder=currentfolder;

% Update handles structure
guidata(hObject, handles);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parallel coordinates
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
groups=handles.groups;
ResultsPath=handles.ResultsPath;
folder=handles.ResultsPath;
currentfolder=handles.currentfolder;
LabelsDefFil;
cd(ResultsPath);
if numel(groups)>7
    colors=jet(numel(groups));
else
    colors=get(gca,'colororder');
end
% close all
Data=[];
for i=1:length(groups)
    load(groups(i).name)
    GroupsNAmes{i}=strrep(strrep(strrep(groups(i).name,'.mat',''),'FilopMetrics',''),'_','-');
    clear data
    data(:,1)=[FilopMetrics.FilamentNumber]';
    data(:,2)=[FilopMetrics.Density]';
    data(:,3)=[FilopMetrics.Length]';
    data(:,4)=[FilopMetrics.Bendiness]';
    data(:,5)=[FilopMetrics.BendingEnergy]';
    data(:,6)=[FilopMetrics.FilamentAreaAVG]';
    data(:,7)=[FilopMetrics.FilamentDensityAVG]';
    data(:,8)=[FilopMetrics.FilamentVolFrAVG]';
    data(:,9)=[FilopMetrics.FilamentAreaNumAVG]';
    data(:,10)=[FilopMetrics.SegRoundnessAVG]';
    data(:,11)=[FilopMetrics.SegIxAVG]';
    data(:,12)=[FilopMetrics.SegIzAVG]';
    data(:,13)=[FilopMetrics.FilamentAnisotropy]';
    data(:,14)=[FilopMetrics.FilamentLengthAVG]';
    data(:,15)=[FilopMetrics.FilamentAngleAVG]';
    data(:,16)=[FilopMetrics.FilamentBendinessAVG]';
    data(:,17)=[FilopMetrics.FilamentBendingEnergyAVG]';
    data(:,18)=[FilopMetrics.FilamentAreaBtT]';
    data(:,19)=[FilopMetrics.FilamentDensityBtT]';
    data(:,20)=[FilopMetrics.FilamentAreaNumBtT]';
    data(:,21)=[FilopMetrics.SegIxBtT]';
    data(:,22)=[FilopMetrics.SegIzBtT]';
    data(:,23)=[FilopMetrics.FilamentSeglengthBtT]';
    data(:,24)=[FilopMetrics.FilamentSegAngleBtT]';
    data(:,25)=[FilopMetrics.FilamentSegBendinessBtT]';
    data(:,26)=[FilopMetrics.FilamentSegBendingEnergyBtT]';
    NumofFilop(i)=size(data,1);
    Data=[Data;data];
end
% Select desired parameters from a list
list=Labels;
[indx,tf]=listdlg('ListString',list,'ListSize',[300,400]);
wb = waitbar(0,'please wait...');
if tf==1
    Labels=Labels(indx);
    Data=Data(:,indx);
end
%
Data=zscore(Data);
figure
hold on
x=-[1:(size(Data,2)+1)];
xx=x(1:end-1);
sumNumofFilop=cumsum(NumofFilop);
data=[];
for i=1:length(groups)
    if i==1
        data=Data(1:NumofFilop(i),:);
    elseif i>1
        data=Data(sumNumofFilop(i-1)+1:sumNumofFilop(i),:);
    end
    for j=1:NumofFilop(i)
        tmp=zeros(size(x));
        tmp(1:length(x)-1)=data(j,1:end);
        tmp(end)=NaN;
        patch(tmp,x,colors(i,:),'EdgeColor',colors(i,:),...
            'FaceVertexAlphaData',0.2*ones(size(x))','AlphaDataMapping','none',...
            'EdgeAlpha','interp')
    end
    if size(data,1)==1
        errBar=zeros(size(data));
    else
        errBar=nanstd(data,1)/sqrt(NumofFilop(i));
    end
    e=errorbar(nanmean(data,1),xx+0.001,errBar,'o','horizontal');
    e.Color = colors(i,:);
    e.LineWidth = 1.5;
    e.MarkerSize=3;
    
    md(i)=plot(nanmean(data,1),xx,'-o','color',colors(i,:),...
        'LineWidth',2,...
        'MarkerSize',5,...
        'MarkerEdgeColor',colors(i,:),...
        'MarkerFaceColor',[0.5,0.5,0.5]);
end
axis tight
ax=gca;
ax.YTickMode = 'manual';
ax.YTick = fliplr(xx);
ax.YLim=[x(end) 0];
ax.YTickLabel = fliplr(Labels);
xlabel('Normalized parameters (centered to mean=0 and SD=1)')
set(gcf,'units','normalized', 'Position',  [0.1 0.1 0.5 0.8])
legend(md,GroupsNAmes)
% corrolations
for ij=1:length(groups)
    if ij==1
        data=Data(1:NumofFilop(ij),:);
    elseif ij>1
        data=Data(sumNumofFilop(ij-1)+1:sumNumofFilop(ij),:);
    end
    [R,P]=corrcoef(data);
    %%%% polygon correlation diagram
    TR=0.5; %Threshold above which the correlations appears in these graphs
    figure
    COR=R;
    rrr=[1:length(R)]*2*pi/length(R);
    polarplot(rrr,15*ones(1,length(R)),'ko','MarkerFaceColor','k')
    hold on
    for i=1:18
        for j=i+1:length(R)
            if abs(P(j,i))<0.02
                h=polarplot([rrr(i) rrr(j)],15*ones(1,2));
                set(h,'color',[max(0.1,-sign(COR(j,i))) 0.5 max(0.1,sign(COR(j,i)))],'linewidth',(max(0.1,10*(abs(COR(j,i))-0.6))))
            end
        end
    end
    ax=gca;
    ax.ThetaTick = rrr*180/pi;
    ax.ThetaTickLabelMode='manual';
    ax.ThetaTickLabel = Labels;
    ax.RTick=[];
    set(gcf,'units','normalized', 'Position',  [0.1 0.1 0.8 0.8])
    title(GroupsNAmes(ij));
    %
end
close(wb);
cd(currentfolder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bar charts
function pushbutton7_Callback(hObject, eventdata, handles)
groups=handles.groups;
ResultsPath=handles.ResultsPath;
folder=handles.ResultsPath;
currentfolder=handles.currentfolder;
LabelsDefFilUnits;
cd(ResultsPath);
[indxs,tf] = listdlg('PromptString','Select parameters for bar chart',...
    'SelectionMode','single','ListString',Labels,'ListSize',[350,450]);
wb = waitbar(0,'please wait...');
% select the desired parameters
Data=[];
Databar=[];Databarer=[];
for i=1:length(groups)
    clear data
    load(groups(i).name)
    XLabels{i}=strrep(strrep(strrep(groups(i).name,'.mat',''),'FilopMetrics',''),'_','-');
    data(:,1)=[FilopMetrics.FilamentNumber]';
    data(:,2)=[FilopMetrics.Density]';
    data(:,3)=[FilopMetrics.Length]';
    data(:,4)=[FilopMetrics.Bendiness]';
    data(:,5)=[FilopMetrics.BendingEnergy]';
    data(:,6)=[FilopMetrics.FilamentAreaAVG]';
    data(:,7)=[FilopMetrics.FilamentDensityAVG]';
    data(:,8)=[FilopMetrics.FilamentVolFrAVG]';
    data(:,9)=[FilopMetrics.FilamentAreaNumAVG]';
    data(:,10)=[FilopMetrics.SegRoundnessAVG]';
    data(:,11)=[FilopMetrics.SegIxAVG]';
    data(:,12)=[FilopMetrics.SegIzAVG]';
    data(:,13)=[FilopMetrics.FilamentAnisotropy]';
    data(:,14)=[FilopMetrics.FilamentLengthAVG]';
    data(:,15)=[FilopMetrics.FilamentAngleAVG]';
    data(:,16)=[FilopMetrics.FilamentBendinessAVG]';
    data(:,17)=[FilopMetrics.FilamentBendingEnergyAVG]';
    data(:,18)=[FilopMetrics.FilamentAreaBtT]';
    data(:,19)=[FilopMetrics.FilamentDensityBtT]';
    data(:,20)=[FilopMetrics.FilamentAreaNumBtT]';
    data(:,21)=[FilopMetrics.SegIxBtT]';
    data(:,22)=[FilopMetrics.SegIzBtT]';
    data(:,23)=[FilopMetrics.FilamentSeglengthBtT]';
    data(:,24)=[FilopMetrics.FilamentSegAngleBtT]';
    data(:,25)=[FilopMetrics.FilamentSegBendinessBtT]';
    data(:,26)=[FilopMetrics.FilamentSegBendingEnergyBtT]';
    databar=[nanmean(data(:,indxs))]';
    Databar=[Databar databar];
    databarer=[nanstd(data(:,indxs))]'/sqrt(size(data,1));
    Databarer=[Databarer databarer];
end
figure
if numel(groups)>7
    colors=jet(numel(groups));
else
    colors=get(gca,'colororder');
end
hold on
for i = 1:length(Databar)
    bb(i)=bar(i,Databar(i),0.7);
    set(bb(i),'FaceColor',colors(i,:));
end
ee=errorbar([1:length(Databar)],Databar,Databarer,'.','color','k','LineWidth',1.5);
ax=gca;ax.FontSize=14; ylabel( Labels(indxs)) ;
ax.XTickMode = 'manual'; ax.XTick =[1:length(Databar)]; ax.XTickLabel = XLabels;
% ax.LabelFontSizeMultiplier=1.3;ax.FontName='Times New Roman';
ax.XLim=[0 length(Databar)+1];
box on
axf=gcf;
set(gcf, 'Units', 'normalized');
set(gcf, 'Position',[0.2 0.2 0.6 0.6]);
close(wb);
cd(currentfolder)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histograms
function pushbutton8_Callback(hObject, eventdata, handles)
groups=handles.groups;
ResultsPath=handles.ResultsPath;
folder=handles.ResultsPath;
currentfolder=handles.currentfolder;
list={'Length of filaments','Angle of filaments','Bendiness of filaments','Bending energy density of filaments','Z-Angle of filaments'};
whichHist=listdlg('PromptString','Plot the histograms of:','ListString',list,'SelectionMode','single');
if whichHist==1
    definputs2={'40','40','600'};
elseif whichHist==2
    definputs2={'1','2','50'};
elseif whichHist==3
    definputs2={'1.01','0.005','1.07'};
elseif whichHist==4
    definputs2={'0.000002','0.000005','0.00005'};
elseif whichHist==5
    definputs2={'1','5','89'};
end
prompt2={'First Edge of the binning','Width of the binning','Last Edge of the binning'}; % filtering out filaments
inputs=inputdlg(prompt2,'Inputs',[1 35],definputs2);
wb = waitbar(0,'please wait...');
firstedge=str2num(inputs{1});steps=str2num(inputs{2});lastedge=str2num(inputs{3});
cd(ResultsPath);
figure
hold on;
edges=[0,firstedge:steps:lastedge,inf];
Xs=1:2:2*(size(edges,2)-1);
if numel(groups)>7
    colors=jet(numel(groups));
else
    colors=get(gca,'colororder');
end
barwidth=0.4/length(groups);
barwidthExt=1.4/length(groups);
for i=1:length(groups)
    load(groups(i).name)
    GroupsNAmes{i}=strrep(strrep(strrep(groups(i).name,'.mat',''),'FilopMetrics',''),'_','-');
    clear hpercentage
    for j=1:length(FilopMetrics)
        if whichHist==1
            N=histcounts([FilopMetrics(j).FilamentLength],edges);
        elseif whichHist==2
            N=histcounts([FilopMetrics(j).FilamentAngle],edges);
        elseif whichHist==3
            N=histcounts([FilopMetrics(j).FilamentBendiness],edges);
        elseif whichHist==4
            N=histcounts([FilopMetrics(j).FilamentBendingEnergy],edges);
        elseif whichHist==5
            N=histcounts([FilopMetrics(j).FilamentAngletoZ],edges);
        end
        hpercentage(j,:)=100*N/sum(N);
    end
    Hpercentage(i,:)=nanmean(hpercentage,1);
    errors(i,:)=nanstd(hpercentage,1)/sqrt(length(FilopMetrics));
    b(i)=bar(Xs+(i-1)*barwidthExt,Hpercentage(i,:),barwidth,'FaceColor',colors(i,:),'EdgeColor',colors(i,:));
    %     er = errorbar(Xs+(i-1)*0.3,Hpercentage(i,:),errors(i,:),'Color','k', 'LineStyle','none','CapSize',3);
    er = errorbar(Xs+(i-1)*barwidthExt,Hpercentage(i,:),errors(i,:),'Color','k', 'LineStyle','none');
end
legend(b,GroupsNAmes)
Edges=(edges(1:end-1)+edges(2:end))/2;
Labels{1}=['<' num2str(firstedge)];
for k=2:length(edges)-2
    Labels{k}=num2str(edges(k));
end
Labels{length(edges)-1}=['>' num2str(edges(end-1))];
ax=gca; ax.XTick=Xs+(length(groups)-1)*0.3; ax.XTickLabel=Labels;
ax.YLim(1)=0;
list={'Length of filaments (nm)','Angle of filaments (deg)','Bendiness of filaments','Bending energy density of filaments (nm^{-2})','Angle of filaments to Z-axis (deg)'};
xlabel(list(whichHist));
ylabel('% of filaments');
ax=gca;
ax.XAxis.TickLength=[0 0];
box on
axf=gcf;axf.Units='normal'; axf.Position=[0.1 0.3 0.8 0.4];
close(wb);
cd(currentfolder);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Properties along axis
function pushbutton9_Callback(hObject, eventdata, handles)
groups=handles.groups;
ResultsPath=handles.ResultsPath;
folder=handles.ResultsPath;
currentfolder=handles.currentfolder;
list={'Length of filaments','Angle of filaments','Bendiness of filaments','Bending energy density of filaments','Cross-sectional area','Cross-sectional density','Cross-sectional volume fraction','Cross-sectional number of filaments',...
    'Cross-sectional circularity','Lateral bending stiffness','Vertical bending stiffness','Barbed/Pointed Ends'};
whichHist=listdlg('PromptString','Select the property to plot in dependence on the distance along the axis of the structure:','ListString',list,'SelectionMode','single','ListSize',[450 250]);
if whichHist<12
    dlgtitle='Edge margin';
    prompt={'Removed margin around edges of the structure (% of the structure length)','Number of cross-sections along the structures'};
    defaultanswer={'10','50'};
    anss=inputdlg(prompt, dlgtitle, [1 75], defaultanswer );
    wb = waitbar(0,'please wait...');
    RemoveMargin=str2num(anss{1});
    numOfCrs=str2num(anss{2});
    figure
    fignum=get(gcf,'Number');
    if numel(groups)>7
        colors=jet(numel(groups));
    else
        colors=get(gca,'colororder');
    end
    hold on;
    options.handle     = figure(fignum);
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'sem';
    for i=1:length(groups)
        cd(ResultsPath);
        load(groups(i).name)
        GroupsNAmes{i}=strrep(strrep(strrep(groups(i).name,'.mat',''),'FilopMetrics',''),'_','-');
        cd(currentfolder)
        options.color_line=colors(i,:);
        options.color_area=1-(1-colors(i,:))/2.5;
        NumofSections=length(FilopMetrics(1).WhichLen);
        reshapedData=[];reshapedWhichLen=[];
        switch whichHist
            case 1
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).Seglength)],FilopMetrics(k).Seglength,linspace(1,numel(FilopMetrics(k).Seglength),numOfCrs));
                end
            case 2
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).SegAngle)],FilopMetrics(k).SegAngle,linspace(1,numel(FilopMetrics(k).SegAngle),numOfCrs));
                end
            case 3
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).SegBendiness)],FilopMetrics(k).SegBendiness,linspace(1,numel(FilopMetrics(k).SegBendiness),numOfCrs));
                end
            case 4
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).SegBendingEnergy)],FilopMetrics(k).SegBendingEnergy,linspace(1,numel(FilopMetrics(k).SegBendingEnergy),numOfCrs));
                end
            case 5
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).SegArea)],FilopMetrics(k).SegArea,linspace(1,numel(FilopMetrics(k).SegArea),numOfCrs));
                end
            case 6
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).SegDensity)],FilopMetrics(k).SegDensity,linspace(1,numel(FilopMetrics(k).SegDensity),numOfCrs));
                end
            case 7
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).SegVolFraction)],FilopMetrics(k).SegVolFraction,linspace(1,numel(FilopMetrics(k).SegVolFraction),numOfCrs));
                end
            case 8
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).SegFilNum)],FilopMetrics(k).SegFilNum,linspace(1,numel(FilopMetrics(k).SegFilNum),numOfCrs));
                end
            case 9
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).SegRoundness)],FilopMetrics(k).SegRoundness,linspace(1,numel(FilopMetrics(k).SegRoundness),numOfCrs));
                end
            case 10
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).Ix)],FilopMetrics(k).Ix,linspace(1,numel(FilopMetrics(k).Ix),numOfCrs));
                end
            case 11
                for k=1:numel(FilopMetrics)
                    reshapedData(k,:)=interp1([1:numel(FilopMetrics(k).Iz)],FilopMetrics(k).Iz,linspace(1,numel(FilopMetrics(k).Iz),numOfCrs));
                end
        end
        for k=1:numel(FilopMetrics)
            reshapedWhichLen(k,:)=interp1([1:numel(FilopMetrics(k).WhichLen)],FilopMetrics(k).WhichLen,linspace(1,numel(FilopMetrics(k).WhichLen),numOfCrs));
        end
        reshapedData(:,1:round(numOfCrs/100*RemoveMargin))=[];    reshapedData(:,size(reshapedData,2)-round(numOfCrs/100*RemoveMargin)+1:end)=[];
        reshapedWhichLen(:,1:round(numOfCrs/100*RemoveMargin))=[];reshapedWhichLen(:,size(reshapedWhichLen,2)-round(numOfCrs/100*RemoveMargin)+1:end)=[];
        hold on
        p(i)=plot_areaerrorbar(reshapedData,options);
        hold on
    end
    figtmp=figure('visible','off');
    plot(reshapedWhichLen,ones(size(reshapedWhichLen)));ax=gca;
    ax.XLim=[reshapedWhichLen(1) reshapedWhichLen(end)];
    xtick=ax.XTick; xticklabel=ax.XTickLabels;close(figtmp)
    figure(fignum)
    ax=gca;
    ax.XTick=interp1([reshapedWhichLen(1) reshapedWhichLen(end)],[1 length(reshapedWhichLen)],xtick);
    ax.XTickLabel=xticklabel;
    xlabel('Normalized distance along the axis of the structure')
    list={'Length of filaments (nm)','Angle of filaments (deg)','Bendiness of filaments','Bending energy density of filaments (nm^{-2})','Cross-sectional area (nm^{2})','Cross-sectional density (nm^{-2})','Cross-sectional volume fraction','Cross-sectional number of filaments',...
        'Cross-sectional circularity','Lateral bending stiffness (nm^{2})','Vertical bending stiffness (nm^{2})'};
    ylabel(list(whichHist))
    cd(currentfolder);
    ax.XLim=[1 length(reshapedWhichLen)];
    box on
    legend(p,GroupsNAmes)
    
else % plotting barbed/pointed ends
    dlgtitle='Inputs';
    prompt={'Removed margin around edges of the structure (% of the structure length)','Number of bins'};
    defaultanswer={'5','10'};
    anss=inputdlg(prompt, dlgtitle, [1 75], defaultanswer );
    RemoveMargin=str2num(anss{1});nbins=str2num(anss{2}); edges=linspace((RemoveMargin)/100,1-((RemoveMargin+1)/100),nbins+1);
    
    figtmp=figure('visible','off');
    plot(edges,ones(size(edges)));ax=gca;
    ax.XLim=[edges(1)+0.01 edges(end)];axis tight;
    xtick=ax.XTick; xticklabel=ax.XTickLabels;close(figtmp)
    
    distancetype=questdlg('Plot filament ends distribution in dependence on:','Distance','normalized distance along axis','real distance from base','real distance from tip','normalized distance along axis');
    wb = waitbar(0,'please wait...');
    cols=lines(2);colsPale=[0.55 0.8 0.9;0.9 0.75 0.65];
    for i=1:length(groups)
        cd(ResultsPath);
        load(groups(i).name)
        GroupsNAmes{i}=strrep(strrep(strrep(groups(i).name,'.mat',''),'FilopMetrics',''),'_','-');
        for k=1:numel(FilopMetrics)
            SEGarea(k,:)=interp1([1:numel(FilopMetrics(k).SegArea)],FilopMetrics(k).SegArea,linspace(1,numel(FilopMetrics(k).SegArea),100));
        end
        cd(currentfolder)
        figure;hold on
        switch distancetype
            case 'normalized distance along axis'
                for j=1:numel(FilopMetrics)
                    tmpB=FilopMetrics(j).FilamentBarbedEndsNorm;  tmpP=FilopMetrics(j).FilamentPointedEndsNorm;
                    tmpB(tmpB<RemoveMargin/100)=[]; tmpB(tmpB>1-((RemoveMargin+1)/100))=[];
                    tmpP(tmpP<RemoveMargin/100)=[]; tmpP(tmpP>1-((RemoveMargin+1)/100))=[];
                    lengthofStructure=sum(FilopMetrics(j).FilamentBarbedEnds)/sum(FilopMetrics(j).FilamentBarbedEndsNorm);
                    InterpArea=interp1([1:size(SEGarea,2)],SEGarea(j,:),linspace(1,size(SEGarea,2),nbins+1));
                    InterpVol=(InterpArea(1:end-1)+InterpArea(2:end))/2* lengthofStructure/nbins;
                    NBEnd(j,:)=10^9*histcounts(tmpB,edges)./InterpVol;
                    NPEnd(j,:)=10^9*histcounts(tmpP,edges)./InterpVol;
                    plot( NBEnd(j,:),'color',colsPale(1,:),'linewidth',0.1);
                    plot( NPEnd(j,:),'color',colsPale(2,:),'linewidth',0.1);
                end
                pB = plot(nanmean(NBEnd),'-o','color',cols(1,:),'linewidth',1.5,'MarkerFaceColor',cols(1,:));
                pP = plot(nanmean(NPEnd),'-o','color',cols(2,:),'linewidth',1.5,'MarkerFaceColor',cols(2,:));
                title(GroupsNAmes{i})
                legend([pB ,pP],{'Barbed ends','Pointed ends'})
                ylabel(['Density of filament ends (' char(181) 'm^{-3})'])
                xlabel('Normalized distance along the axis of the structure')
                ax=gca;ax.XLim=[1 length(nanmean(NBEnd))];
                ax.XTick=interp1( [edges(1) edges(end)],[1 length(nanmean(NBEnd))],xtick);
                ax.XTickLabel=xticklabel;
                box on
            case 'real distance from base'
                for j=1:numel(FilopMetrics)
                    ttmpB=FilopMetrics(j).FilamentBarbedEndsNorm;  ttmpP=FilopMetrics(j).FilamentPointedEndsNorm;
                    tttmpB=find(ttmpB<RemoveMargin/100 | ttmpB>1-((RemoveMargin+1)/100));
                    tttmpP=find(ttmpP<RemoveMargin/100 | ttmpP>1-((RemoveMargin+1)/100));
                    tmpB=FilopMetrics(j).FilamentBarbedEnds;  tmpP=FilopMetrics(j).FilamentPointedEnds;
                    tmpB(tttmpB)=[];tmpP(tttmpP)=[];
                    lengthofStructure=sum(FilopMetrics(j).FilamentBarbedEnds)/sum(FilopMetrics(j).FilamentBarbedEndsNorm);
                    InterpArea=interp1([1:size(SEGarea,2)],SEGarea(j,:),linspace(1,size(SEGarea,2),nbins+1));
                    InterpVol=(InterpArea(1:end-1)+InterpArea(2:end))/2* lengthofStructure/nbins;
                    Ends(j).B=10^9*histcounts(tmpB,nbins)./InterpVol;
                    Ends(j).P=10^9*histcounts(tmpP,nbins)./InterpVol;
                    Ends(j).len=linspace(0,lengthofStructure*(1-2*RemoveMargin/100),nbins);
                    plot( Ends(j).len,Ends(j).B,'color',colsPale(1,:),'linewidth',0.1);
                    plot( Ends(j).len,Ends(j).P,'color',colsPale(2,:),'linewidth',0.1);
                end
                for j=1:numel(FilopMetrics)
                    interpEndB(j,:)=interp1(Ends(j).len,Ends(j).B,linspace(0,max([Ends.len]),nbins));
                    interpEndP(j,:)=interp1(Ends(j).len,Ends(j).P,linspace(0,max([Ends.len]),nbins));
                end
                pB = plot(linspace(0,max([Ends.len]),nbins),nanmean(interpEndB),'-o','color',cols(1,:),'linewidth',1.5,'MarkerFaceColor',cols(1,:));
                pP = plot(linspace(0,max([Ends.len]),nbins),nanmean(interpEndP),'-o','color',cols(2,:),'linewidth',1.5,'MarkerFaceColor',cols(2,:));
                title(GroupsNAmes{i})
                legend([pB ,pP],{'Barbed ends','Pointed ends'})
                ylabel(['Density of filament ends (' char(181) 'm^{-3})'])
                xlabel('Distance from the base along the axis of the structure (nm)')
                ax=gca;ax.XLim=[0 max([Ends.len])];
                box on
                
            case 'real distance from tip'
                for j=1:numel(FilopMetrics)
                    ttmpB=FilopMetrics(j).FilamentBarbedEndsNorm;  ttmpP=FilopMetrics(j).FilamentPointedEndsNorm;
                    tttmpB=find(ttmpB<RemoveMargin/100 | ttmpB>1-((RemoveMargin+1)/100));
                    tttmpP=find(ttmpP<RemoveMargin/100 | ttmpP>1-((RemoveMargin+1)/100));
                    tmpB=FilopMetrics(j).FilamentBarbedEnds;  tmpP=FilopMetrics(j).FilamentPointedEnds;
                    tmpB(tttmpB)=[];tmpP(tttmpP)=[];
                    
                    lengthofStructure=sum(FilopMetrics(j).FilamentBarbedEnds)/sum(FilopMetrics(j).FilamentBarbedEndsNorm);
                    InterpArea=interp1([1:size(SEGarea,2)],SEGarea(j,:),linspace(1,size(SEGarea,2),nbins+1));
                    InterpVol=(InterpArea(1:end-1)+InterpArea(2:end))/2* lengthofStructure/nbins;
                    Ends(j).B=10^9*histcounts(tmpB,nbins)./InterpVol;
                    Ends(j).P=10^9*histcounts(tmpP,nbins)./InterpVol;
                    Ends(j).len=linspace(-lengthofStructure*(1-2*RemoveMargin/100),0,nbins);
                    plot( -Ends(j).len,Ends(j).B,'color',colsPale(1,:),'linewidth',0.1);
                    plot( -Ends(j).len,Ends(j).P,'color',colsPale(2,:),'linewidth',0.1);
                end
                for j=1:numel(FilopMetrics)
                    interpEndB(j,:)=interp1(Ends(j).len,Ends(j).B,linspace(0,min([Ends.len]),nbins));
                    interpEndP(j,:)=interp1(Ends(j).len,Ends(j).P,linspace(0,min([Ends.len]),nbins));
                end
                pB = plot(-linspace(0,min([Ends.len]),nbins),nanmean(interpEndB),'-o','color',cols(1,:),'linewidth',1.5,'MarkerFaceColor',cols(1,:));
                pP = plot(-linspace(0,min([Ends.len]),nbins),nanmean(interpEndP),'-o','color',cols(2,:),'linewidth',1.5,'MarkerFaceColor',cols(2,:));
                title(GroupsNAmes{i})
                legend([pB ,pP],{'Barbed ends','Pointed ends'})
                ylabel(['Density of filament ends (' char(181) 'm^{-3})'])
                xlabel('Distance from the tip along the axis of the structure (nm)')
                ax=gca;ax.XLim=[  0 -min([Ends.len])];
                box on
        end
    end
end
close(wb)
cd(currentfolder)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration of filaments
function pushbutton10_Callback(hObject, eventdata, handles)
groups=handles.groups;
ResultsPath=handles.ResultsPath;
folder=handles.ResultsPath;
currentfolder=handles.currentfolder;
prompt2={'Min distance between filaments (nm)','Max distance between filaments (nm)','Min interfilament angle (deg)','Max interfilament angle (deg)','Number of bins for distance','Number of bins for angle' }; % filtering out filaments
definputs={'0','300','0','90','30','15'};
filterings=inputdlg(prompt2,'Inputs for the intensity map',[1 70],definputs);
minDist=str2num(filterings{1});
maxDist=str2num(filterings{2});
minAng=str2num(filterings{3});
maxAng=str2num(filterings{4});
NbinsDist=str2num(filterings{5});
NbinsAng=str2num(filterings{6});
wb = waitbar(0,'please wait...');
for i=1:length(groups)
    clear FilopMetrics
    cd(ResultsPath);
    load(groups(i).name)
    GroupsNAmes{i}=strrep(strrep(strrep(groups(i).name,'.mat',''),'FilopMetrics',''),'_','-');
    cd(currentfolder)
    AllAngles=[]; AllDistance=[];
    for ii=1:length(FilopMetrics)
        AllAngles=[AllAngles ; FilopMetrics(ii).PairwiseFilsAngles];
        AllDistance=[ AllDistance ; FilopMetrics(ii).PairwiseFilsDistances];
    end
    Ntmp=hist3([AllDistance AllAngles ],'Edges',{linspace(minDist,maxDist,NbinsDist) linspace(minAng,maxAng,NbinsAng)});
    NN=reshape(Ntmp',[NbinsAng,NbinsDist]);
    
    figure
    pc=pcolor(linspace(minDist,maxDist,NbinsDist),linspace(minAng,maxAng,NbinsAng),NN);
    pc.EdgeColor='none';
    pc.FaceColor='interp';
    ylabel('relative orientation (deg)');
    xlabel('interfilament distance (nm)')
    title(GroupsNAmes{i});
    box on
    figure
    histogram(AllAngles,NbinsAng,'Normalization','probability');
    text(0.8,0.85,['mean=' num2str(round(10*mean(AllAngles))/10)],'Units','normalized')
    text(0.8,0.92,['median=' num2str(round(10*median(AllAngles))/10)],'Units','normalized')
    xlabel('relative orientation (deg)');
    ylabel('probability');
    title(GroupsNAmes{i});
    box on
    figure
    histogram(AllDistance,NbinsDist,'Normalization','probability');
    text(0.8,0.85,['mean=' num2str(round(10*mean(AllDistance))/10)],'Units','normalized')
    text(0.8,0.92,['median=' num2str(round(10*median(AllDistance))/10)],'Units','normalized')
    xlabel('interfilament distance (nm)');
    ylabel('probability');
    title(GroupsNAmes{i});
    box on
    
    figure
    filteroutliers=find(AllDistance>maxDist | AllDistance<minDist | AllAngles>maxAng | AllAngles<minAng);
    filteredDist=AllDistance;filteredDist(filteroutliers)=[]; filteredAngle=AllAngles;filteredAngle(filteroutliers)=[];
    a=scatterhist(filteredDist,filteredAngle,'kernel','on','Location','NorthEast','Direction','Out','Color','k','Marker','.');
    a(1).XAxisLocation='bottom';a(1).YAxisLocation='left';
    a(1).XLim=[minDist maxDist]; a(1).YLim=[minAng maxAng];
    a(1).Children.Color='none';
    a(1).YLabel.String='relative orientation (deg)';
    a(1).XLabel.String='interfilament distance (nm)';
    copyobj(pc,a(1))
    a(2).Position(2)=a(1).Position(2)+a(1).Position(4)+0.01;
    a(3).Position(1)=a(1).Position(1)+a(1).Position(3)+0.015;
    mybox=uicontrol('style','text');
    pause(0.1)
    set(mybox,'string',GroupsNAmes{i},'Position',[a(3).Position(1) a(2).Position(2) 0.2 0.2],'Units','Normalized','FontSize',13);
    set(mybox,'string',GroupsNAmes{i},'Position',[a(3).Position(1) a(2).Position(2) 0.2 0.2],'Units','Normalized','FontSize',13);
    pause(0.1)
    waitbar((i)/length(groups),wb);
end
close(wb);
cd(currentfolder)
