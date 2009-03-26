function handles = findEstimates(handles)
l=0; %index of structure TGA
lines = handles.var.lines;
handles.E.checkbox=[];
handles.E.a = 1; %default
handles.E.line = 0; %line number of chosen TGA data
for i=1:lines
   if strcmp(handles.EXPDATA(i).type,'TGA')
      l=l+1;
      handles.E.TGA(l).rate = handles.EXPDATA(i).rate;
      handles.E.TGA(l).temperature = handles.EXPDATA(i).temperature;
      handles.E.TGA(l).time = handles.EXPDATA(i).time;
      handles.E.TGA(l).TGA = handles.EXPDATA(i).TGA;
      handles.E.TGA(l).material = handles.EXPDATA(i).material;
      handles.E.TGA(l).gas = handles.EXPDATA(i).gas;
      handles.E.TGA(l).name = handles.EXPDATA(i).name;
      handles.E.TGA(l).gradient = handles.EXPDATA(i).gradient;
      handles.E.TGA(l).check = 0;
      handles.E.TGA(l).estimates = [0 0];
      handles.E.TGA(l).n = 1;
   end
end
handles.TGA.hTGAEstimates = figure('Visible','off',...
   'Name', 'Find Estimates', ...
   'NumberTitle', 'off', ...
   'Menubar', 'none', ...
   'Toolbar', 'none', ...
   'Color', get(0,...
   'defaultuicontrolbackgroundcolor'), ...
   'Units', 'normalized', ...
   'Position',[0.36,0.3,0.4,0.7]);
hTGAEstimates = handles.TGA.hTGAEstimates;
movegui(hTGAEstimates,'center')

htitleInt1 = uicontrol(hTGAEstimates,'Style','text',...
   'String','MATERIAL',...
   'Units', 'normalized', ...
   'Position',[0.07 0.93 0.15 0.03]);
htitleInt2 = uicontrol(hTGAEstimates,'Style','text',...
   'String','RATE',...
   'Units', 'normalized', ...
   'Position',[0.24 0.93 0.08 0.03]);
htitleInt3 = uicontrol(hTGAEstimates,'Style','text',...
   'String','GAS',...
   'Units', 'normalized', ...
   'Position',[0.32 0.93 0.08 0.03]);
htitleInt4 = uicontrol(hTGAEstimates, 'Style', 'text', ...
   'string', 'Select', ...
   'Units', 'normalized', ...
   'Position', [0.40 0.93 0.14 0.03]);

% Create the button group.
h = uibuttongroup('visible','off','Position',[0.55 0.95 0.42 0.04]);
% Create three radio buttons in the button group.
u1 = uicontrol('Style','Radio','String','a=m/m(0)',...
    'Units', 'normalized', ...
    'pos',[0.01 0.4 0.4 0.5],'parent',h,'HandleVisibility','off');
u2 = uicontrol('Style','Radio','String','a=(m-(inf))/m(0)',...
    'Units', 'normalized', ...
    'pos',[0.4 0.4 0.55 0.5],'parent',h,'HandleVisibility','off');
% Initialize some button group properties. 
handles.E.radio = [u1 u2];
set(h,'SelectionChangeFcn',@alfa_cb);
set(h,'Visible','on');

hfindE=uicontrol(hTGAEstimates, 'Style', 'pushbutton', ...
   'String', 'Find Estimates', ...
   'Callback', @findE_cb, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.73 0.2 0.04]);
hshowE =uicontrol(hTGAEstimates, 'Style', 'pushbutton', ...
   'String', 'Show Estimates', ...
   'Callback', @showE_cb, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.67 0.2 0.04]);

hclosega=uicontrol(hTGAEstimates, 'Style', 'pushbutton', ...
   'String', 'Close', ...
   'Callback', @closeE_cb, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.61 0.2 0.04]);

for i=1:l
   hlinenrInt = uicontrol(hTGAEstimates, 'Style', 'text', ...
      'String', i, ...
      'Units', 'normalized', ...
      'Position', [0.02 0.93-i*0.036 0.04 0.029]);
   hMaterialInt = uicontrol(hTGAEstimates,'Style','text',...
      'String',{handles.E.TGA(i).material},...
      'Units', 'normalized', ...
      'Position',[0.07 0.93-i*0.036 0.15 0.029]);
   hRateInt = uicontrol(hTGAEstimates,'Style','text',...
      'String',{handles.E.TGA(i).rate},...
      'Units', 'normalized', ...
      'Position',[0.24 0.93-i*0.036 0.08 0.029]);
   hGasInt = uicontrol(hTGAEstimates,'Style','text',...
      'String',{handles.E.TGA(i).gas},...
      'Units', 'normalized', ...
      'Position',[0.32 0.93-i*0.036 0.08 0.029]);
   hcheck = uicontrol(hTGAEstimates,'Style','checkbox',...
      'CallBack', {@Echeck_cb,i},...
      'Units', 'normalized', ...
      'Value',0,...
      'Position',[0.46 0.93-i*0.036 0.05 0.04]);
  handles.E.checkbox(i)=hcheck;
end

% Save handles
guidata(handles.TGA.hTGAEstimates,handles);

% Make the GUI visible.
set(hTGAEstimates, 'Visible','on');
uiwait(hTGAEstimates);
% update handles-data
handles = guidata(handles.TGA.hTGAEstimates);
guidata(handles.hPyroPlot, handles);
close(handles.TGA.hTGAEstimates);
end

function alfa_cb(hObject, eventdata)
str=get(get(hObject, 'SelectedObject'), 'String');

switch str
    case 'a=m/m(0)'
        handles.E.a = 1;
    case 'a=(m-(inf))/m(0)'
        handles.E.a = 2;
end
end

function Echeck_cb(hObject, eventdata,nItem)
handles = guidata(hObject);
handles.E.TGA(nItem).check = get(hObject,'Value');
%it is possible to choose only one data series by time
if get(hObject, 'Value')==1
    handles.E.line = nItem;
    for i=1:nItem-1
    set(handles.E.checkbox(i), 'Value', 0);
    handles.E.TGA(i).check = 0;
    end
    for i=nItem+1:length(handles.E.TGA)
    set(handles.E.checkbox(i), 'Value', 0);
    handles.E.TGA(i).check = 0;
    end
    hrorder = uicontrol(handles.TGA.hTGAEstimates, 'Style', 'text', ...
      'String', 'N', ...
      'Units', 'normalized', ...
      'Position', [0.52 0.93-nItem*0.036 0.04 0.029]);
    hordern = uicontrol(handles.TGA.hTGAEstimates, 'Style', 'edit', ...
      'String', handles.E.TGA(nItem).n, ...
      'BackgroundColor', 'w',...
      'Callback', {@rorder_cb,nItem}, ...
      'Units', 'normalized', ...
      'Position', [0.57 0.935-nItem*0.036 0.05 0.029]);
end
guidata(hObject,handles);
end

function rorder_cb(hObject, eventdata, nItem)
handles = guidata(hObject);
handles.E.TGA(nItem).n = str2double(get(hObject, 'String'));
guidata(hObject, handles);
end

%----------------------------------------------------

function findE_cb(hObject, eventdata)
handles = guidata(hObject);
N = handles.E.line;
if N~=0
handles.E.hReaction = figure(...
   'Visible', 'off',...
   'Name', 'Reaction',...
   'toolbar', 'none', ...
   'menubar', 'none', ...
   'NumberTitle', 'off', ...
   'HandleVisibility','callback', ...
   'Units', 'normalized', ...
   'Position',[0.2,0.05,0.6,0.6],...
   'Color', get(0,...
   'defaultuicontrolbackgroundcolor'));
hokR = uicontrol(handles.E.hReaction, 'Style', 'pushbutton', ...
   'string', 'OK', ...
   'Units', 'normalized', ...
   'Position', [0.93 0.5 0.044 0.05], ...
   'Callback', @okReaction_cb);
hclearR = uicontrol(handles.E.hReaction, 'Style', 'pushbutton', ...
   'String', 'Clear', ...
   'Units', 'normalized', ...
   'Position', [0.93 0.4 0.044 0.05], ...
   'Callback', @clearReaction_cb);
set(handles.E.hReaction, 'Visible', 'on');

plot(handles.E.TGA(N).temperature, handles.E.TGA(N).TGA);
title('Choose 2 points to define reaction area');
xlabel('Temperature [\circC]');
ylabel('Mass fraction');
handles = drawArea(handles);
guidata(handles.TGA.hTGAEstimates,handles)
guidata(hObject,handles);
uiwait(handles.E.hReaction);
close(handles.E.hReaction);
else
    msgbox('You must choose one TGA data.');
end

end

function handles = drawArea(handles);
N = handles.E.line;
set(0, 'CurrentFigure', handles.E.hReaction);
currentpointer=get(handles.E.hReaction, 'Pointer');
set(handles.E.hReaction, 'Pointer','fullcrosshair');
pos = ginput(1);
pos(2) = interp1(handles.E.TGA(N).temperature, handles.E.TGA(N).TGA,pos(1));
hold on
handles.E.homark(1) = plot(pos(1), pos(2), 'o');
handles.E.limits(1,:) = [pos(1) pos(2)];

pos = ginput(1);
pos(2) = interp1(handles.E.TGA(N).temperature, handles.E.TGA(N).TGA,pos(1));
set(0, 'CurrentFigure', handles.E.hReaction);
hold on
handles.E.homark(2) = plot(pos(1), pos(2), 'o');
handles.E.limits(2,:) = [pos(1) pos(2)];

if handles.E.limits(1,1) < handles.E.limits(2,1)
   %if first chosen value is smaller than the second
   handles.E.min = handles.E.limits(1,1);
   handles.E.max = handles.E.limits(2,1);
else
   %if second chosen value is smaller than the first
   handles.E.min = handles.E.limits(2,1);
   handles.E.max = handles.E.limits(1,1);
end

handles.E.lb = find(handles.E.TGA(N).temperature<=handles.E.min,1,'last');
handles.E.ub = find(handles.E.TGA(N).temperature>=handles.E.max,1,'first');

%plot box around area
a=zeros(1,handles.E.ub-handles.E.lb+1);
a(1,:)=handles.E.TGA(N).temperature(handles.E.lb);
b=zeros(1,handles.E.ub-handles.E.lb+1);
b(1,:)=handles.E.TGA(N).temperature(handles.E.ub);
c=zeros(1,handles.E.ub-handles.E.lb+1);
c(1,:)=handles.E.TGA(N).TGA(handles.E.lb);
d=zeros(1,handles.E.ub-handles.E.lb+1);
d(1,:)=handles.E.TGA(N).TGA(handles.E.ub);
handles.E.box(1)=plot(a, handles.E.TGA(N).TGA(handles.E.lb:handles.E.ub),'r');
handles.E.box(2)=plot(b, handles.E.TGA(N).TGA(handles.E.lb:handles.E.ub),'r');
handles.E.box(3)=plot(handles.E.TGA(N).temperature(handles.E.lb:handles.E.ub), c,'r');
handles.E.box(4)=plot(handles.E.TGA(N).temperature(handles.E.lb:handles.E.ub), d,'r');
guidata(handles.E.hReaction,handles);

end

function clearReaction_cb(hObject, eventdata)
handles = guidata(hObject);
if isfield(handles.E,'homark')
   delete(handles.E.homark);
   handles.E=rmfield(handles.E, 'homark');
end
if isfield(handles.E, 'box')
   for i=1:4
   delete(handles.E.box(i));
   end
   handles.E = rmfield(handles.E, 'box');
end
handles = drawArea(handles);
guidata(hObject,handles);
end
%----------------------------------
function okReaction_cb(hObject, eventdata)
handles = guidata(hObject);
if isempty(handles)
    msgbox('You must choose two points.');
    return;
end
uiresume(handles.E.hReaction);
N = handles.E.line;
handles.E.hFitLine = figure(...
   'Visible', 'off',...
   'Name', 'Fit Line',...
   'toolbar', 'none', ...
   'menubar', 'none', ...
   'NumberTitle', 'off', ...
   'HandleVisibility','callback', ...
   'Units', 'normalized', ...
   'Position',[0.2,0.05,0.6,0.6],...
   'Color', get(0,...
   'defaultuicontrolbackgroundcolor'));
hFit = uicontrol(handles.E.hFitLine, 'Style', 'pushbutton', ...
   'string', 'Fit', ...
   'Units', 'normalized', ...
   'Position', [0.93 0.5 0.044 0.05], ...
   'Callback', @fitLine_cb);
hOKFit = uicontrol(handles.E.hFitLine, 'Style', 'pushbutton', ...
   'String', 'OK', ...
   'Units', 'normalized', ...
   'Position', [0.93 0.4 0.044 0.05], ...
   'Callback', @okFit_cb);

T = handles.E.TGA(N).temperature(handles.E.lb:handles.E.ub)+273.15;
handles.E.T = T;
alfa = [];
if handles.E.a==1
alfa = handles.E.TGA(N).TGA(handles.E.lb:handles.E.ub)./handles.E.TGA(N).TGA(handles.E.lb);
elseif handles.E.a ==2
alfa = (handles.E.TGA(N).TGA(handles.E.lb:handles.E.ub)-handles.E.TGA(N).TGA(handles.E.ub))./handles.E.TGA(N).TGA(handles.E.lb);
end
handles.E.alfa = alfa;
handles.E.F_a = log(-gradient(alfa,(handles.E.TGA(N).time(2)-handles.E.TGA(N).time(1))))-handles.E.TGA(N).n.*log(alfa);
warning off all
plot(1./T, handles.E.F_a);
xlabel('1/T [1/K]');
ylabel('ln(-da/dt)-Nln(a)');
title('Choose linear area');
warning on
set(handles.E.hFitLine, 'Visible', 'on');

guidata(hObject,handles);

currentpointer=get(handles.E.hFitLine, 'Pointer');
set(handles.E.hFitLine, 'Pointer','fullcrosshair');
pos2 = ginput(1);
pos2(2) = interp1(1./T, handles.E.F_a,pos2(1));
hold on
handles.E.homark2(1) = plot(pos2(1), pos2(2), 'o');
handles.E.limits2(1,:) = [pos2(1) pos2(2)];

pos2 = ginput(1);
pos2(2) = interp1(1./T, handles.E.F_a,pos2(1));
set(0, 'CurrentFigure', handles.E.hFitLine);
hold on
handles.E.homark2(2) = plot(pos2(1), pos2(2), 'o');
handles.E.limits2(2,:) = [pos2(1) pos2(2)];

if handles.E.limits2(1,1) > handles.E.limits2(2,1)
   %if first chosen value is in axis before the second
   handles.E.min2 = handles.E.limits2(1,1);
   handles.E.max2 = handles.E.limits2(2,1);
else
   %if second chosen value is in axis before the second
   handles.E.min2 = handles.E.limits2(2,1);
   handles.E.max2 = handles.E.limits2(1,1);
end

handles.E.lb2 = find(1./T<=handles.E.min2,1,'first');
handles.E.ub2 = find(1./T>=handles.E.max2,1,'last');

%plot box around area
a=zeros(1,handles.E.ub2-handles.E.lb2+1);
a(1,:)=1./T(handles.E.lb2);
b=zeros(1,handles.E.ub2-handles.E.lb2+1);
b(1,:)=1./T(handles.E.ub2);
c=zeros(1,handles.E.ub2-handles.E.lb2+1);
c(1,:)=handles.E.F_a(handles.E.lb2);
d=zeros(1,handles.E.ub2-handles.E.lb2+1);
d(1,:)=handles.E.F_a(handles.E.ub2);
handles.E.box(1)=plot(a, handles.E.F_a(handles.E.lb2:handles.E.ub2),'r');
handles.E.box(2)=plot(b, handles.E.F_a(handles.E.lb2:handles.E.ub2),'r');
handles.E.box(3)=plot(1./T(handles.E.lb2:handles.E.ub2), c,'r');
handles.E.box(4)=plot(1./T(handles.E.lb2:handles.E.ub2), d,'r');

guidata(hObject, handles);
guidata(handles.E.hFitLine,handles);

uiwait(handles.E.hFitLine);
close(handles.E.hFitLine);

end

function fitLine_cb(hObject, eventdata)
handles = guidata(hObject);

if isempty(handles)
    msgbox('You must choose two points.');
    return;
end
if isfield(handles.E, 'fittedLine')
    delete(handles.E.fittedLine);
end

[handles.E.estimates, model] = fitline(handles.E.T(handles.E.lb2:handles.E.ub2), handles.E.F_a(handles.E.lb2:handles.E.ub2));
[sse, FittedCurve] = model(handles.E.estimates);
set(0, 'CurrentFigure', handles.E.hFitLine);
warning off all
handles.E.fittedLine = plot(1./handles.E.T(handles.E.lb2:handles.E.ub2), FittedCurve,'g');
s = ['A = ' num2str(handles.E.estimates(1)) ', E = ' num2str(handles.E.estimates(2))];
warning on
if isfield(handles.E, 'textbox')
    delete(handles.E.textbox);
end
handles.E.textbox = annotation(handles.E.hFitLine,'textbox','String', s,...
    'Position',[0.6 0.8 0.1607 0.0746],...
    'FitHeightToText',...
    'on');
guidata(hObject,handles);
end

function okFit_cb(hObject, eventdata)
handles = guidata(hObject);
if isempty(handles) || isfield(handles.E, 'fittedLine')==0
    msgbox('You must fit line first');
    return;
end
handles.E.TGA(handles.E.line).estimates=handles.E.estimates;
handles.E = rmfield(handles.E, 'fittedLine');
handles.E = rmfield(handles.E, 'textbox');
handles.E = rmfield(handles.E, 'limits');
handles.E = rmfield(handles.E, 'limits2');
handles.E = rmfield(handles.E, 'min');
handles.E = rmfield(handles.E, 'min2');
handles.E = rmfield(handles.E, 'max');
handles.E = rmfield(handles.E, 'max2');
handles.E = rmfield(handles.E, 'lb');
handles.E = rmfield(handles.E, 'lb2');
handles.E = rmfield(handles.E, 'ub');
handles.E = rmfield(handles.E, 'ub2');
guidata(hObject, handles);
guidata(handles.TGA.hTGAEstimates, handles);
uiresume(handles.E.hFitLine);
end

function closeE_cb(hObject, eventdata)
handles = guidata(hObject);
uiresume(handles.TGA.hTGAEstimates);
end
%-----------------------------------------
function showE_cb(hObject, eventdata)
handles = guidata(hObject);
name = [cd, '\Estimates.txt'];
fid = fopen(name, 'w');
fprintf(fid, '%s\n\n', 'ESTIMATES');
fprintf(fid, '%s\t\t\t\t', 'NAME');
fprintf(fid, '%10s\t\t', 'A');
fprintf(fid, '%10s\n', 'E');
for i=1:length(handles.E.TGA)
    fprintf(fid, '%20s\t', handles.E.TGA(i).name);
    if handles.E.TGA(i).estimates(1)~= 0
    fprintf(fid, '%10.2f\t', handles.E.TGA(i).estimates(1));
    fprintf(fid, '%10.2f\n', handles.E.TGA(i).estimates(2));
    else
        fprintf(fid, '%10s\n', 'No estimates');
    end
end
fclose(fid);
open(name);
guidata(hObject, handles);
end


%-----------------------------------------

function [estimates, model] = fitline(xdata, ydata)
start_point = [5+rand*10 (1+rand)];
start_point(1)=10^start_point(1);
start_point(2)=1E5*start_point(2);
%start_point = [1e10 1e5];
model = @lin;
options = optimset('MaxFunEvals',5000);
options = optimset('MaxIter',5000);
estimates = fminsearch(model, start_point, options);

function [sse, FittedCurve] = lin(params)
        A = params(1);
        E = params(2);
        FittedCurve = log(A)-E./(8.3145.*xdata);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end