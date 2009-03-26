function handles = ppdscint(handles)
% PPDSCINT       DSC integration tool for PyroPlot
%
% handles = ppdscint(handles) 

% Initialize
handles.DSC.deltaHvector=[];
handles.DSC.homark = -1;
handles.DSC.hlimitline = -1;
handles.DSC.checkbox = [];
%create new object of only DSC datatype
l=0; %index of structure DSC
lines = handles.var.lines;
for i=1:lines
   if strcmp(handles.EXPDATA(i).type,'DSC')
      l=l+1;
      handles.DSC.DSC(l).rate = handles.EXPDATA(i).rate;
      handles.DSC.DSC(l).temperature = handles.EXPDATA(i).temperature;
      handles.DSC.DSC(l).time = handles.EXPDATA(i).time;
      handles.DSC.DSC(l).DSC = handles.EXPDATA(i).DSC;
      handles.DSC.DSC(l).material = handles.EXPDATA(i).material;
      handles.DSC.DSC(l).gas = handles.EXPDATA(i).gas;
      handles.DSC.DSC(l).name = handles.EXPDATA(i).name;
      handles.DSC.DSC(l).pair = handles.EXPDATA(i).pair; %this shows the index to TGA in structure EXPDATA
      handles.DSC.DSC(l).check = 0;
   end
end

handles.DSC.hDSCInt = figure('Visible','off',...
   'Name', 'DSC-integration', ...
   'NumberTitle', 'off', ...
   'Menubar', 'none', ...
   'Toolbar', 'none', ...
   'Color', get(0,...
   'defaultuicontrolbackgroundcolor'), ...
   'Units', 'normalized', ...
   'Position',[0.36,0.3,0.4,0.7]);
hDSCInt = handles.DSC.hDSCInt;
movegui(hDSCInt,'center')

htitleInt1 = uicontrol(hDSCInt,'Style','text',...
   'String','MATERIAL',...
   'Units', 'normalized', ...
   'Position',[0.07 0.93 0.15 0.03]);
htitleInt2 = uicontrol(hDSCInt,'Style','text',...
   'String','RATE',...
   'Units', 'normalized', ...
   'Position',[0.24 0.93 0.08 0.03]);
htitleInt3 = uicontrol(hDSCInt,'Style','text',...
   'String','GAS',...
   'Units', 'normalized', ...
   'Position',[0.32 0.93 0.08 0.03]);
htitleInt4 = uicontrol(hDSCInt, 'Style', 'text', ...
   'string', 'Select', ...
   'Units', 'normalized', ...
   'Position', [0.40 0.93 0.14 0.03]);
htitleInt5  = uicontrol(hDSCInt, 'Style', 'text', ...
   'String', 'Delta H (kJ/kg)', ...
   'Units', 'normalized', ...
   'Position', [0.50 0.93 0.2 0.03] );

for i=1:l
   hlinenrInt = uicontrol(hDSCInt, 'Style', 'text', ...
      'String', i, ...
      'Units', 'normalized', ...
      'Position', [0.02 0.93-i*0.036 0.04 0.029]);
   hMaterialInt = uicontrol(hDSCInt,'Style','text',...
      'String',{handles.DSC.DSC(i).material},...
      'Units', 'normalized', ...
      'Position',[0.07 0.93-i*0.036 0.15 0.029]);
   hRateInt = uicontrol(hDSCInt,'Style','text',...
      'String',{handles.DSC.DSC(i).rate},...
      'Units', 'normalized', ...
      'Position',[0.24 0.93-i*0.036 0.08 0.029]);
   hGasInt = uicontrol(hDSCInt,'Style','text',...
      'String',{handles.DSC.DSC(i).gas},...
      'Units', 'normalized', ...
      'Position',[0.32 0.93-i*0.036 0.08 0.029]);
   hcheck = uicontrol(hDSCInt,'Style','checkbox',...
      'CallBack', {@DSCcheck_cb,i},...
      'Units', 'normalized', ...
      'Value',0,...
      'Position',[0.46 0.93-i*0.036 0.05 0.04]);
  handles.DSC.checkbox(i) = hcheck;
end

% add popup for reference mass
hMassReftitle = uicontrol(hDSCInt,'Style','text',...
   'String','Reference mass',...
   'Units', 'normalized', ...
   'Position',[0.76, 0.92, 0.2, 0.03]);
handles.DSC.hMRefpopup = uicontrol(hDSCInt,'Style','popupmenu',...
      'String',{'Original','Just before','Consumed'},...
      'BackgroundColor', 'w',...
      'Units', 'normalized', ...
      'Value',3,...
      'Position',[0.76, 0.89, 0.2, 0.03]);

% add pushbuttons 'Integrate' and 'close'
% and text title 'Average'

hIntegrate=uicontrol(hDSCInt, 'Style', 'pushbutton', ...
   'String', 'Integrate', ...
   'Callback', @plotint_cb, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.73 0.1 0.04]);
hclearInt = uicontrol(hDSCInt, 'Style', 'pushbutton', ...
   'String', 'Clear', ...
   'Callback', @clear_cb, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.67 0.1 0.04]);
hcloseint=uicontrol(hDSCInt, 'Style', 'pushbutton', ...
   'String', 'Close', ...
   'Callback', @closeInt_cb, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.61 0.1 0.04]);
haverage = uicontrol(hDSCInt, 'Style', 'text',...
   'string', 'AVERAGE',...
   'Units', 'normalized',...
   'Position', [0.50 0.92-(l+1)*0.036 0.2 0.04]);
% Save handles
guidata(handles.DSC.hDSCInt,handles);
% Make the GUI visible.
set(hDSCInt, 'Visible','on');
% Wait for close-button
uiwait(hDSCInt)
% update handles-data
handles = guidata(handles.DSC.hDSCInt);
close(handles.DSC.hDSCInt);
if (isfield(handles.DSC,'DSC'))
   handles.DSC = rmfield(handles.DSC,'DSC');
end
end

%----------------------------------
% Select DSC data for integration
function DSCcheck_cb(hObject, eventdata,nItem)
handles = guidata(hObject);
handles.DSC.DSC(nItem).check = get(hObject,'Value');
guidata(hObject,handles);
end

% callback of 'Integrate' pushbutton
% plots one selected DSC data and asks user to set integration limits
% after that user pushs 'ok' and programm plots next data from the list
% after every plot programm calculates integral and average of
% non-zero integrals
function plotint_cb(hObject, eventdata)
handles = guidata(hObject);
nDSC = length(handles.DSC.DSC);
%
if ishandle(handles.DSC.deltaHvector)
   for i=1:length(handles.DSC.deltaHvector)
       if handles.DSC.deltaHvector(i)~=0
            delete(handles.DSC.deltaHvector(i))
       end
   end
end
%
handles.DSC.hDSCPlot = figure(...
   'Visible', 'off',...
   'Name', 'DSC',...
   'toolbar', 'none', ...
   'menubar', 'none', ...
   'NumberTitle', 'off', ...
   'HandleVisibility','callback', ...
   'Units', 'normalized', ...
   'Position',[0.2,0.05,0.6,0.6],...
   'Color', get(0,...
   'defaultuicontrolbackgroundcolor'));
hDSCPlot = handles.DSC.hDSCPlot;
hnext = uicontrol(hDSCPlot, 'Style', 'pushbutton', ...
   'string', 'Next', ...
   'Units', 'normalized', ...
   'Position', [0.93 0.5 0.044 0.05], ...
   'Callback', @next_cb);
hclearpoints = uicontrol(hDSCPlot, 'Style', 'pushbutton', ...
   'String', 'Clear', ...
   'Units', 'normalized', ...
   'Position', [0.93 0.4 0.044 0.05], ...
   'Callback', @clearpoints_cb);
set(hDSCPlot, 'Visible', 'on');
%
% Loop over the DSC data
handles.DSC.IntAx = -1;
for i = 1:nDSC
   if (handles.DSC.DSC(i).check)
      pair = handles.DSC.DSC(i).pair;
      handles.DSC.DSCIntIndx=i;
      if ishandle(handles.DSC.IntAx)
         delete(handles.DSC.IntAx)
      end
      if (pair == 0)
         h1 = plot(handles.DSC.DSC(i).temperature, handles.DSC.DSC(i).DSC);
         legend(handles.DSC.DSC(i).name);
         handles.DSC.IntAx = gca;
         %if DSC data doesnt have a pair, programm asks the mass loss from user
         RefMass = get(handles.DSC.hMRefpopup,'value');
         switch RefMass
            case 1
            case 2
               prompt = 'Enter mass fraction prior to reaction [0,1]';
            case 3
               prompt = 'Enter mass fraction consumed [0,1]';
         end
         name = 'Mass loss specification';
         answer = inputdlg(prompt,name);
         handles.DSC.DSC(i).massloss = str2double(answer);
      else
         [ax, h1, h2]=plotyy(handles.EXPDATA(pair).temperature,...
            handles.EXPDATA(pair).TGA*100, ...
            handles.DSC.DSC(i).temperature, handles.DSC.DSC(i).DSC);
         legend([h1 h2], [handles.EXPDATA(pair).name; handles.DSC.DSC(i).name]);
         handles.DSC.IntAx = ax;
      end
      handles = drawIntLine(handles);
      guidata(hDSCPlot,handles);
      uiwait(hDSCPlot);
      handles = guidata(hDSCPlot);
   else
      handles.DSC.DSC(i).intvalue = 0;
   end
end
close(hDSCPlot);
%
% Calculate averages and print individual values
%
nInt = 0;
deltaHave = 0;
for i = 1:nDSC
   if (handles.DSC.DSC(i).check)
      % update average
      nInt = nInt + 1;
      deltaHave = deltaHave + handles.DSC.DSC(i).intvalue;
      dHstr = sprintf('%0.5g',handles.DSC.DSC(i).intvalue/1000); %kJ/kg
      hdeltaH = uicontrol(handles.DSC.hDSCInt, 'Style', 'text', ...
         'String', dHstr, ...
         'Units', 'normalized', ...
         'Position', [0.5 0.93-(i)*0.036 0.2 0.029]);
      handles.DSC.deltaHvector(i)=hdeltaH;
   end
end
if (nInt>0)
   deltaHave = deltaHave / nInt;
end   
handles.DSC.deltaHave = deltaHave;
%
dHstr = sprintf('%0.5g',deltaHave/1000); %kJ/kg
haverage = uicontrol(handles.DSC.hDSCInt, 'Style', 'text', ...
   'String', dHstr,...
   'Units', 'normalized', ...
   'Position', [0.5 0.93-(nDSC+2)*0.036 0.2 0.04]);
%
% create Log
%
handles = createLogfile(handles); 
handles.DSC.analysis = handles.DSC.analysis +1;
%
guidata(hObject,handles);
end

%--------------------------------------------------
% clear callback
%
function clear_cb(hObject, eventdata)
%set DSC(i).check = 0
%set delta H values to zero
%set average value to zero
handles = guidata(hObject);
nDSC = length(handles.DSC.DSC);
for i=1:nDSC
   handles.DSC.DSC(i).check = 0;
   set(handles.DSC.checkbox(i), 'Value', 0);
   handles.DSC.DSC(i).intvalue = 0;
   handles.DSC.DSC(i).min = 0;
   handles.DSC.DSC(i).max = 0;
   handles.DSC.DSCIntIndx = 0;
   handles.DSC.DSCline = 0;
end
if ishandle(handles.DSC.deltaHvector)
   for i=1:length(handles.DSC.deltaHvector)
       if handles.DSC.deltaHvector(i)~=0
            delete(handles.DSC.deltaHvector(i))
       end
   end
end
handles.DSC.average = uicontrol(handles.DSC.hDSCInt, 'Style', 'text', ...
   'String', mat2str(0), ...
   'Units', 'normalized', ...
   'Position', [0.5 0.93-(nDSC+2)*0.036 0.2 0.04]);
guidata(hObject,handles);
end

%--------------------------------------------------
% Callback of 'Close' pushbutton
%
function closeInt_cb(hObject, eventdata)
handles = guidata(hObject);
uiresume(handles.DSC.hDSCInt);
end


% clearpoints_cb
function clearpoints_cb(hObject, eventdata)
handles = guidata(hObject);
if ishandle(handles.DSC.homark)
   delete(handles.DSC.homark);
end
if ishandle(handles.DSC.hlimitline)
   delete(handles.DSC.hlimitline);
end
handles = drawIntLine(handles);
guidata(hObject,handles);
end


% DRAWINTLINE
function handles = drawIntLine(handles)
%
i = handles.DSC.DSCIntIndx;
ax = handles.DSC.IntAx(1);
hDSCPlot = handles.DSC.hDSCPlot;
set(0, 'CurrentFigure', hDSCPlot);
title(ax,'Choose two points of the line');
currentpointer=get(hDSCPlot, 'Pointer');
set(hDSCPlot, 'Pointer','fullcrosshair');
%
pos = ginput(1);
pos(2) = interp1(handles.DSC.DSC(i).temperature, handles.DSC.DSC(i).DSC,pos(1));
hold on
handles.DSC.homark(1) = plot(pos(1), pos(2), 'o');
handles.DSC.limits(1,:) = [pos(1) pos(2)];
%
pos = ginput(1);
pos(2) = interp1(handles.DSC.DSC(i).temperature, handles.DSC.DSC(i).DSC,pos(1));
set(0, 'CurrentFigure', hDSCPlot);
hold on
handles.DSC.homark(2) = plot(pos(1), pos(2), 'o');
handles.DSC.limits(2,:) = [pos(1) pos(2)];
%
handles.DSC.hlimitline = plot(handles.DSC.limits(:,1),handles.DSC.limits(:,2));
set(hDSCPlot, 'Pointer',currentpointer);
hold off
set(0, 'CurrentFigure', hDSCPlot);
title(ax,'Press Next to continue or Clear to choose new points.');
end

%--------------------------------------------------
% Callback of 'next' pushbutton
% The actual integration
%
function next_cb(hObject, eventdata)
handles=guidata(hObject);
i = handles.DSC.DSCIntIndx;
limits = handles.DSC.limits;
[p,v] = size(limits);
if p<2
   herror = msgbox('Error','Choose two points!','error');
   return;
end
handles.DSC.DSC(i).intvalue = 0;
pair = handles.DSC.DSC(i).pair;
if limits(1,1) < limits(2,1)
   %if first chosen value is smaller than the second
   handles.DSC.DSC(i).min = limits(1,1);
   handles.DSC.DSC(i).max = limits(2,1);
else
   %if second chosen value is smaller than the first
   handles.DSC.DSC(i).min = limits(2,1);
   handles.DSC.DSC(i).max = limits(1,1);
end
if (handles.DSC.DSC(i).min ==0 || handles.DSC.DSC(i).max ==0)
   herror = msgbox('Error','Bad integration points.','error');
   return;
end
% Find integration limits
j = find(handles.DSC.DSC(i).temperature<=handles.DSC.DSC(i).min,1,'last');
k = find(handles.DSC.DSC(i).temperature>=handles.DSC.DSC(i).max,1,'first');
Dtime = handles.DSC.DSC(i).time(k)-handles.DSC.DSC(i).time(j);
if (Dtime ==0)
   herror = msgbox('Error','Bad integration points.','error');
   return;
end
% Find straight line limiting the integration area
A = (handles.DSC.DSC(i).DSC(k)-handles.DSC.DSC(i).DSC(j))/Dtime;
B = handles.DSC.DSC(i).DSC(j) - A * handles.DSC.DSC(i).time(j);
Sline = B + A * handles.DSC.DSC(i).time(j:k);
% Find massloss from the TGA-pair
if pair>0
   MassMax = interp1(handles.EXPDATA(pair).temperature,handles.EXPDATA(pair).TGA,handles.DSC.DSC(i).min);
   MassMin = interp1(handles.EXPDATA(pair).temperature,handles.EXPDATA(pair).TGA,handles.DSC.DSC(i).max);
   handles.DSC.DSC(i).massloss = (MassMax-MassMin);
end
% Check for zero massloss
if (handles.DSC.DSC(i).massloss == 0)
   herror = msgbox('Error','Mass loss is zero.','error');
   return;
end
% Integrate
RefMass = get(handles.DSC.hMRefpopup,'value');
switch RefMass
   case 1 % Original
      MRef = 1;
   case 2 % Just before
      if pair>0
         MRef = MassMax;
      else
         MRef = handles.DSC.DSC(i).massloss;
      end
   case 3 % Consumed
      MRef = handles.DSC.DSC(i).massloss;
end
handles.DSC.DSC(i).intvalue =trapz(handles.DSC.DSC(i).time(j:k),...
   (handles.DSC.DSC(i).DSC(j:k)-Sline)) / MRef;
guidata(hObject,handles);
uiresume(handles.DSC.hDSCPlot);
end

%--------------------------------------------------
% createLogFile
% Function creates new Log.txt file or adds new analysis after the end
% of file
function handles = createLogfile(handles)
nDSC = length(handles.DSC.DSC);
analysis = handles.DSC.analysis;
fname = [cd, '\DSC_Log.txt'];
if (analysis == 1 || ~exist(fname,'file'))
   %if first analysis, open new file or overwrite old one
   fid = fopen(fname, 'w');
   fprintf(fid, '%s\n', '% DSC-integration Log');
   fprintf(fid, '%s\n', '% If you want to save this file, change the file name');
   fprintf(fid, '%s\n', ' ');
else
   %file allready exists, write after eof
   fid = fopen(fname, 'a');
end
s = ['Analysis #' mat2str(analysis)];
fprintf(fid, '%10s\n\n', s);
fprintf(fid, '%7s\t',  'Sample');
fprintf(fid, '%25s\t', 'Name');
fprintf(fid, '%10s\t', 'Tmin (C)');
fprintf(fid, '%10s\t', 'Tmax (C)');
RefMass = get(handles.DSC.hMRefpopup,'value');
switch RefMass
   case 1
      fprintf(fid, '%10s\n', 'Delta H (kJ/kg original mass)');
   case 2
      fprintf(fid, '%10s\n', 'Delta H (kJ/kg) (mass just before reaction)');
   case 3
      fprintf(fid, '%10s\n', 'Delta H (kJ/kg consumed mass)');
end
j = 0;
for i=1:nDSC
   if (handles.DSC.DSC(i).check)
      j = j + 1;
      fprintf(fid, '%7d\t', j);
      fprintf(fid, '%25s\t', handles.DSC.DSC(i).name);
      fprintf(fid, '%10.2f\t',  handles.DSC.DSC(i).min);
      fprintf(fid, '%10.2f\t',  handles.DSC.DSC(i).max);
      fprintf(fid, '%0.5g\n',   handles.DSC.DSC(i).intvalue/1000.0);
   end
end
fprintf(fid, '%s\t\t', 'AVERAGE (kJ/kg)');
fprintf(fid, '%0.5g\n', handles.DSC.deltaHave/1000.0);
fprintf(fid, '%40s\n\n', '------------------------------------------------');
fclose(fid);
open(fname);
end
