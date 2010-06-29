function handles = ga(handles, type)

l=0; %index of structure TGA
lines = handles.var.lines;
handles.GA.hvariables = -1;
handles.GA.hparameters = -1;
handles.hvar = -1;
handles.hpar = -1;
handles.hNvar2 = -1;
handles.GA.var = [];
handles.GA.par = [];
handles.GA.estimates = [];
handles.GA.varstr = 'N';
handles.GA.Aw = 15;
handles.GA.P = zeros(1,4); %parameters
handles.GA.parInd = [];

%algorithm parameters
handles.GA.NIND = 20;
handles.GA.MAXGEN = 1200;
handles.GA.MIGGEN = 20;
handles.GA.MUTR = 0.15;
handles.GA.SUBPOP = 4;
handles.GA.GGAP = 0.8;
handles.GA.XOVR = 0.7;
handles.GA.INSR = 0.9;
handles.GA.MIGR = 0.2;

if (strcmp(type,'TGA') || strcmp(type,'TGAFDS'))
for i=1:lines
   if (strcmp(handles.EXPDATA(i).type,'TGA') || strcmp(handles.EXPDATA(i).type,'TGAFDS'))
      l=l+1;
      handles.GA.data(l).rate = handles.EXPDATA(i).rate;
      handles.GA.data(l).temperature = handles.EXPDATA(i).temperature;
      handles.GA.data(l).time = handles.EXPDATA(i).time;
      handles.GA.data(l).TGA = handles.EXPDATA(i).TGA;
      handles.GA.data(l).material = handles.EXPDATA(i).material;
      handles.GA.data(l).gas = handles.EXPDATA(i).gas;
      handles.GA.data(l).name = handles.EXPDATA(i).name;
      handles.GA.data(l).check = 0;
   end
end
elseif (strcmp(type, 'Cone') || strcmp(type, 'ConeFDS'))
    for i=1:lines
        if (strcmp(handles.EXPDATA(i).type,'Cone') || strcmp(handles.EXPDATA(i).type,'ConeFDS'))
      l=l+1;
      handles.GA.data(l).rate = handles.EXPDATA(i).rate;
      handles.GA.data(l).time = handles.EXPDATA(i).time;
      handles.GA.data(l).TGA = handles.EXPDATA(i).TGA;
      handles.GA.data(l).material = handles.EXPDATA(i).material;
      handles.GA.data(l).gas = handles.EXPDATA(i).gas;
      handles.GA.data(l).name = handles.EXPDATA(i).name;
      handles.GA.data(l).MLR = handles.EXPDATA(i).ConeMLR;
      handles.GA.data(l).HRR = handles.EXPDATA(i).ConeHRR;
      handles.GA.data(l).A = 100E-4; %surface area of sample
      handles.GA.data(l).check = 0;
        end
    end
end
handles.GA.hTGAga = figure('Visible','off',...
   'Name', 'Genetic Algorithm', ...
   'NumberTitle', 'off', ...
   'Menubar', 'none', ...
   'Toolbar', 'none', ...
   'Color', get(0,...
   'defaultuicontrolbackgroundcolor'), ...
   'Units', 'normalized', ...
   'Position',[0.36,0.3,0.4,0.7]);
hTGAga = handles.GA.hTGAga;
movegui(hTGAga,'center')

htitleInt1 = uicontrol(hTGAga,'Style','text',...
   'String','MATERIAL',...
   'Units', 'normalized', ...
   'Position',[0.07 0.93 0.15 0.03]);
htitleInt2 = uicontrol(hTGAga,'Style','text',...
   'String','RATE',...
   'Units', 'normalized', ...
   'Position',[0.24 0.93 0.08 0.03]);
htitleInt3 = uicontrol(hTGAga,'Style','text',...
   'String','GAS',...
   'Units', 'normalized', ...
   'Position',[0.32 0.93 0.08 0.03]);
htitleInt4 = uicontrol(hTGAga, 'Style', 'text', ...
   'string', 'Select', ...
   'Units', 'normalized', ...
   'Position', [0.40 0.93 0.14 0.03]);

if strcmp(type,'Cone')
    htitleInt5 = uicontrol(hTGAga, 'Style', 'text', ...
   'string', 'Surface area (cm^2)', ...
   'Units', 'normalized', ...
   'Position', [0.51 0.93 0.14 0.05]);
end
handles.hvar = uicontrol(hTGAga, 'Style', 'pushbutton', ...
   'String', 'Variables', ...
   'Callback', {@variables_cb,type}, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.86 0.15 0.04]);

%if strcmp(type,'TGA')
handles.hpar = uicontrol(hTGAga, 'Style', 'pushbutton', ...
   'String', 'Parameters', ...
   'Callback', {@parameters_cb,type}, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.8 0.15 0.04]);
%end

%if strcmp(type,'TGA')
%handles.hpar = uicontrol(hTGAga, 'Style', 'pushbutton', ...
%   'String', 'Initial values', ...
%   'Callback', @initialEstimates_cb, ...
%   'Units', 'normalized', ...
%   'Position', [0.76 0.73 0.15 0.04]);
%end

if strcmp(type,'Cone')
    hChoose = uicontrol(hTGAga, 'Style', 'text', ...
        'String', 'Choose ', ...
        'Units', 'normalized', ...
        'Position', [0.74 0.73 0.1 0.04]);
    handles.hDataType = uicontrol(hTGAga, 'Style', 'popupmenu', ...
        'String',{'BR','SD'},...
        'BackgroundColor', 'w',...
        'Units', 'normalized', ...
        'Value',1,...
        'Position',[0.84 0.735 0.1 0.04]);
    hHelp =uicontrol(hTGAga, 'Style', 'pushbutton', ...
        'String','?',...
        'Callback', @help_cone_GA, ...
        'Units', 'normalized', ...
        'Position',[0.95 0.745 0.04 0.03]);
end

hAlgorithm = uicontrol(hTGAga, 'Style', 'pushbutton', ...
   'String', 'Algorithm', ...
   'Callback', {@GA_param_cb,type}, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.67 0.1 0.04]);

hrunGA=uicontrol(hTGAga, 'Style', 'pushbutton', ...
   'String', 'Run GA', ...
   'Callback', {@runGA, type}, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.61 0.1 0.04]);

hstopGA=uicontrol(hTGAga, 'Style', 'pushbutton', ...
   'String', 'STOP', ...
   'Callback', @stopGA_cb, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.55 0.1 0.04]);

hclosega=uicontrol(hTGAga, 'Style', 'pushbutton', ...
   'String', 'Close', ...
   'Callback', @closeGA_cb, ...
   'Units', 'normalized', ...
   'Position', [0.76 0.49 0.1 0.04]);

%hfocus = uicontrol(hTGAga, 'Style', 'pushbutton', ...
%   'String', 'Focus', ...
%   'Callback', @grGA_cb, ...
%   'Units', 'normalized', ...
%   'Position', [0.76 0.49 0.1 0.04]);

for i=1:l
   hlinenrInt = uicontrol(hTGAga, 'Style', 'text', ...
      'String', i, ...
      'Units', 'normalized', ...
      'Position', [0.02 0.93-i*0.036 0.04 0.029]);
   hMaterialInt = uicontrol(hTGAga,'Style','text',...
      'String',{handles.GA.data(i).material},...
      'Units', 'normalized', ...
      'Position',[0.07 0.93-i*0.036 0.15 0.029]);
   hRateInt = uicontrol(hTGAga,'Style','text',...
      'String',{handles.GA.data(i).rate},...
      'Units', 'normalized', ...
      'Position',[0.24 0.93-i*0.036 0.08 0.029]);
   hGasInt = uicontrol(hTGAga,'Style','text',...
      'String',{handles.GA.data(i).gas},...
      'Units', 'normalized', ...
      'Position',[0.32 0.93-i*0.036 0.08 0.029]);
   hcheck = uicontrol(hTGAga,'Style','checkbox',...
      'CallBack', {@GAcheck_cb,i},...
      'Units', 'normalized', ...
      'Value',0,...
      'Position',[0.46 0.93-i*0.036 0.05 0.04]);
  if strcmp(type,'Cone')
     A = handles.GA.data(i).A*10^4;
     hArea = uicontrol(hTGAga,'Style','edit',...
     'BackgroundColor', 'w',...
     'String', A, ...
     'Units', 'normalized', ...
     'Callback', {@area_cb, i}, ...
     'Position',[0.55 0.93-i*0.036 0.07 0.03]);
  end
end

% Save handles
guidata(handles.GA.hTGAga,handles);

% Make the GUI visible.
set(hTGAga, 'Visible','on');
uiwait(hTGAga);
% update handles-data
handles = guidata(handles.GA.hTGAga);
guidata(handles.hPyroPlot, handles);
close(handles.GA.hTGAga);
handles.GA = rmfield(handles.GA,'hTGAga');
if isfield(handles.GA, 'hvariables')
handles.GA = rmfield(handles.GA, 'hvariables');
end
if isfield(handles.GA, 'hparameters')
handles.TGA = rmfield(handles.GA, 'hparameters');
end
if isfield(handles.TGA, 'hTfin')
handles.GA = rmfield(handles.GA, 'hTfin');
handles.GA = rmfield(handles.GA, 'hTfinE');
handles.GA = rmfield(handles.GA, 'hTini');
handles.GA = rmfield(handles.GA, 'hTiniE');
end

if isfield(handles.GA, 'data')
    handles.GA = rmfield(handles.GA, 'data');
end

if isfield(handles.GA, 'variables')
handles.TGA = rmfield(handles.GA, 'variables');
end
if isfield(handles.GA, 'parameters')
handles.GA = rmfield(handles.GA, 'parameters');
end
if isfield(handles.GA, 'var')
    handles.GA = rmfield(handles.GA, 'var');
end
if isfield(handles.GA, 'par')
    handles.GA = rmfield(handles.GA, 'par');
end
if isfield(handles.GA, 'Tindex')
    handles.GA = rmfield(handles.GA, 'Tindex');
end
if isfield(handles.GA, 'Findex')
    handles.GA = rmfield(handles.GA, 'Findex');
end
if isfield(handles.GA, 'Rindex')
    handles.GA = rmfield(handles.GA, 'Rindex');
end
if isfield(handles.GA, 'LogScaling')
    handles.GA = rmfield(handles.GA, 'LogScaling');
end
if isfield(handles.GA, 'log')
    handles.GA = rmfield(handles.GA, 'log');
end
if isfield(handles.GA, 'limits')
    handles.GA = rmfield(handles.GA, 'limits');
end
if isfield(handles.GA, 'parInd')
    handles.GA = rmfield(handles.GA, 'parInd');
end
if isfield(handles.GA, 'Pindex')
    handles.GA = rmfield(handles.GA, 'Pindex');
end
if isfield(handles.GA, 'Hres')
    handles.GA = rmfield(handles.GA, 'Hres');
end
if isfield(handles.GA, 'checkbox')
    handles.GA = rmfield(handles.GA, 'checkbox');
end
if isfield(handles.GA, 'Aindex')
    handles.GA = rmfield(handles.GA, 'Aindex');
end
if isfield(handles.GA, 'P')
    handles.GA = rmfield(handles.GA, 'P');
end
if isfield(handles.GA, 'estimates')
    handles.GA = rmfield(handles.GA, 'estimates');
end

guidata(handles.hPyroPlot, handles.GA);
end

%----------------------------------

function area_cb(hObject, eventdata,N)
handles = guidata(hObject);
str = get(hObject, 'String');
A = str2double(str);
handles.GA.data(N).A = A/10000;
guidata(hObject,handles);
end

function help_cone_GA(hObject, eventdata)
open('Help_cone_GA.txt');
end

function variables_cb(hObject, eventdata,type)
handles = guidata(hObject);
[FDS_file, handles.template_path]=uigetfile('*.fds', 'Choose FDS Input File', handles.template_path);
if (~isequal(FDS_file,0) && ~isequal(handles.template_path,0))
handles.template = fullfile(handles.template_path, FDS_file);
open(handles.template);
%[Exe_file, Exe_path]=uigetfile('*.exe', 'Choose Executive File');
handles.FdsExe = 'fds5.exe';

if strcmp(type, 'TGA')
    handles.GA.weights = [492/1000,493/1000, 15/1000];
else
    handles.GA.weights = [1 0 0];
end
guidata(hObject,handles);

handles.GA.hvariables = figure('Visible','off',...
      'Name', 'Variables for GA', ...
      'Menubar', 'none', ...
      'NumberTitle', 'off', ...
      'Toolbar', 'none', ...
      'Color', get(0,'defaultuicontrolbackgroundcolor'), ...
      'Units', 'normalized', ...
      'Position',[0.4,0.6,0.4,0.7]);
% Add items for GUI
handles = AddVar(handles,type);
guidata(hObject, handles);
% save handles
   guidata(handles.GA.hvariables,handles)
   % Assign the GUI a name to appear in the window title.
   movegui(handles.GA.hvariables,'center')
   % Make the GUI visible.
   set(handles.GA.hvariables,'Visible','on');
   % Wait for ok-button
   uiwait(handles.GA.hvariables)
   % update handles-data
   handles = guidata(handles.GA.hvariables);
   guidata(handles.GA.hTGAga, handles);
   close(handles.GA.hvariables);
end
end

function stopGA_cb(hObject, eventdata)
fid = fopen('ga.stop', 'w');
fclose(fid);
msgbox('Wait, algorithm stops after this generation');
end

function closeGA_cb(hObject, eventdata)
handles = guidata(hObject);
uiresume(handles.GA.hTGAga);
end

function GAcheck_cb(hObject, eventdata,nItem)
handles = guidata(hObject);
handles.GA.data(nItem).check = get(hObject,'Value');
guidata(hObject,handles);
end

function handles = AddVar(handles,type)
    h = 0.035;
    w = 0.05;
    d = 0.04;
    hNvar = uicontrol(handles.GA.hvariables,'Style','text',...
   'String', 'Number of variables',...
   'Units', 'normalized', ...
   'Position',[0.001 0.965 0.25 0.02]);
    handles.hNvar2 = uicontrol(handles.GA.hvariables,'Style','edit',...
   'String',handles.GA.varstr,...
   'BackgroundColor', 'w',...
   'Units', 'normalized', ...
   'Callback', @var_cb, ...
   'Position',[0.25 0.96 w h]);
    if strcmp(type, 'TGA')
       
        hAw = uicontrol(handles.GA.hvariables,'Style','text',...
           'String', 'Weights (%)',...
           'Units', 'normalized', ...
           'Position',[0.31 0.965 0.2 0.02]);
        
       text1 = uicontrol(handles.GA.hvariables,'Style','text',...
           'String', 'Mass',...
           'Units', 'normalized', ...
           'Position',[0.5 0.965 w+0.02 0.02]);
       hAw2 = uicontrol(handles.GA.hvariables,'Style','edit',...
           'String',num2str(handles.GA.weights(1)*100),...
           'BackgroundColor', 'w',...
           'Units', 'normalized', ...
           'Callback', @w_mass_cb, ...
           'Position',[0.6 0.96 w+0.02 h]);
        text2 = uicontrol(handles.GA.hvariables,'Style','text',...
           'String', 'Grad.',...
           'Units', 'normalized', ...
           'Position',[0.68 0.965 w+0.02 0.02]);
       hAw3 = uicontrol(handles.GA.hvariables,'Style','edit',...
           'String',num2str(handles.GA.weights(2)*100),...
           'BackgroundColor', 'w',...
           'Units', 'normalized', ...
           'Callback', @w_grad_cb, ...
           'Position',[0.75 0.96 w+0.02 h]);
       text1 = uicontrol(handles.GA.hvariables,'Style','text',...
           'String', 'A',...
           'Units', 'normalized', ...
           'Position',[0.82 0.965 w+0.02 0.02]);
       hAw4 = uicontrol(handles.GA.hvariables,'Style','edit',...
           'String',num2str(handles.GA.weights(3)*100),...
           'BackgroundColor', 'w',...
           'Units', 'normalized', ...
           'Callback', @w_A_cb, ...
           'Position',[0.89 0.96 w+0.02 h]);
        
    end
    hisA = uicontrol(handles.GA.hvariables,'Style','text',...
   'String', 'A',...
   'Units', 'normalized', ...
   'Position',[0.1 0.96-2*d 0.1 0.03]);
    hlogSc = uicontrol(handles.GA.hvariables,'Style','text',...
   'String', 'Log-scale',...
   'Units', 'normalized', ...
   'Position',[0.25 0.96-2*d 0.1 0.03]);
    hmin = uicontrol(handles.GA.hvariables,'Style','text',...
   'String', 'Min',...
   'Units', 'normalized', ...
   'Position',[0.4 0.96-2*d 0.1 0.03]);
    hmax= uicontrol(handles.GA.hvariables,'Style','text',...
   'String', 'Max',...
   'Units', 'normalized', ...
   'Position',[0.55 0.96-2*d 0.1 0.03]);
    hID = uicontrol(handles.GA.hvariables,'Style','text',...
   'String', 'ID',...
   'Units', 'normalized', ...
   'Position',[0.72 0.96-2*d 0.1 0.03]);
    str = ['Executive file: ' handles.FdsExe];
    handles.GA.hexe = uicontrol(handles.GA.hvariables,'Style','text',...
   'String', str,...
   'Units', 'normalized', ...
   'Position',[0.01 0.965-d 0.4 0.02]);
    hexechange = uicontrol(handles.GA.hvariables, 'Style', 'pushbutton', ...
   'String', 'Change', ...
   'Callback', @changeExe_cb, ...
   'Units', 'normalized', ...
   'Position', [0.4 0.96-d 0.1 h]);
    hOkvar = uicontrol(handles.GA.hvariables, 'Style', 'pushbutton', ...
   'String', 'OK', ...
   'Callback', @okvar_cb, ...
   'Units', 'normalized', ...
   'Position', [0.35 0.01 0.1 0.03]);
    hcancelVar = uicontrol(handles.GA.hvariables,  'Style', 'pushbutton', ...
   'String', 'Cancel', ...
   'Callback', @cancelvar_cb, ...
   'Units', 'normalized', ...
   'Position', [0.55 0.01 0.1 0.03]);
end

function var_cb(hObject, eventdata)
handles = guidata(hObject);
var = str2double(get(hObject, 'string'));
handles.GA.varstr = get(hObject, 'string');
handles.GA.variables = var;
handles.GA.log = [];
if isfield(handles.GA.var, 'check')==0
    for i=1:var  
    handles.GA.var(i).check = 0;
    handles.GA.var(i).log = 0;
    handles.GA.var(i).min = 0;
    handles.GA.var(i).max = 0;   
    s = ['VAR' num2str(i)];
    handles.GA.var(i).id = s;
    end
    guidata(hObject,handles);
end
    h = 0.03;
    w = 0.05;
    d = 0.04;
for i=1:var 
    s = ['VAR' num2str(i)];
    hvari = uicontrol(handles.GA.hvariables,'Style','text',...
   'String', s,...
   'Units', 'normalized', ...
   'Position',[0.001 0.875-i*d 0.1 h]);
    hAch = uicontrol(handles.GA.hvariables,'Style','checkbox',...
      'CallBack', {@Ach_cb,i},...
      'Units', 'normalized', ...
      'Value',handles.GA.var(i).check,...
      'Position',[0.135 0.88-i*d 0.8 h]);
    hlogch = uicontrol(handles.GA.hvariables,'Style','checkbox',...
      'CallBack', {@log_cb,i},...
      'Units', 'normalized', ...
      'Value',handles.GA.var(i).log,...
      'Position',[0.285 0.88-i*d 0.8 h]);
     hminedit = uicontrol(handles.GA.hvariables,'Style','edit',...
     'BackgroundColor', 'w',...
     'String', num2str(handles.GA.var(i).min), ...
     'Units', 'normalized', ...
     'Callback', {@min_cb, i}, ...
     'Position',[0.4 0.88-i*d 0.1 h]);
     hmaxedit = uicontrol(handles.GA.hvariables,'Style','edit',...
     'BackgroundColor', 'w',...
     'String', num2str(handles.GA.var(i).max), ...
     'Units', 'normalized', ...
     'Callback', {@max_cb,i}, ...
     'Position',[0.55 0.88-i*d 0.1 h]);
    hIDedit = uicontrol(handles.GA.hvariables,'Style','edit',...
     'String', handles.GA.var(i).id, ...
     'BackgroundColor', 'w',...
     'Units', 'normalized', ...
     'Callback', {@id_cb,i}, ...
     'Position',[0.7 0.88-i*d 0.15 h]);
 handles.GA.log(i) = hlogch;
end
guidata(hObject,handles);
end

function w_mass_cb(hObject, eventdata)
handles = guidata(hObject);
w = str2double(get(hObject, 'string'));
handles.GA.weights(1) = w./100;
guidata(hObject,handles);
end

function w_grad_cb(hObject, eventdata)
handles = guidata(hObject);
w = str2double(get(hObject, 'string'));
handles.GA.weights(2) = w./100;
guidata(hObject,handles);
end

function w_A_cb(hObject, eventdata)
handles = guidata(hObject);
w = str2double(get(hObject, 'string'));
handles.GA.weights(3) = w./100;
guidata(hObject,handles);
end

function Ach_cb(hObject, eventdata, N)
handles = guidata(hObject);
handles.GA.var(N).check = get(hObject,'Value');
if get(hObject,'Value') == 1
    handles.GA.var(N).log = 1;
    set(handles.GA.log(N), 'Value', 1);
end
guidata(hObject,handles);
end

function log_cb(hObject, eventdata, N)
handles = guidata(hObject);
handles.GA.var(N).log = get(hObject,'Value');
guidata(hObject,handles);
end

function min_cb(hObject, eventdata,N)
handles = guidata(hObject);
handles.GA.var(N).min = str2double(get(hObject, 'string'));
guidata(hObject,handles);
end

function max_cb(hObject, eventdata,N)
handles = guidata(hObject);
handles.GA.var(N).max = str2double(get(hObject, 'string'));
guidata(hObject,handles);
end

function id_cb(hObject, eventdata,N)
handles = guidata(hObject);
handles.GA.var(N).id = get(hObject, 'string');
guidata(hObject,handles);
end

function changeExe_cb(hObject, evetdata)
handles = guidata(hObject);
[v, cdir] = system('cd');
[Exe_file, handles.FdsExe_path]=uigetfile('*.exe', 'Choose Executive File', handles.FdsExe_path);
handles.FdsExe = fullfile(handles.FdsExe_path,Exe_file);
% k=1;
% for i=1:min(length(handles.FdsExe),length(cdir))
%     if handles.FdsExe(i)==cdir(i)
%         k=k+1;
%     else
%         break;
%     end
% end
% l=0;
% if length(cdir) > k
%     l=1;
%     for i=k:length(cdir)
%         if cdir(i)=='\'
%             l=l+1;
%         end
%     end
% end
% s = [];
% for i=1:l
%     s = [s '..\'];
% end
% if handles.FdsExe(k)=='\'
%     k=k+1;
% end
%handles.FdsExe = handles.FdsExe(k:length(handles.FdsExe));
%handles.FdsExe = ['"' s handles.FdsExe '"'];
set(handles.GA.hexe, 'String', ['Executive file: ' handles.FdsExe]);
guidata(hObject,handles);
end

function okvar_cb(hObject, eventdata)
handles = guidata(hObject);
handles.GA.limits=[];
handles.GA.LogScaling = [];
handles.GA.Aindex = [];
k=1;
for i=1:handles.GA.variables
    if handles.GA.var(i).check == 1
       handles.GA.Aindex(k)=i;
       k=k+1;
    end
    handles.GA.limits(1,i)=handles.GA.var(i).min;
    handles.GA.limits(2,i)=handles.GA.var(i).max;
    handles.GA.LogScaling(i) = handles.GA.var(i).log;
    
end
guidata(hObject,handles);
uiresume(handles.GA.hvariables);
end

function cancelvar_cb(hObject, eventdata)
handles = guidata(hObject);
uiresume(handles.GA.hvariables);
end
%-----------------------------------
function parameters_cb(hObject, eventdata, type)
handles = guidata(hObject);
if isfield(handles.GA, 'variables')
handles.GA.par = [];
guidata(hObject, handles);
handles.GA.hparameters = figure('Visible','off',...
      'Name', 'Parameters for GA', ...
      'Menubar', 'none', ...
      'NumberTitle', 'off', ...
      'Toolbar', 'none', ...
      'Color', get(0,'defaultuicontrolbackgroundcolor'), ...
      'Units', 'normalized', ...
      'Position',[0.4,0.6,0.4,0.5]);
    % Add items for GUI
   
    handles = AddPar(handles,type);
    
    % save handles
   guidata(handles.GA.hparameters,handles);
   % Assign the GUI a name to appear in the window title.
   movegui(handles.GA.hparameters,'center');
   % Make the GUI visible.
   set(handles.GA.hparameters,'Visible','on');
   % Wait for ok-button
   uiwait(handles.GA.hparameters);
   % update handles-data
   handles = guidata(handles.GA.hparameters);
   guidata(handles.GA.hTGAga, handles);
   % update handles-data
   close(handles.GA.hparameters);
else
   msgbox('Define first the variables.'); 
end
end

function handles = AddPar(handles, type)
    hNpar = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 'Number of Parameters',...
   'Units', 'normalized', ...
   'Position',[0.001 0.95 0.25 0.02]);
    hNpar2 = uicontrol(handles.GA.hparameters,'Style','edit',...
   'String','N',...
   'BackgroundColor', 'w',...
   'Units', 'normalized', ...
   'Callback', {@par_cb,type}, ...
   'Position',[0.25 0.93 0.08 0.05]);
    hHelp = uicontrol(handles.GA.hparameters, 'Style', 'pushbutton', ...
   'String', 'HELP!', ...
   'Callback', @helpPar_cb, ...
   'Units', 'normalized', ...
   'Position', [0.5 0.93 0.1 0.04]);
    
if strcmp(type, 'TGA')
    hrampT = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 't(T_MAX)',...
   'Units', 'normalized', ...
   'Position',[0.1 0.85 0.1 0.05]);
    htwfin = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 'TWFIN',...
   'Units', 'normalized', ...
   'Position',[0.2 0.85 0.1 0.05]);
else
    hrampT = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 'HR',...
   'Units', 'normalized', ...
   'Position',[0.1 0.85 0.1 0.05]);
end
    hnufuel = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 'NU_FUEL',...
   'Units', 'normalized', ...
   'Position',[0.3 0.85 0.1 0.05]);
    hramp_value = uicontrol(handles.GA.hparameters,'Style','text',...
       'String', 'Ramp value',...
       'Units', 'normalized', ...
       'Position',[0.4 0.85 0.1 0.05]);
%     hmassfraction = uicontrol(handles.GA.hparameters,'Style','text',...
%    'String', 'Mass fraction',...
%    'Units', 'normalized', ...
%    'Position',[0.4 0.85 0.1 0.05]);

    % Add parameter for cone heater radiation level!
    
   %hCone = uicontrol(handles.GA.hparameters,'Style','text',...
   %'String', 'TWFIN',...
   %'Units', 'normalized', ...
   %'Position',[0.2 0.85 0.1 0.05]);
    
    hOkpar = uicontrol(handles.GA.hparameters, 'Style', 'pushbutton', ...
   'String', 'OK', ...
   'Callback', {@okpar_cb,type}, ...
   'Units', 'normalized', ...
   'Position', [0.35 0.01 0.1 0.04]);
    hcancelPar = uicontrol(handles.GA.hparameters, 'Style', 'pushbutton', ...
   'String', 'Cancel', ...
   'Callback', @cancelpar_cb, ...
   'Units', 'normalized', ...
   'Position', [0.55 0.01 0.1 0.04]);
    guidata(handles.GA.hparameters, handles);
    handles = guidata(handles.GA.hparameters);
end

function par_cb (hObject, eventdata, type)
handles = guidata(hObject);
par = str2double(get(hObject, 'string'));
handles.GA.parameters = par;

handles.GA.Tini = 20;
handles.GA.checkbox = [];
handles.GA.Hres = zeros(par, 2);
handles.GA.moisture = 0;
% handles.GA.RampT = 0;
handles.GA.Tindex = 0; 
handles.GA.Findex = 0;
handles.GA.Rindex = []; %residue
handles.GA.Mindex = []; %mass fraction
handles.GA.HR = 0;

guidata(hObject, handles);
for i=1:par
 
    handles.GA.par(i).T = 0; % 0 if not checked
    handles.GA.par(i).F = 0;
    handles.GA.par(i).R = 0;
    handles.GA.par(i).res = 0; %0 if not nu_fuel, otherwise index of variable
    handles.GA.par(i).ramp = 0; %0 if not ramp_value, otherwise index of variable
    handles.GA.par(i).HR = 0;
%     handles.GA.par(i).M = 0;
%     handles.GA.par(i).massfraction = []; %empty if none, otherwise indexes
    s = ['PAR' num2str(i)];
    hpari = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', s,...
   'Units', 'normalized', ...
   'Position',[0.001 0.85-i*0.06 0.1 0.05]);
    hT = uicontrol(handles.GA.hparameters,'Style','checkbox',...
      'CallBack', {@rampT_cb,i},...
      'Units', 'normalized', ...
      'Value',0,...
      'Position',[0.135 0.86-i*0.06 0.8 0.05]);
    hF = uicontrol(handles.GA.hparameters,'Style','checkbox',...
      'CallBack', {@twfin_cb,i},...
      'Units', 'normalized', ...
      'Value',0,...
      'Position',[0.235 0.86-i*0.06 0.8 0.05]);
    hR = uicontrol(handles.GA.hparameters,'Style','checkbox',...
      'CallBack', {@residue_cb,i},...
      'Units', 'normalized', ...
      'Value',0,...
      'Position',[0.335 0.86-i*0.06 0.8 0.05]);
    hRV = uicontrol(handles.GA.hparameters,'Style','checkbox',...
      'CallBack', {@ramp_value_cb,i},...
      'Units', 'normalized', ...
      'Value',0,...
      'Position',[0.435 0.86-i*0.06 0.8 0.05]);
%     hM = uicontrol(handles.GA.hparameters,'Style','checkbox',...
%       'CallBack', {@massfraction_cb,i},...
%       'Units', 'normalized', ...
%       'Value',0,...
%       'Position',[0.435 0.86-i*0.06 0.8 0.05]);
  
    handles.GA.checkbox(i,:) = [hT hF hR hRV];
    
  
end
guidata(hObject, handles);
%handles = guidata(hObject);

%guidata(hObject, handles);
end

function helpPar_cb(hObject, eventdata)
f = [cd, '\Help_par.txt'];
open(f);
end

function rampT_cb(hObject, eventdata, N)
handles = guidata(hObject);
handles.GA.par(N).T = get(hObject,'Value');
if get(hObject, 'Value')==1
    handles.GA.Tindex = N;
    if handles.GA.par(N).R == 1
        delete(handles.GA.Hres(N,:));
    end
    handles.GA.par(N).F = 0;
    handles.GA.par(N).R = 0;
    handles.HA.par(N).ramp = 0;
    set(handles.GA.checkbox(N,2), 'Value', 0);
    set(handles.GA.checkbox(N,3), 'Value', 0);
    set(handles.GA.checkbox(N,4), 'Value', 0); 
    for i =1:N-1
        handles.GA.par(i).T = 0;
        set(handles.GA.checkbox(i,1), 'Value', 0);
        if isfield(handles.GA, 'hTfin')
            handles.GA
            delete(handles.GA.hTfin);
            delete(handles.GA.hTfinE);
            delete(handles.GA.hTini);
            delete(handles.GA.hTiniE);
            handles.GA= rmfield(handles.GA, 'hTfin');
            handles.GA = rmfield(handles.GA, 'hTfinE');
            handles.GA = rmfield(handles.GA, 'hTini');
            handles.GA = rmfield(handles.GA, 'hTiniE');
           % handles.GA.RampT = 0;
        end
    end
    for i=N+1:handles.GA.parameters
        handles.GA.par(i).T = 0;
        set(handles.GA.checkbox(i,1), 'Value', 0);
        if isfield(handles.GA, 'hTfin')
            handles.GA
            delete(handles.GA.hTfin);
            delete(handles.GA.hTfinE);
            delete(handles.GA.hTini);
            delete(handles.GA.hTiniE);
            handles.GA = rmfield(handles.GA, 'hTfin');
            handles.GA = rmfield(handles.GA, 'hTfinE');
            handles.GA = rmfield(handles.GA, 'hTini');
            handles.GA = rmfield(handles.GA, 'hTiniE');
        end
    end
    handles.GA.hTfin = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 'Max TMP',...
   'Units', 'normalized', ...
   'Position',[0.55 0.85-N*0.06 0.1 0.05]);
    handles.GA.hTfinE = uicontrol(handles.GA.hparameters,'Style','edit',...
   'BackgroundColor', 'w',...
   'Units', 'normalized', ...
   'Position',[0.65 0.86-N*0.06 0.09 0.04]);
    handles.GA.hTini = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 'Initial TMP',...
   'Units', 'normalized', ...
   'Position',[0.75 0.85-N*0.06 0.1 0.05]);
    handles.GA.hTiniE = uicontrol(handles.GA.hparameters,'Style','edit',...
        'String', '20', ...
   'BackgroundColor', 'w',...
   'Units', 'normalized', ... 
   'Position',[0.85 0.86-N*0.06 0.09 0.04]);
   %handles.GA.RampT = N;
   %'Callback', @Tini_cb, ...
   %'Callback', @Tfin_cb, ...
end
%handles = guidata(handles.GA.hparameters);
%guidata(handles.GA.hTGAga, handles);
guidata(hObject,handles);
end

function Tfin_cb(hObject, eventdata)
handles = guidata(hObject);
handles.GA.Tfin = str2double(get(hObject, 'String'));
guidata(handles.GA.hparameters,handles);
end

function Tini_cb(hObject, eventdata)
handles = guidata(hObject);
handles.GA.Tini = str2double(get(hObject, 'String'));
guidata(handles.GA.hparameters,handles);
end

function twfin_cb(hObject, eventdata, N)
handles = guidata(hObject);
handles.GA.par(N).F = get(hObject,'Value');
if get(hObject, 'Value')==1
    handles.GA.Findex = N;
    if handles.GA.par(N).R == 1
        delete(handles.GA.Hres(N,:));
    end
    if handles.GA.par(N).T==1
            delete(handles.GA.hTfin);
            delete(handles.GA.hTfinE);
            delete(handles.GA.hTini);
            delete(handles.GA.hTiniE);
            handles.GA = rmfield(handles.GA, 'hTfin');
            handles.GA = rmfield(handles.GA, 'hTfinE');
            handles.GA = rmfield(handles.GA, 'hTini');
            handles.GA = rmfield(handles.GA, 'hTiniE');
    end
    handles.GA.par(N).T = 0;
    handles.GA.par(N).R = 0;
    handles.GA.par(N).ramp = 0;
    set(handles.GA.checkbox(N,1), 'Value', 0);
    set(handles.GA.checkbox(N,3), 'Value', 0);
    set(handles.GA.checkbox(N,4), 'Value', 0); 
    for i =1:N-1
        handles.GA.par(i).F = 0;
        set(handles.GA.checkbox(i,2), 'Value', 0);        
    end
    for i=N+1:handles.GA.parameters
        handles.GA.par(i).F = 0;
        set(handles.GA.checkbox(i,2), 'Value', 0);
    end
end
guidata(hObject,handles);
end

function residue_cb(hObject, eventdata, N,type)
handles = guidata(hObject);
handles.GA.par(N).R = get(hObject,'Value');

 if get(hObject, 'Value')==1
    
    handles.GA.Rindex(length(handles.GA.Rindex)+1)=N;
    if handles.GA.par(N).T==1
            delete(handles.GA.hTfin);
            delete(handles.GA.hTfinE);
            delete(handles.GA.hTini);
            delete(handles.GA.hTiniE);
            handles.GA = rmfield(handles.GA, 'hTfin');
            handles.GA = rmfield(handles.GA, 'hTfinE');
            handles.GA = rmfield(handles.GA, 'hTini');
            handles.GA = rmfield(handles.GA, 'hTiniE');
    end
    handles.GA.par(N).T = 0;
    handles.GA.par(N).F = 0;
    handles.GA.par(N).ramp = 0;
    set(handles.GA.checkbox(N,1), 'Value', 0);
    set(handles.GA.checkbox(N,2), 'Value', 0);    
    set(handles.GA.checkbox(N,4), 'Value', 0); 
   
    hresF = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 'Corresponding variable: VAR ',...
   'Units', 'normalized', ...
   'Position',[0.48 0.85-N*0.06 0.35 0.05]);
    
if isfield(handles.GA, 'variables')
    s = ['{'];
    for i=1:handles.GA.variables
    s = [s '''' handles.GA.var(i).id ''''];
    if i~=handles.GA.variables
        s = [s ','];
    end
    end
    s =[s '}'];
    hresFE = uicontrol(handles.GA.hparameters,'Style','popupmenu',...
   'String', eval(s),...
   'Units', 'normalized', ...
   'BackgroundColor', 'w',...
   'Callback', {@resIpop_cb,N}, ...
   'Position', [0.8 0.87-N*0.06 0.15 0.04]);
    handles.GA.Hres(N,:) = [hresF, hresFE];
else
    
    hresFE = uicontrol(handles.GA.hparameters,'Style','edit',...
   'BackgroundColor', 'w',...
   'Units', 'normalized', ...
   'Callback', {@resI_cb,N}, ...
   'Position', [0.87 0.87-N*0.06 0.07 0.04]);
    handles.GA.Hres(N,:) = [hresF, hresFE];
end
end
%handles = guidata(handles.GA.hparameters);
guidata(hObject, handles);
%guidata(handles.GA.hTGAga,handles);
end

function resIpop_cb(hObject, eventdata, N)
handles = guidata(hObject);
str = get(hObject, 'String');
val = get(hObject,'Value');
     for i=1:handles.GA.variables
         if strcmp(str{val},handles.GA.var(i).id)
             handles.GA.par(N).res = i;
         end
     end
%guidata(handles.GA.hTGAga,handles);
guidata(hObject, handles);
end

function resI_cb(hObject, eventdata, N)
handles = guidata(hObject);
handles.GA.par(N).res = str2double(get(hObject, 'String'));
guidata(handles.GA.hTGAga,handles);
end

function ramp_value_cb(hObject, eventdata, N)
handles = guidata(hObject);

handles.GA.par(N).ramp = get(hObject, 'Value');
if get(hObject, 'Value')==1
    
    if handles.GA.par(N).T==1
            delete(handles.GA.hTfin);
            delete(handles.GA.hTfinE);
            delete(handles.GA.hTini);
            delete(handles.GA.hTiniE);
            handles.GA = rmfield(handles.GA, 'hTfin');
            handles.GA = rmfield(handles.GA, 'hTfinE');
            handles.GA = rmfield(handles.GA, 'hTini');
            handles.GA = rmfield(handles.GA, 'hTiniE');
    end
    handles.GA.par(N).T = 0;
    handles.GA.par(N).F = 0;
    handles.GA.par(N).R = 0;
    set(handles.GA.checkbox(N,1), 'Value', 0);
    set(handles.GA.checkbox(N,2), 'Value', 0); 
    set(handles.GA.checkbox(N,3), 'Value', 0); 
    %choose corresponding variable
    s = ['{'];
    for i=1:handles.GA.variables
    s = [s '''' handles.GA.var(i).id ''''];
    if i~=handles.GA.variables
        s = [s ','];
    end
    end
    s =[s '}'];
    hrv_text = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 'VAR ',...
   'Units', 'normalized', ...
   'Position',[0.46 0.85-N*0.06 0.07 0.05]);
    hrv = uicontrol(handles.GA.hparameters,'Style','popupmenu',...
   'String', eval(s),...
   'Units', 'normalized', ...
   'BackgroundColor', 'w',...
   'Position', [0.54 0.87-N*0.06 0.1 0.04]);

    hUB_text = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 'Upper boundary: ',...
   'Units', 'normalized', ...
   'Position',[0.64 0.85-N*0.06 0.18 0.05]);
    
    hrv_UB = uicontrol(handles.GA.hparameters,'Style','edit',...
   'BackgroundColor', 'w',...
   'Units', 'normalized', ...
   'Position', [0.82 0.87-N*0.06 0.07 0.04]);

    hLog = uicontrol(handles.GA.hparameters,'Style','text',...
   'String', 'Log ',...
   'Units', 'normalized', ...
   'Position',[0.91 0.85-N*0.06 0.04 0.05]);
    hLog_ch = uicontrol(handles.GA.hparameters,'Style','checkbox',...
      'Units', 'normalized', ...
      'Value',0,...
      'Position',[0.96 0.86-N*0.06 0.05 0.05]);
    handles.GA.Hramp(N,:) = [hrv_text, hrv, hUB_text, hrv_UB, hLog, hLog_ch];
end

   guidata(handles.GA.hTGAga, handles);
   % update handles-data
guidata(hObject,handles);
end

% function massfraction_cb(hObject, eventdata, N)
% handles = guidata(hObject);
% 
% handles.GA.par(N).M = get(hObject, 'Value');
% if get(hObject, 'Value')==1
%     handles.GA.Mindex=N;
%     if handles.GA.par(N).T==1
%             delete(handles.GA.hTfin);
%             delete(handles.GA.hTfinE);
%             delete(handles.GA.hTini);
%             delete(handles.GA.hTiniE);
%             handles.GA = rmfield(handles.GA, 'hTfin');
%             handles.GA = rmfield(handles.GA, 'hTfinE');
%             handles.GA = rmfield(handles.GA, 'hTini');
%             handles.GA = rmfield(handles.GA, 'hTiniE');
%     end
%     handles.GA.par(N).T = 0;
%     handles.GA.par(N).F = 0;
%     handles.GA.par(N).R = 0;
%     set(handles.GA.checkbox(N,1), 'Value', 0);
%     set(handles.GA.checkbox(N,2), 'Value', 0); 
%     set(handles.GA.checkbox(N,3), 'Value', 0); 
% end
% handles.GA.hMassFraction = figure('Visible','off',...
%       'Name', 'Mass fraction', ...
%       'Menubar', 'none', ...
%       'NumberTitle', 'off', ...
%       'Toolbar', 'none', ...
%       'Color', get(0,'defaultuicontrolbackgroundcolor'), ...
%       'Units', 'normalized', ...
%       'Position',[0.4,0.6,0.25,0.13]);
% guidata(handles.GA.hTGAga,handles);
% handles = addMassF(handles,N);
% 
% set(handles.GA.hMassFraction, 'Visible', 'on');
% guidata(hObject, handles);
% uiwait(handles.GA.hMassFraction);
%    % update handles-data
%    handles = guidata(handles.GA.hMassFraction);
%    guidata(handles.GA.hTGAga, handles);
%    % update handles-data
%    close(handles.GA.hMassFraction);
% guidata(hObject,handles);
% end

% function handles = addMassF(handles,N)
% 
% hMF = uicontrol(handles.GA.hMassFraction,'Style','text',...
%    'String', 'Corresponding variables',...
%    'Units', 'normalized', ...
%    'Position',[0.1 0.9 0.5 0.1]);
% if isfield(handles.GA, 'variables')
%     s = ['{'];
%     for i=1:handles.GA.variables
%     s = [s '''' handles.GA.var(i).id ''''];
%     if i~=handles.GA.variables
%         s = [s ','];
%     end
%     end
%     s =[s '}'];
%     hMFE = uicontrol(handles.GA.hMassFraction,'Style','listbox',...
%    'String', eval(s),...
%    'Units', 'normalized', ...
%    'BackgroundColor', 'w',...
%    'Callback', @Mlist_cb, ...
%    'Min', 0,...
%    'Max',handles.GA.variables, ...
%    'Value', 1, ...
%    'Position', [0.1 0.1 0.3 0.8]);
%     handles.GA.Hres(N,:) = [hMF, hMFE];
% else
%     
%     hMFE = uicontrol(handles.GA.hMassFraction,'Style','edit',...
%    'BackgroundColor', 'w',...
%    'Units', 'normalized', ...
%    'Callback', {@massI_cb}, ...
%    'Position', [0.1 0.8 0.3 0.1]);
%     handles.GA.Hres(N,:) = [hMF, hMFE];
% end
% hMoisttext = uicontrol(handles.GA.hMassFraction,'Style','text', ...
%     'String', 'Moisture-%', ...
%     'Units', 'normalized', ...
%     'Position', [0.45 0.7 0.2 0.15]);
% handles.GA.hMoist = uicontrol(handles.GA.hMassFraction,'Style','edit',...
%    'BackgroundColor', 'w',...
%    'String', '0', ...
%    'Units', 'normalized', ...
%    'Position', [0.45 0.4 0.2 0.15]);
% 
% handles.GA.hOKMass = uicontrol(handles.GA.hMassFraction,'Style','pushbutton',...
%    'String', 'OK', ...
%    'Units', 'normalized', ...
%    'Callback', {@massOk_cb}, ...
%    'Position', [0.75 0.5 0.15 0.15]);
% 
% guidata(handles.GA.hMassFraction, handles);
% handles = guidata(handles.GA.hMassFraction);
% end

% function Mlist_cb(hObject, eventdata)
% handles = guidata(hObject);
% val = get(hObject, 'Value');
% handles.GA.massfraction = val;
% guidata(handles.GA.hMassFraction, handles);
% end
% 
% function massI_cb(hObject, eventdata)
% handles = guidata(hObject);
% str = get(hObject, 'String');
% handles.GA.massfraction = str2double(str);
% guidata(handles.GA.hMassFraction,handles);
% end
% 
% function massOk_cb(hObject, eventdata)
% handles = guidata(hObject);
% str= get(handles.GA.hMoist, 'String');
% handles.GA.moisture = str2double(str)./100;
% guidata(handles.GA.hMassFraction,handles);
% uiresume(handles.GA.hMassFraction);
% end

function okpar_cb(hObject, eventdata,type)
handles = guidata(hObject);
    
for i =1:handles.GA.parameters
if ~any(cell2mat(get(handles.GA.checkbox(i,:), 'Value')))
    msgbox('You must define parameters.');
    return
end
end

if exist('hTfinE')
if strcmp(get(handles.GA.hTfinE, 'String'), '') && any(cell2mat(get(handles.GA.checkbox(:,1), 'Value')))
    msgbox('Don''t forget to add the final temperature.');
    return
end
end

%-----------------------------------------------------
%new parameter structure

handles.GA.Par_struct = []; %index, values
handles.GA.Par_struct.nu_fuel.index = [];
handles.GA.Par_struct.ramp_value.index = [];
for i = 1:handles.GA.parameters

    if get(handles.GA.checkbox(i,1), 'Value')==1 %t(T_MAX)
       Tfin = get(handles.GA.hTfinE, 'String');
       Tini = get(handles.GA.hTiniE, 'String');
       TRamp = str2double(Tfin)-str2double(Tini);
        
       handles.GA.Par_struct.Tmax.index = i; %can be only one
       handles.GA.Par_struct.Tmax.values = TRamp; % temperature interval
     
    elseif get(handles.GA.checkbox(i,2), 'Value')==1 % t_fin
       handles.GA.Par_struct.tfin.index = i; %can be only one
       handles.GA.Par_struct.tfin.values = handles.GA.Par_struct.Tmax.values;
    elseif get(handles.GA.checkbox(i,3), 'Value')==1 %nu_fuel
        len = length(handles.GA.Par_struct.nu_fuel.index);
        handles.GA.Par_struct.nu_fuel.index(len+1) = i;
        handles.GA.Par_struct.nu_fuel.values(len+1) = handles.GA.par(i).res;
%     elseif get(handles.GA.checkbox(i,4), 'Value')==1 %mass fraction
%         len = length(handles.GA.Par_struct.mass_fraction.index);
%         handles.GA.Par_struct.mass_fraction.index(len+1) = i;
%         handles.GA.Par_struct.mass_fraction.M(len+1) = handles.GA.moisture;
    elseif get(handles.GA.checkbox(i,4), 'Value')==1 %ramp_value
        %this value hast to be greater or equal than the reference value
        rv  = get(handles.GA.Hramp(i,2),'Value');
        UB_val = get(handles.GA.Hramp(i,4),'String');
        UB = str2double(UB_val);
        Log = get(handles.GA.Hramp(i,6),'Value');
        len = length(handles.GA.Par_struct.ramp_value.index);
        handles.GA.Par_struct.ramp_value.index(len+1) = i;
        handles.GA.Par_struct.ramp_value.values(len+1:len+3) = [UB,rv, Log];
        %in groups of three: UB, rv,Log, UB, rv, Log...etc.
    end
end

%-----------------------------------------------------


%     TRamp = 0;
%     k=0;
%     NR=0;
%     for i =1:handles.GA.parameters
%     %value of TRamp
%     if get(handles.GA.checkbox(i,1), 'Value')==1
%     Tfin = get(handles.GA.hTfinE, 'String');
%     Tini = get(handles.GA.hTiniE, 'String');
%     TRamp = str2double(Tfin)-str2double(Tini);
%     end
%     NR = NR + handles.GA.par(i).R; % number of residue variables
%     end
%     if isfield(handles.GA, 'massfraction')
%     NM = length(handles.GA.massfraction);
%     end
%     handles.GA.P(5:(4+NR))=0;
%    % if NR ~= 0
%        handles.GA.P(5)=NR; 
%    % end
%     
%     
% for i =1:handles.GA.parameters
%     
%     % if ramp T is parameter
%     if get(handles.GA.checkbox(i,1), 'Value')==1
%     
%     handles.GA.P(1)=1;
%     handles.GA.P(2)=TRamp;  
%     %if final time is parameter
%     elseif get(handles.GA.checkbox(i,2), 'Value')==1
%     handles.GA.P(3)=1;
%     handles.GA.P(4)=TRamp;
%     %if residue mass(es) is variable
%     elseif get(handles.GA.checkbox(i,3), 'Value')==1
%     k=k+1;
%     handles.GA.P(5+k) = handles.GA.par(i).res;
%     %if mass fraction is variable
%     else
%     handles.GA.P(6+NR:5+NR+NM) = handles.GA.massfraction;
%     end
% end
%     
% if handles.GA.Tindex ~= 0
%     if handles.GA.Findex ~= 0
%        handles.GA.Pindex = [handles.GA.Tindex, handles.GA.Findex, handles.GA.Rindex, handles.GA.Mindex];
%     else
%        handles.GA.Pindex = [handles.GA.Tindex, handles.GA.Rindex, handles.GA.Mindex]; 
%     end
% else
%     handles.GA.Pindex = [handles.GA.Rindex, handles.GA.Mindex];
% end
% [p, handles.GA.parInd] = sort(handles.GA.Pindex);
guidata(hObject,handles);
uiresume(handles.GA.hparameters);
end

%-----------------

function GA_param_cb(hObject, eventdata,type) %GA parameters

handles = guidata(hObject);
if isfield(handles.GA, 'variables')
handles = Add_GA_param(handles, type);

guidata(handles.GA.hTGAga, handles);
else
    msgbox('Choose first the experimental data.');
end
end

function initialEstimates_cb(hObject, eventdata)
    handles = guidata(hObject);
    handles.GA.estimates = zeros(1,handles.GA.variables);
   if isfield(handles.GA, 'variables')==0
    msgbox('You must determine variables');
    return;
   end
    handles.GA.hestimates = figure('Visible','off',...
      'Name', 'Initial estimates', ...
      'Menubar', 'none', ...
      'NumberTitle', 'off', ...
      'Toolbar', 'none', ...
      'Color', get(0,'defaultuicontrolbackgroundcolor'), ...
      'Units', 'normalized', ...
      'Position',[0.4,0.6,0.2,0.5]);
    guidata(hObject, handles);
    handles = AddEst(handles);
    set(handles.GA.hestimates, 'Visible', 'on');
    guidata(handles.GA.hestimates,handles);
    uiwait(handles.GA.hestimates);
    handles = guidata(handles.GA.hestimates);
    guidata(handles.GA.hTGAga, handles);
    close(handles.GA.hestimates);
end

function handles = AddEst(handles)
htitle = uicontrol(handles.GA.hestimates, 'Style', 'text', ...
        'String', 'Set initial estimates', ...
        'Units', 'normalized', ...
        'Position', [0.1 0.8 0.7 0.1]);
    
    k=0.055;
    for i=1:handles.GA.variables
        hvar = uicontrol(handles.GA.hestimates, 'Style', 'text', ...
        'String', handles.GA.var(i).id, ...
        'Units', 'normalized', ...
        'Position', [0.2 0.8-i*k-0.004 0.2 0.05]);
        hset = uicontrol(handles.GA.hestimates, 'Style', 'edit', ...
            'String', handles.GA.limits(1,i), ...
           'BackgroundColor', 'w',...
            'Units', 'normalized', ...
            'Callback', {@est_cb,i}, ...
            'Position', [0.5 0.8-i*k 0.2 0.05]); 
    end
    hok = uicontrol(handles.GA.hestimates, 'Style', 'pushbutton', ...
        'String', 'OK', ...
        'Units', 'normalized', ...
        'Callback', @estok_cb, ...
        'Position', [0.35 0.05 0.3 0.05]);
end

function est_cb(hObject, eventdata, N)
    handles = guidata(hObject);
    estimate = str2double(get(hObject, 'String'));
    handles.GA.estimates(N)=estimate;
    guidata(hObject, handles);
end

function estok_cb(hObject, eventdata)
handles = guidata(hObject);
for i = 1:handles.GA.variables
   if handles.GA.estimates(i)==0
      handles.GA.estimates(i)=handles.GA.limits(1,i); 
   end
end
guidata(hObject, handles);
uiresume(handles.GA.hestimates);
end

%------------------
function cancelpar_cb(hObject, eventdata)
handles = guidata(hObject);
uiresume(handles.GA.hparameters);
end

%function grGA_cb(hObject, eventdata)
%handles = guidata(hObject);
%if isfield(handles.GA, 'bestChrom')==0
%    msgbox('You must first use the genetic algorithm');
%    return;
%end
%n=0;
%handles.template = 'bestInput.fds';
%for i=1:length(handles.GA.data)
%    if (handles.GA.data(i).check)
%        n = n+1;
%        data(n).Type = 'TGA';
%        data(n).T = handles.GA.data(i).temperature;
%       data(n).Time = handles.GA.data(i).time;
%        d = data(n).Time(2)-data(n).Time(1);
%        data(n).M = handles.GA.data(i).TGA;
%        data(n).dMdt = gradient(data(n).M./max(data(n).M),d);
%        data(n).Rate = handles.GA.data(i).rate;
%    end
%end
%handles.GA.best = gradientMethod(handles.GA.bestChrom,handles.GA.oV, data, handles.GA.LogScaling, handles.template, handles.FdsExe, handles.GA.P, handles.GA.parInd);
%guidata(hObject, handles);
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runGA(hObject, eventdata, type)
handles = guidata(hObject);
if isfield(handles.GA, 'variables')==0
    msgbox('You must determine variables');
    return;
end

if ~exist('recombin.m','file')
   msgbox('Add genetic algorithm toolbox to your Matlab path');
    return
end

if ~isfield(handles.GA, 'Par_struct')
    handles.GA.Par_struct = [];
    handles.GA.Par_struct.nu_fuel.index = [];
    handles.GA.Par_struct.ramp_value.index = [];
end

%create TGA structure
n = 0;
TGA = [];
% if isfield(handles.GA, 'parameters')== 0
%     handles.GA.P = [0 0 0 0];
%     handles.GA.parInd = [];
% end
if strcmp(type,'TGA')
for i=1:length(handles.GA.data)
    if (handles.GA.data(i).check)
        n = n+1;
        data(n).Type = 'TGA';
        data(n).T = handles.GA.data(i).temperature;
        data(n).Time = handles.GA.data(i).time;
        d = data(n).Time(2)-data(n).Time(1);
        data(n).M = handles.GA.data(i).TGA;
        data(n).dMdt = gradient(data(n).M./max(data(n).M),d);
        data(n).Rate = handles.GA.data(i).rate;
    end
end
else
    ox = 0; % if oxygen limited cone exists, ox = 1
    for i=1:length(handles.GA.data)
    if (handles.GA.data(i).check)
        %HRR or MLR ??
        X = get(handles.hDataType, 'Value');
        Mass = [];
        dataType = 1;
        
                
        switch X
            case 1
                %Mass_h = handles.GA.data(i).HRR;
                dataType = 1;
            case 2
                %Mass_m = handles.GA.data(i).MLR;
                dataType = 2;
        end
        
        if strcmp(handles.GA.data(i).gas,'Air')
            Mass_h = handles.GA.data(i).HRR;
            Mass_m = handles.GA.data(i).MLR;
            if (handles.Options.NFilter(1)>1)
            y(:,1)=filtNs(Mass_h,handles.Options.NFilter(1));
            y(:,2)=filtNs(Mass_m,handles.Options.NFilter(1));
            else
            y(:,1) = Mass_h;
            y(:,2) = Mass_m;
            end
        elseif strcmp(handles.GA.data(i).gas,'N2')
            ox = 1;
            if (handles.Options.NFilter(1)>1)
               y(:,1) =  filtNs(handles.GA.data(i).MLR,handles.Options.NFilter(1));
            else
               y(:,1) = handles.GA.data(i).MLR;
            end
        else
            msgbox('Error');
            
        end
        
        n = n+1;
        data(n).Type = 'Cone';
        data(n).T = handles.GA.data(i).time;
        data(n).M = y;
        data(n).A = handles.GA.data(i).A;
        data(n).dataType = dataType;
        data(n).gas = handles.GA.data(i).gas;
    end
    handles.GA.moisture = 0;
    end
end

if strcmp(data(1).Type,'Cone') && ox == 1 %write oxygen cone fds-input
    new_file = change_input_ambient(handles.template);
end

if n==0
    msgbox('Choose at least one TGA data series');
    return;
end

%parameters of algorithm
alg_param = [];
if isfield(handles.GA, 'weight_index')
   if ~isempty(handles.GA.weight_index(:,1))
       alg_param.index = handles.GA.weight_index;
       alg_param.weight = handles.GA.exp_weight;
   end
end
alg_param.NIND = handles.GA.NIND;
alg_param.GGAP = handles.GA.GGAP;
alg_param.XOVR = handles.GA.XOVR;
alg_param.MUTR = handles.GA.MUTR;
alg_param.MAXGEN = handles.GA.MAXGEN;
alg_param.INSR = handles.GA.INSR;
alg_param.SUBPOP = handles.GA.SUBPOP;
alg_param.MIGR = handles.GA.MIGR;
alg_param.MIGGEN = handles.GA.MIGGEN;

guidata(hObject, handles);

% Best Chrom - log
fid = fopen('bestChrom.csv', 'w');
fprintf(fid, '%s\n', 'Best Chrom of the generation');
fprintf(fid, '%s\n', '');
str = 'Generation, ';
for i = 1:handles.GA.variables
    str = [str, handles.GA.var(i).id, ','];
end
if ~isempty(handles.GA.Par_struct.ramp_value.index)
   for i = 1:length(handles.GA.Par_struct.ramp_value.index)
       str = [str, 'Ramp value ', mat2str(i),','];
   end
end
str = [str, 'Fitness value'];
fprintf(fid, '%s\n', str);

fclose(fid);

[handles.GA.bestChrom, handles.GA.oV, parameters] = fdsga(handles.GA.estimates, data, handles.template, handles.FdsExe, handles.GA.weights, handles.GA.Aindex, handles.GA.limits, handles.GA.LogScaling, handles.GA.Par_struct, handles.GA.moisture,alg_param);
 fid = fopen('bestChrom.csv', 'a');
% fprintf(fid, '%s\n', 'Parameters');
% fprintf(fid, '%s\n', '');
% for i = 1:handles.GA.variables
%     fprintf(fid, '%s\t', handles.GA.var(i).id);
%     fprintf(fid, '%f\n', handles.GA.bestChrom(i));
% end
% for i = 1:length(parameters)
%     fprintf(fid, '%s\t', ['PAR', mat2str(i)]);
%     fprintf(fid, '%f\n', parameters(i));
% end
% fclose(fid);
open('bestChrom.csv');
guidata(handles.GA.hTGAga,handles);
end
