function pyroplot

% TODO
% 
% - cone data selection of HRR, MLR columns.


% Documentation
%
% handles.EXPDATA().temperature       Temperature (C)
% handles.EXPDATA().TGA               TGA mass ratio (0...1)
% handles.EXPDATA().DSC               DSC power (W/kg)
% handles.EXPDATA().EXO               Direction of exothermic peaks in original data file
% handles.EXPDATA().ConeHRR           Cone heat release rate(W/m2)
% handles.EXPDATA().ConeMLR           Cone mass loss rate (kg/s/m2)
% handles.EXPDATA().type              Type of data ('TGA','DSC', 'Cone')
% handles.EXPDATA().name              Data name
% handles.EXPDATA().path              Full file name of data
% handles.EXPDATA().check             Plot options
% handles.EXPDATA().pair              Index of data pair for STA data
% handles.EXPDATA().p_set_l           Settings for plotting (left)
% handles.EXPDATA().p_set_r           Settings for plotting (right)

%  Initialize and hide the GUI as it is being constructed.
hPyroPlot = figure(...
   'Visible', 'off',...
   'Name', 'Pyro Plot',...
   'MenuBar','none', ...
   'Toolbar','none', ...
   'NumberTitle', 'off', ...
   'HandleVisibility','callback', ...
   'Units', 'normalized', ...
   'Position',[0.1,0.2,0.8,0.7],...
   'Color', get(0,'defaultuicontrolbackgroundcolor'));

% Add menus
handles.hPyroPlot = hPyroPlot;
%Initialize
handles = InitializePyroPlot(handles);
%Add menus
handles = AddMenus(handles);
%Add static titles and panels
handles = AddStatic(handles);

% Move the GUI to the center of the screen.
movegui(handles.hPyroPlot,'center')
% Make the GUI visible.
set(handles.hPyroPlot,'Visible','on');

guidata(hPyroPlot,handles)

return
end

% CALLBACKS

%-------------------------------------------------------
% Add Menus
%-------------------------------------------------------
function handles = AddMenus(handles)

% add file menu
handles.hMenu.hFileMenu = uimenu(...       
                        'Parent',handles.hPyroPlot,...
                        'HandleVisibility','callback', ...
                        'Label','&File');
% add 'Load project' menu item
handles.hMenu.hLoadProjectitem  =   uimenu(...       
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','Load &Project',...
                        'HandleVisibility','callback', ...
                        'Callback', @loadProjectitemCallback);
                    
% add 'Save Project' menu item
handles.hMenu.hSaveProjectitem  =   uimenu(... 
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','&Save Project',...
                        'HandleVisibility','callback', ...
                        'Callback', @saveProjectitemCallback);

% add 'Save As Project' menu item
handles.hMenu.hSaveAsProjectitem  =   uimenu(...      
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','Save &As Project',...
                        'HandleVisibility','callback', ...
                        'Callback', @saveAsProjectitemCallback);

% add 'Read STA' menu item
handles.hMenu.hReadDTAitem  =   uimenu(...       
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','&Read STA',...
                        'Separator','on',...
                        'HandleVisibility','callback', ...
                        'Callback', {@readDATACallback,'DTA'});
                    
% add 'Read TGA' menu item
handles.hMenu.hReadTGAitem  =   uimenu(...       
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','Read &TGA',...
                        'Separator','on',...
                        'HandleVisibility','callback', ...
                        'Callback', {@readDATACallback,'TGA'});
handles.hMenu.hReadTGAFDSitem  =   uimenu(...       
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','Read T&GA (FDS)',...
                        'HandleVisibility','callback', ...
                        'Callback', {@readDATACallback,'TGAFDS'});
% add 'Read DSC' menu item
handles.hMenu.hReadDSCitem  =   uimenu(...       
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','Read &DSC',...
                        'HandleVisibility','callback', ...
                        'Callback', {@readDATACallback,'DSC'});
% add 'Read Cone'
handles.hMenu.hReadConeitem  =   uimenu(...      
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','Read &Cone type 1',...
                        'HandleVisibility','callback', ...
                        'Callback', {@readDATACallback,'Cone'});
handles.hMenu.hReadConeitem2  =   uimenu(...
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','Read Cone type &2',...
                        'HandleVisibility','callback', ...
                        'Callback', {@readDATACallback,'Cone2'});
handles.hMenu.hReadConeitem  =   uimenu(...      
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','Read Cone (&FDS)',...
                        'HandleVisibility','callback', ...
                        'Callback', {@readDATACallback,'ConeFDS'});           
% add 'Close' menu item
handles.hMenu.hCloseMenuitem  =  uimenu(...       
                        'Parent',handles.hMenu.hFileMenu,...
                        'Label','Clos&e',...
                        'Separator','on',...
                        'HandleVisibility','callback', ...
                        'Callback', @hCloseMenuitemCallback);
                    
% add Settings menu
handles.hMenu.hSettingsMenu      =   uimenu(...       
                        'Parent',handles.hPyroPlot,...
                        'HandleVisibility','callback', ...
                        'Label','&Settings');
  
% add Settings menu items
handles.hMenu.hFilter  =   uimenu(...       
                        'Parent',handles.hMenu.hSettingsMenu,...
                        'Label','&Filter',...
                        'HandleVisibility','callback', ...
                        'Callback', @filter_Callback);

handles.hMenu.hPlot_set  =   uimenu(...       
                        'Parent',handles.hMenu.hSettingsMenu,...
                        'Label','&Plot Settings',...
                        'HandleVisibility','callback', ...
                        'Callback', @plot_set_Callback);
                    
handles.hMenu.hPlot_set  =   uimenu(...       
                        'Parent',handles.hMenu.hSettingsMenu,...
                        'Label','&Restore default plot settings',...
                        'HandleVisibility','callback', ...
                        'Callback', @restore_set_Callback);
                    
                    
% add Tools menu
handles.hMenu.hToolsMenu      =   uimenu(...       
                        'Parent',handles.hPyroPlot,...
                        'HandleVisibility','callback', ...
                        'Label','&Tools');
  
% add Tools menu items
handles.hMenu.hSubPlot  =   uimenu(...       
                        'Parent',handles.hMenu.hToolsMenu,...
                        'Label','&SubPlot',...
                        'HandleVisibility','callback', ...
                        'Callback', @subplot_Callback);
handles.hMenu.hIntegration  =   uimenu(...       
                        'Parent',handles.hMenu.hToolsMenu,...
                        'Label','&DSC-integration',...
                        'HandleVisibility','callback', ...
                        'Callback', @DSC_Integration_Callback);
handles.hMenu.hEstimates  =   uimenu(...       
                        'Parent',handles.hMenu.hToolsMenu,...
                        'Label','&Initial estimates',...
                        'HandleVisibility','callback', ...
                        'Callback', @Estimates_Callback);
handles.hMenu.hGA = uimenu(...       
                        'Parent',handles.hMenu.hToolsMenu,...
                        'Label','&Genetic Algorithm (TGA)',...
                        'HandleVisibility','callback', ...
                        'Callback', @GA_Callback);
handles.hMenu.hGAC = uimenu(...       
                        'Parent',handles.hMenu.hToolsMenu,...
                        'Label','Genetic Algorithm (&Cone)',...
                        'HandleVisibility','callback', ...
                        'Callback', @GAC_Callback);   
handles.hMenu.hRP = uimenu(...       
                        'Parent',handles.hMenu.hToolsMenu,...
                        'Label','Reaction parameters &without estimation',...
                        'HandleVisibility','callback', ...
                        'Callback', @PR_Callback);   
handles.hMenu.hFDS = uimenu(...       
                        'Parent',handles.hMenu.hToolsMenu,...
                        'Label','&Run FDS',...
                        'HandleVisibility','callback', ...
                        'Callback', @fds_Callback);                
end %AddMenus

%----------------------------------------------------------
% Add static titles, figure panels and buttons
%
function handles = AddStatic(handles)
%
%add static text titles             
handles.hStatic.htitle1 = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','MATERIAL',...
   'Units', 'normalized', ...
   'Position',[0.025 0.92 0.05 0.0286]);
         
handles.hStatic.htitle2 = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','TYPE',...
   'Units', 'normalized', ...
   'Position',[0.10 0.92 0.05 0.0286]);

handles.hStatic.htitle2 = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','RATE',...
   'Units', 'normalized', ...
   'Position',[0.1545 0.92 0.05 0.0286]);

handles.hStatic.htitle3 = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','GAS',...
   'Units', 'normalized', ...
   'Position',[0.21 0.92 0.05 0.0286]);
handles.hStatic.htitle4 = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','|',...
   'Units', 'normalized', ...
   'Position',[0.265 0.92 0.007 0.0286]);
handles.hStatic.htitle5 = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','|',...
   'Units', 'normalized', ...
   'Position',[0.48 0.92 0.007 0.0286]);

%figure 1

handles.hStatic.htitlefig1 = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','FIGURE 1',...
   'Units', 'normalized', ...
   'Position',[0.27 0.957 0.07 0.0286]);
handles.hStatic.hpopupfig1 = uicontrol(handles.hPyroPlot, ...
   'Style','popupmenu',...
   'String',{'Temperature','Time'},...
   'BackgroundColor', 'w',...
   'Units', 'normalized', ...
   'Position',[0.348 0.957 0.1 0.0357],...
   'Callback',{@hpopfig_cb,1});

handles.hStatic.hleftfig1 =uicontrol(handles.hPyroPlot,'Style','text',...
   'String','left',...
   'Units', 'normalized', ...
   'Position',[0.27 0.92 0.03 0.0286]);
handles.hStatic.hrightfig1=uicontrol(handles.hPyroPlot,'Style','text',...
   'String','right',...
   'Units', 'normalized', ...
   'Position',[0.3 0.92 0.03 0.0286]);
            
%figure 2
handles.hStatic.htitlefig2 = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','FIGURE 2',...
   'Units', 'normalized', ...
   'Position',[0.482 0.957 0.07 0.0286]);
handles.hStatic.hpopupfig2 = uicontrol(handles.hPyroPlot, ...
   'Style','popupmenu',...
   'String',{'Temperature','Time'},...
   'BackgroundColor', 'w',...
   'Units', 'normalized', ...
   'Position',[0.564 0.957 0.1 0.0357],...
   'Callback',{@hpopfig_cb,2});
handles.hStatic.hleftfig2 =uicontrol(handles.hPyroPlot,'Style','text',...
   'String','left',...
   'Units', 'normalized', ...
   'Position',[0.485 0.92 0.03 0.0286]);
handles.hStatic.hrightfig2=uicontrol(handles.hPyroPlot,'Style','text',...
   'String','right',...
   'Units', 'normalized', ...
   'Position',[0.515 0.92 0.03 0.0286]);
    
%Panel for FIGURE 1
handles.hStatic.hf1 = uipanel(handles.hPyroPlot, 'Position',[.71 .7 .25 .25]);
handles.hStatic.htitlef1 = uicontrol(handles.hStatic.hf1,'Style','text',...
   'String','FIGURE 1',...
   'Units', 'normalized', ...
   'Position',[0.04 0.85 0.28 0.11]);
handles.hStatic.htextf1 =  uicontrol(handles.hStatic.hf1,'Style','text',...
   'String','Sync Axis',...
   'Units', 'normalized', ...
   'Position',[0.04 0.71 0.28 0.11]);
handles.hStatic.hRL1 = uicontrol(handles.hStatic.hf1,'Style','pushbutton','String','R to L',...
   'Units', 'normalized', ...
   'Position',[0.04 0.57 0.2 0.11],...
   'Callback',{@RL1_Callback});
handles.hStatic.hLR1 = uicontrol(handles.hStatic.hf1,'Style','pushbutton','String','L to R',...
   'Units', 'normalized', ...
   'Position',[0.28 0.57 0.2 0.11],...
   'Callback',{@LR1_Callback});
handles.hStatic.hylabel1l = uicontrol(handles.hStatic.hf1,'Style','text',...
   'String','Y-Label Left',...
   'Units', 'normalized', ...
   'Position',[0.04 0.4 0.28 0.11]);
handles.hStatic.hlabel1l = uicontrol(handles.hStatic.hf1,'Style','edit',...
   'String','TGA (%)',...
   'BackgroundColor', 'w',...
   'CallBack', {@L1l_Callback}, ...
   'Units', 'normalized', ...
   'Position',[0.04 0.3 0.6 0.11]);
handles.hStatic.hylabel1r = uicontrol(handles.hStatic.hf1,'Style','text',...
   'String','Y-Label Right',...
   'Units', 'normalized', ...
   'Position',[0.04 0.17 0.28 0.11]);
handles.hStatic.hlabel1r = uicontrol(handles.hStatic.hf1,'Style','edit',...
   'String','DSC (kW/kg)',...
   'BackgroundColor', 'w',...
   'CallBack', {@L1r_Callback}, ...
   'Units', 'normalized', ...
   'Position',[0.04 0.075 0.6 0.11]);

%Panel for FIGURE 2
handles.hStatic.hf2 = uipanel(handles.hPyroPlot, 'Position',[.71 .4 .25 .25]);
handles.hStatic.htitlef2 = uicontrol(handles.hStatic.hf2,'Style','text',...
   'String','FIGURE 2',...
   'Units', 'normalized', ...
   'Position',[0.04 0.85 0.28 0.11]);
handles.hStatic.htextf2 =  uicontrol(handles.hStatic.hf2,'Style','text',...
   'String','Sync Axis',...
   'Units', 'normalized', ...
   'Position',[0.04 0.71 0.28 0.11]);
handles.hStatic.hRL2 = uicontrol(handles.hStatic.hf2,'Style','pushbutton','String','R to L',...
   'Units', 'normalized', ...
   'Position',[0.04 0.57 0.2 0.11],...
   'Callback',{@RL2_Callback});
handles.hStatic.hLR2 = uicontrol(handles.hStatic.hf2,'Style','pushbutton','String','L to R',...
   'Units', 'normalized', ...
   'Position',[0.28 0.57 0.2 0.11],...
   'Callback',{@LR2_Callback});
handles.hStatic.hylabel2l = uicontrol(handles.hStatic.hf2,'Style','text',...
   'String','Y-Label Left',...
   'Units', 'normalized', ...
   'Position',[0.04 0.4 0.28 0.11]);
handles.hStatic.hlabel2l = uicontrol(handles.hStatic.hf2,'Style','edit',...
   'String','TGA (%)',...
   'BackgroundColor', 'w',...
   'CallBack', @L2l_Callback, ...
   'Units', 'normalized', ...
   'Position',[0.04 0.3 0.6 0.11]);
handles.hStatic.hylabel2r = uicontrol(handles.hStatic.hf2,'Style','text',...
   'String','Y-Label Right',...
   'Units', 'normalized', ...
   'Position',[0.04 0.17 0.28 0.11]);
handles.hStatic.hlabel2r = uicontrol(handles.hStatic.hf2,'Style','edit',...
   'String','DSC (kW/kg)',...
   'BackgroundColor', 'w',...
   'CallBack', @L2r_Callback, ...
   'Units', 'normalized', ...
   'Position',[0.04 0.075 0.6 0.11]);

align([handles.hStatic.htitle1],'Center','None');
set([handles.hPyroPlot,handles.hStatic.htitle1],'Units','normalized');

%delete line
handles.hStatic.hdelpanel = uipanel(handles.hPyroPlot, ...
   'Position',[.71 .2 .25 .12]);
handles.hStatic.hdeltitle = uicontrol(handles.hStatic.hdelpanel, ...
   'Style', 'text', ...
   'string', 'DELETE LINE', ...
   'Units', 'normalized', ...
   'Position', [0.04 0.71 0.3 0.24]);
handles.hStatic.hdeltitle2 = uicontrol(handles.hStatic.hdelpanel, ...
   'Style', 'text', ...
   'string', 'Enter line number', ...
   'Units', 'normalized', ...
   'Position', [0.04 0.35 0.36 0.24]);
handles.hStatic.hdelline = uicontrol(handles.hStatic.hdelpanel,...
   'Style','edit',...
   'String','line',...
   'BackgroundColor', 'w',...
   'Units', 'normalized', ...
   'Position',[0.04 0.11 0.2 0.24]);
handles.hStatic.hdelok = uicontrol(handles.hStatic.hdelpanel,...
   'Style','pushbutton','String','OK',...
   'Units', 'normalized', ...
   'Position',[0.4 0.11 0.2 0.24],...
   'Callback',{@delok_Callback});

%Plot -pushbutton (data is plotted only after this is clicked)
% handles.hStatic.hdelok = uicontrol(handles.hPyroPlot,...
%    'Style','pushbutton','String','PLOT',...
%    'Units', 'normalized', ...
%    'Position',[0.8 0.1 0.05 0.05],...
%    'Callback',{@plot_data_Callback});

end %AddStatic

%-------------------------------------------------------
% Load Project   
%-------------------------------------------------------
function loadProjectitemCallback(hObject, eventdata)
handles=guidata(hObject);
[file, path]=uigetfile('*.mat', 'Load Project', handles.var.save_path);    
if (~isequal(file,0) && ~isequal(path,0))
    handles.var.save_path = path;
    A=load(fullfile(path, file));
    B=loadobj(A);
    struct=B.obj;
    handles.var.lines=length(struct);
    for k=1:length(struct)
       handles.EXPDATA(k)=struct(k);
       handles = update_line(handles,k);
       handles.Options.rb(k,1:3)=0;
    end
    handles = updatePlot(handles);
end
guidata(hObject,handles)
return
end

%at this moment 'save' and 'save as' are the same

%-------------------------------------------------------
% Save Project
%-------------------------------------------------------
function saveProjectitemCallback(hObject, eventdata)
handles=guidata(hObject);
[file, path]=uiputfile('*.mat', 'Save Project', handles.var.save_path);
if (~strcmp(file,'') && exist(path,'dir'))
   handles.var.save_path = path;
   obj = saveobj(handles.EXPDATA);
   save(fullfile(path, file), 'obj');
   guidata(hObject,handles)
end
end

%-------------------------------------------------------
% Save As
%-------------------------------------------------------
function saveAsProjectitemCallback(hObject, eventdata)
handles=guidata(hObject);
[file, path]=uiputfile('*.mat', 'Save As', handles.var.save_path);
if (file==0)
   return
end
if (~strcmp(file,''))
   if (exist(path,'dir'))
      handles.var.save_path = path;
      obj = saveobj(handles.EXPDATA);
      save(fullfile(path, file), 'obj');
      guidata(hObject,handles)
   end
end
end


%----------------------------------------------------------
% Read DATA
%
function readDATACallback(hObject,eventdata,DType)
%
handles = guidata(hObject);
lines_old = handles.var.lines;
handles = ppread(handles,DType);
if (handles.var.lines > lines_old)
   if strcmp(DType,'DTA')
      handles = update_line(handles,handles.var.lines-1);
   end
   handles = update_line(handles,handles.var.lines);
end
guidata(hObject,handles)
end

        
%-------------------------------------------------
%close callback
function hCloseMenuitemCallback(hObject, eventdata)
    close all;
end

%---------------------------------------------------
% Settings menu callbacks
%
% Filter
function filter_Callback(hObject, eventdata)
handles = guidata(hObject);

if ~exist('filtNs.m', 'file')
    msgbox('Add Matlab path to Palo\FIRAS\Tutkimus\Osatehtävä 1\PyroPlot\MatlabTools_Fire.');
    return;
end

filter_ok = 0;
while (~filter_ok)
   %open window to ask how many point average user wants
   p = {'\fontsize{12} Figure 1 Left',...
      '\fontsize{12} Figure 1 Right',...
      '\fontsize{12} Figure 2 Left',...
      '\fontsize{12} Figure 2 Right'};
   name = 'Filter Width (odd integers)';
   numlines = 1;
   for i = 1:4
      defaultanswer{i} = num2str(handles.Options.NFilter(i));
   end
   options.Interpreter='tex';
   answer=inputdlg(p,name,numlines,defaultanswer,options);
   if isempty(answer)
      return;
   end
   %
   for i = 1:4
      handles.Options.NFilter(i) = round(str2double(answer{i}));
   end
   if (any(mod(handles.Options.NFilter,2)==0))
      uiwait(msgbox('Filter widths must be odd integers','Error','error','modal'));
   else
      filter_ok = 1;
   end
end   
guidata(hObject,handles);
end

function plot_set_Callback(hObject, eventdata)
handles = guidata(hObject);

% user set defaults
 p = {'\fontsize{10} Position',...
      '\fontsize{10} FontSize - left',...
      '\fontsize{10} Line width - left',...
      '\fontsize{10} Line style - left (1 - solid, 2 - dash, 3 - dot, 4 - dashdot', ...
      '\fontsize{10} FontSize - right',...
      '\fontsize{10} Line width - right',...
      '\fontsize{10} Line style - right (1 - solid, 2 - dash, 3 - dot, 4 - dashdot'};
 name = 'Plot settings';
  
  numlines = 1;
  load plot_settings.mat
  
  defaultanswer{1} = num2str(settings.Position);
  defaultanswer{2} = num2str(settings.left.FontSize);
  defaultanswer{3} = num2str(settings.left.Width);
  defaultanswer{4} = num2str(settings.left.Style);
  defaultanswer{5} = num2str(settings.right.FontSize);
  defaultanswer{6} = num2str(settings.right.Width);
  defaultanswer{7} = num2str(settings.right.Style);
  
  options.Interpreter='tex';
  answer=inputdlg(p,name,numlines,defaultanswer,options);
   
  if isempty(answer)
     return;
  end
   
  str = str2mat(answer(1));
  k = 0;
  num = '';
   
  for i = 1:length(str)
    if ~isnan(str2double(str(i))) || strcmp(str(i), '.')
        num = [num, str(i)];
    else
        if ~isempty(num)
            k = k+1;
            settings.Position(k) = str2double(num);
            num = '';
        end
    end
  end
  if ~isempty(num)
      k = k+1;
      settings.Position(k) = str2double(num);
      clear num
  end
  settings.left.FontSize = round(str2double(answer{2}));
  settings.left.Width = round(str2double(answer{3}));
  settings.left.Style = round(str2double(answer{4}));
  settings.right.FontSize = round(str2double(answer{5}));
  settings.right.Width = round(str2double(answer{6}));
  settings.right.Style = round(str2double(answer{7}));
  save plot_settings.mat settings
  
guidata(hObject,handles);
end

% restore defaults
function restore_set_Callback(hObject, eventdata)
handles = guidata(hObject);

settings.Position = [0.1, 0.1, 0.65, 0.6];
settings.left.FontSize = 14;
settings.left.Width = 2;
settings.left.Style = 1;
settings.right.FontSize = 14;
settings.right.Width = 2;
settings.right.Style = 2;
save plot_settings.mat settings

msgbox('Default settings have been restored');

guidata(hObject,handles);
end
%------------------------------------------------
%Tools callbacks

function subplot_Callback(hObject, eventdata)
handles = guidata(hObject);
handles = plot_subplot(handles);
guidata(hObject,handles);
end

function DSC_Integration_Callback(hObject, eventdata)
handles = guidata(hObject);
handles = ppdscint(handles);
guidata(hObject,handles);
end

function Estimates_Callback(hObject, eventdata)
handles = guidata(hObject);
handles = findEstimates(handles);
guidata(hObject,handles);
end

function GA_Callback(hObject, eventdata)
handles = guidata(hObject);
handles = ga(handles, 'TGA');
guidata(hObject, handles);
end

function GAC_Callback(hObject, eventdata)
handles = guidata(hObject);
handles = ga(handles, 'Cone');
guidata(hObject, handles);
end

function PR_Callback(hObject, eventdata)
handles = guidata(hObject);
handles = find_reac_param(handles);
guidata(hObject, handles);
end

function fds_Callback(hObject, eventdata)
handles = guidata(hObject);
[v, cdir] = system('cd');
[fds_file, handles.template_path]=uigetfile('*.fds', 'Choose FDS Input File', handles.template_path);
handles.template = fullfile(handles.template_path, fds_file);
k=1;
for i=1:min(length(handles.template),length(cdir))
    if handles.template(i)==cdir(i)
        k=k+1;
    else
        break;
    end
end
l=0;
if length(cdir) > k
    l=1;
    for i=k:length(cdir)
        if cdir(i)=='\'
            l=l+1;
        end
    end
end
s = [];
for i=1:l
    s = [s '..\'];
end
if handles.template(k)=='\'
    k=k+1;
end
handles.template = handles.template(k:length(handles.template));
handles.template = ['"' s handles.template '"'];
[Exe_file, handles.FdsExe_path]=uigetfile('*.exe', 'Choose Executive File', handles.FdsExe_path);
handles.FdsExe = fullfile(handles.FdsExe_path, Exe_file);
k=1;
for i=1:min(length(handles.FdsExe),length(cdir))
    if handles.FdsExe(i)==cdir(i)
        k=k+1;
    else
        break;
    end
end
l=0;
if length(cdir) > k
    l=1;
    for i=k:length(cdir)
        if cdir(i)=='\'
            l=l+1;
        end
    end
end
s = [];
for i=1:l
    s = [s '..\'];
end
if handles.FdsExe(k)=='\'
    k=k+1;
end
handles.FdsExe = handles.FdsExe(k:length(handles.FdsExe));
handles.FdsExe = ['"' s handles.FdsExe '"'];
handles.err_file = 'fds.err';
handles.err_file = ['"' handles.err_file '"'];

sysstr=[handles.FdsExe ' ' handles.template ' 2> ' handles.err_file];
system(sysstr);
guidata(hObject, handles);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Callbacks of x-axis
function hpopfig_cb(hObject, eventdata,nFig)
handles = guidata(hObject);
val=get(hObject, 'Value');
str=get(hObject, 'String');
set(0, 'CurrentFigure', handles.hFigure(nFig));
switch str{val};
   case 'Temperature'
      handles.Options.Fig_xtype(nFig)=0;
      handles = updatePlot(handles);
   case 'Time'
      handles.Options.Fig_xtype(nFig)=1;
      handles = updatePlot(handles);
end
guidata(hObject,handles)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks of buttongroups
% If user wants TGA data as mass % or as its gradient
% or if user wants Cone data as HRR, MLR or EHC

function selcbk(hObject, eventdata,group)
handles = guidata(hObject);
N = group-1;
str=get(get(hObject, 'SelectedObject'), 'String');
[a,b]=size(handles.hRBvector);
switch str
   case '%'
      % Code for when '%' is selected.
      for i=1:a
         if get(hObject, 'SelectedObject')==handles.hRBvector(i,1+N*3)
            l=i; % l is linenumber of selected data serie
         end
      end
      handles.Options.rb(l,group)=0;
   case '%/s'
      % Code for when '%/s' is selected.
      for i=1:a
         if get(hObject, 'SelectedObject')==handles.hRBvector(i,2+N*3)
            l=i; % l is linenumber of selected data serie
         end
      end
      handles.Options.rb(l,group)=1;
   case '-%/s'
      for i=1:a
         if get(hObject, 'SelectedObject')==handles.hRBvector(i,3+N*3)
            l=i; % l is linenumber of selected data serie
         end
      end
      handles.Options.rb(l,group)=2;
   % cone cases
   case 'HRR'
      for i=1:a
         if get(hObject, 'SelectedObject')==handles.hRBvector(i,1+N*3)
            l=i; % l is linenumber of selected data serie
         end
      end
      handles.Options.rb(l,group)=0;
   case 'MLR'
      for i=1:a
         if get(hObject, 'SelectedObject')==handles.hRBvector(i,2+N*3)
            l=i; % l is linenumber of selected data serie
         end
      end
      handles.Options.rb(l,group)=1;
   case 'EHC'
      for i=1:a
         if get(hObject, 'SelectedObject')==handles.hRBvector(i,3+N*3)
            l=i; % l is linenumber of selected data serie
         end
      end
      handles.Options.rb(l,group)=2;
end
handles = updatePlot(handles);
guidata(hObject,handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback of checkboxes
function h_cb(hObject, eventdata,indx)
handles = guidata(hObject);
%size of hvector 
[a,b]=size(handles.hCBvector);
% which check box?
for i=1:a 
   if (hObject==handles.hCBvector(i,indx))
      l=i;
      break;
   end
end
%
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked
   handles.EXPDATA(l).check(indx) = 1;
   handles = updatePlot(handles);
else
   % Checkbox is not checked
   handles.EXPDATA(l).check(indx)=0;
end
guidata(hObject,handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Callbacks of Figure1 Panel

function RL1_Callback(hObject, eventdata)
handles = guidata(hObject);
set(handles.hAxes.ax12, 'Position', get(handles.hAxes.ax11, 'Position'));
guidata(hObject,handles);
end

function LR1_Callback(hObject, eventdata)
handles = guidata(hObject);
set(handles.hAxes.ax11, 'Position', get(handles.hAxes.ax12, 'Position'));
guidata(hObject,handles);
end

function L1l_Callback(hObject, eventdata)
handles = guidata(hObject);
user_entry = get(hObject, 'string');
handles.var.yLabel1Left= user_entry;
handles = updatePlot(handles);
guidata(hObject,handles)
end

function L1r_Callback(hObject, eventdata)
handles = guidata(hObject);
user_entry = get(hObject, 'string');
handles.var.yLabel1Right = user_entry;
handles = updatePlot(handles);
guidata(hObject,handles)
end

%Callbacks of Figure2 Panel

function RL2_Callback(hObject, eventdata)
handles = guidata(hObject);
set(handles.hAxes.ax22, 'Position', get(handles.hAxes.ax21, 'Position'));
guidata(hObject,handles);
end

function LR2_Callback(hObject, eventdata)
handles = guidata(hObject);
set(handles.hAxes.ax21, 'Position', get(handles.hAxes.ax22, 'Position'));
guidata(hObject,handles);
end

function L2l_Callback(hObject, eventdata)
handles = guidata(hObject);
user_entry = get(hObject, 'string');
handles.var.yLabel2Left= user_entry;
handles = updatePlot(handles);
guidata(hObject,handles)
end

function L2r_Callback(hObject, eventdata)
handles = guidata(hObject);
user_entry = get(hObject, 'string');
handles.var.yLabel2Right = user_entry;
handles = updatePlot(handles);
guidata(hObject,handles)
end

function delok_Callback(hObject, eventdata)
handles = guidata(hObject);
linestr = get(handles.hStatic.hdelline, 'string');
lineno = str2double(linestr);
if (isnan(lineno) || isempty(lineno))
   msgbox('Enter the line number');
elseif (lineno > length(handles.EXPDATA)) || (lineno < 0)
   msgbox('Not a valid line number');
else
   handles=delete_line(lineno,handles);
end
guidata(hObject,handles)
end

%------------------------------------------------------
% delete line
function handles = delete_line(del_line,handles)
lines = handles.var.lines;
for i=1:lines
   delete(handles.hLine(i).hlinenr);
   delete(handles.hLine(i).hMaterial);
   delete(handles.hLine(i).hType);
   delete(handles.hLine(i).hRate);
   delete(handles.hLine(i).hGas);
   delete(handles.hLine(i).hspace1);
   delete(handles.hLine(i).hspace2);
   if (ishandle(handles.hLine(i).gbg1))
      delete(handles.hLine(i).gbg1);
      delete(handles.hLine(i).gbg2);
   end
   delete(handles.hLine(i).hl1);
   delete(handles.hLine(i).hr1);
   delete(handles.hLine(i).hl2);
   delete(handles.hLine(i).hr2);
end
%
handles.EXPDATA = handles.EXPDATA([1:del_line-1 del_line+1:lines]);
lines=lines - 1;
handles.var.lines = lines;
for i=1:lines
   handles = update_line(handles,i);
end
end

    
% Initialization
function handles = InitializePyroPlot(handles)
%
dummyh = -1;
%variables
handles.var.yLabel1Left = ' ';
handles.var.yLabel1Right = ' ';
handles.var.yLabel2Left =' ';
handles.var.yLabel2Right = ' ';
handles.var.lines = 0;
handles.var.starting_path = ' ';
handles.var.save_name = 'TGA Material 10 Air';
handles.var.save_material = 'Material';
handles.var.save_rate = '10';
handles.var.save_path = ' ';
handles.hAxes.ax11=dummyh;
handles.hAxes.ax12=dummyh;
handles.hAxes.ax21=dummyh;
handles.hAxes.ax22=dummyh;
handles.Options.Fig_xtype=[0 0]; % initialized x-axis as temperature
handles.Options.rb(1,1:2)=0; % initialized TGA as mass-%
handles.Options.NFilter = [1 1 1 1];
%initialize vectors
handles.EXPDATA(1).temperature=0;
handles.EXPDATA(1).time = 0;
handles.EXPDATA(1).TGA = 0;
handles.EXPDATA(1).gradient = 0;
handles.EXPDATA(1).opgradient = 0; % opposite gradient, -1*gradient
handles.EXPDATA(1).DSC = 0;
handles.EXPDATA(1).ConeHRR = 0;
handles.EXPDATA(1).ConeMLR = 0;
handles.EXPDATA(1).ConeEHC=0;
handles.EXPDATA(1).type = '';
handles.EXPDATA(1).rate = 1;
handles.EXPDATA(1).material = '';
handles.EXPDATA(1).gas = '';
handles.EXPDATA(1).name = '';
handles.EXPDATA(1).path = '';
handles.EXPDATA(1).check = [0 0 0 0]; %check boxes not cheked
handles.EXPDATA(1).pair = 0; %index to TGA-DSC pairs, if not pair, this is 0
handles.EXPDATA(1).EXO = 0;
handles.EXPDATA(1).p_set_l = [1, 2]; %solid line, width 2
handles.EXPDATA(1).p_set_l = [2, 2]; %dash line, width 2
handles.A = 100E-4; %surface area of cone sample
handles.hCbvector(1,:)=[0 0 0 0]; % vector for checkbox callbacks
handles.hRBvector(1,:)=[1 1 1 1 1 1]; % vector for radiobutton callbacks
%fds input
handles.template_path = '';
handles.FdsExe_path = '';
%DSC structure
handles.DSC.analysis = 1;
handles.var.checkvector(1)=0;
%TGA stucture
handles.TGA.analysis = 1;
%line handles
handles.hLine(1).hlinenr   = dummyh;
handles.hLine(1).hMaterial = dummyh; 
handles.hLine(1).hType     = dummyh;
handles.hLine(1).hRate     = dummyh;
handles.hLine(1).hGas      = dummyh;
handles.hLine(1).hspace1   = dummyh;
handles.hLine(1).hspace2   = dummyh;
handles.hLine(1).gbg1      = dummyh;
handles.hLine(1).gbg2      = dummyh;
handles.hLine(1).hl1       = dummyh;
handles.hLine(1).hr1       = dummyh;
handles.hLine(1).hl2       = dummyh;
handles.hLine(1).hr2       = dummyh;
handles.hFigure            = [dummyh dummyh];
end % Initialization done

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other functions

    
%
%add new line every time user adds a data serie
function handles = update_line(handles,n_line)
%
% co-ordinates
ypos = 0.92-n_line*0.0357;
ypos2 = 0.925-n_line*0.0357;
ydim = 0.0286;
%
hlinenr = uicontrol(handles.hPyroPlot, 'Style', 'text', ...
   'String', {mat2str(n_line)}, ...
   'HorizontalAlignment','left',...
   'Units', 'normalized', ... 
   'Position', [0.003 0.92-n_line*0.0357 0.03 ydim]);
hMaterial = uicontrol(handles.hPyroPlot,'Style','text',...
   'String',{handles.EXPDATA(n_line).material},...
   'Units', 'normalized', ...
   'Position',[0.02 ypos 0.07 ydim]);
if (~strcmp(handles.EXPDATA(n_line).type, 'ConeFDS') && ~strcmp(handles.EXPDATA(n_line).type, 'TGAFDS'))
    hType = uicontrol(handles.hPyroPlot,'Style','text',...
   'String',{handles.EXPDATA(n_line).type},...
   'Units', 'normalized', ...
   'Position',[0.10 ypos 0.05 ydim]);
elseif strcmp(handles.EXPDATA(n_line).type, 'ConeFDS')
    hType = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','Cone',...
   'Units', 'normalized', ...
   'Position',[0.10 ypos 0.05 ydim]);
elseif strcmp(handles.EXPDATA(n_line).type, 'TGAFDS')
    hType = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','TGA',...
   'Units', 'normalized', ...
   'Position',[0.10 ypos 0.05 ydim]);
end
hRate = uicontrol(handles.hPyroPlot,'Style','text',...
   'String',{handles.EXPDATA(n_line).rate},...
   'Units', 'normalized', ...
   'Position',[0.145 ypos 0.05 ydim]);
hGas = uicontrol(handles.hPyroPlot,'Style','text',...
   'String',{handles.EXPDATA(n_line).gas},...
   'Units', 'normalized', ...
   'Position',[0.206 ypos 0.05 ydim]);
hspace1 = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','|',...
   'Units', 'normalized', ...
   'Position',[0.265 ypos 0.007 ydim]);
hspace2 = uicontrol(handles.hPyroPlot,'Style','text',...
   'String','|',...
   'Units', 'normalized', ...
   'Position',[0.48 ypos 0.007 ydim]);

% radio buttons
radio_set = 0;
if (strcmp(handles.EXPDATA(n_line).type,'TGA') || strcmp(handles.EXPDATA(n_line).type,'TGAFDS'))
   radio_set = 1;
   %button group 1
   gbg1 = uibuttongroup('Parent',handles.hPyroPlot,...
      'Position',[0.34 ypos 0.135 0.033]);
   u11 = uicontrol('Style', 'Radio', 'String', '%',...
      'Units', 'normalized', ...
      'pos',[0.007 0.04 0.3 0.87], 'parent', gbg1, 'HandleVisibility', 'off', 'Tag', 'r11');
   u12 = uicontrol('Style', 'Radio', 'String', '%/s',...
      'Units', 'normalized', ...
      'pos',[0.3 0.04 0.3 0.87], 'parent', gbg1, 'HandleVisibility', 'off', 'Tag', 'r12');
   u13 = uicontrol('Style', 'Radio', 'String', '-%/s',...
      'Units', 'normalized', ...
      'pos',[0.63 0.04 0.33 0.87], 'parent', gbg1, 'HandleVisibility', 'off', 'Tag', 'r12');
   %button group 2
   gbg2 = uibuttongroup('Parent',handles.hPyroPlot,...
      'Position',[0.56 ypos 0.135 0.033]);
   u21 = uicontrol('Style', 'Radio', 'String', '%',...
      'Units', 'normalized', ...
      'pos',[0.007 0.04 0.3 0.87], 'parent', gbg2, 'HandleVisibility', 'off', 'Tag', 'r11');
   u22 = uicontrol('Style', 'Radio', 'String', '%/s',...
      'Units', 'normalized', ...
      'pos',[0.3 0.04 0.3 0.87], 'parent', gbg2, 'HandleVisibility', 'off', 'Tag', 'r12');
   u23 = uicontrol('Style', 'Radio', 'String', '-%/s',...
      'Units', 'normalized', ...
      'pos',[0.63 0.04 0.33 0.87], 'parent', gbg2, 'HandleVisibility', 'off', 'Tag', 'r12');
elseif (strcmp(handles.EXPDATA(n_line).type,'Cone') || strcmp(handles.EXPDATA(n_line).type,'ConeFDS'))
   radio_set = 1;
   %button group 1
   gbg1 = uibuttongroup('Parent',handles.hPyroPlot,...
      'Position',[0.34 ypos 0.135 0.033]);
   u11 = uicontrol('Style', 'Radio', 'String', 'HRR',...
      'Units', 'normalized', ...
      'pos',[0.007 0.04 0.3 0.87], 'parent', gbg1, 'HandleVisibility', 'off', 'Tag', 'r11');
   u12 = uicontrol('Style', 'Radio', 'String', 'MLR',...
      'Units', 'normalized', ...
      'pos',[0.3 0.04 0.3 0.87], 'parent', gbg1, 'HandleVisibility', 'off', 'Tag', 'r12');
   u13 = uicontrol('Style', 'Radio', 'String', 'EHC',...
      'Units', 'normalized', ...
      'pos',[0.63 0.04 0.33 0.87], 'parent', gbg1, 'HandleVisibility', 'off', 'Tag', 'r12');
   % Initialize some button group properties.
   %button group 2
   gbg2 = uibuttongroup('Parent',handles.hPyroPlot,...
      'Position',[0.56 ypos 0.135 0.033]);
   u21 = uicontrol('Style', 'Radio', 'String', 'HRR',...
      'Units', 'normalized', ...
      'pos',[0.007 0.04 0.3 0.87], 'parent', gbg2, 'HandleVisibility', 'off', 'Tag', 'r11');
   u22 = uicontrol('Style', 'Radio', 'String', 'MLR',...
      'Units', 'normalized', ...
      'pos',[0.3 0.04 0.3 0.87], 'parent', gbg2, 'HandleVisibility', 'off', 'Tag', 'r12');
   u23 = uicontrol('Style', 'Radio', 'String', 'EHC',...
      'Units', 'normalized', ...
      'pos',[0.63 0.04 0.33 0.87], 'parent', gbg2, 'HandleVisibility', 'off', 'Tag', 'r12');
end
% Initialize some button group properties.
if radio_set
   set(gbg1,'SelectionChangeFcn',{@selcbk,1});
   set(gbg1,'Visible','on');
   set(gbg2,'SelectionChangeFcn',{@selcbk,2});
   set(gbg2,'Visible','on');
   handles.hRBvector(n_line, :) = [u11 u12 u13 u21 u22 u23];
else
   gbg1 = -1;
   gbg2 = -1;
end

%check boxes

hl1 = uicontrol(handles.hPyroPlot,'Style','checkbox',...
   'CallBack', {@h_cb,1},...
   'Value',handles.EXPDATA(n_line).check(1),...
   'Units', 'normalized', ...
   'Position',[0.28 ypos2 0.02 ydim]);
hr1 = uicontrol(handles.hPyroPlot,'Style','checkbox',...
   'CallBack', {@h_cb,2},...
   'Value',handles.EXPDATA(n_line).check(2),...
   'Units', 'normalized', ...
   'Position',[0.31 ypos2 0.02 ydim]);
hl2 = uicontrol(handles.hPyroPlot,'Style','checkbox',...
   'CallBack', {@h_cb,3},...
   'Value',handles.EXPDATA(n_line).check(3),...
   'Units', 'normalized', ...
   'Position',[0.495 ypos2 0.02 ydim]);
hr2 = uicontrol(handles.hPyroPlot,'Style','checkbox',...
   'CallBack', {@h_cb,4},...
   'Value',handles.EXPDATA(n_line).check(4),...
   'Units', 'normalized', ...
   'Position',[0.525 ypos2 0.02 ydim]);

handles.hCBvector(n_line,:)=[hl1 hr1 hl2 hr2];

handles.hLine(n_line).hlinenr = hlinenr;
handles.hLine(n_line).hMaterial = hMaterial;
handles.hLine(n_line).hType = hType;
handles.hLine(n_line).hRate = hRate;
handles.hLine(n_line).hGas = hGas;
handles.hLine(n_line).hspace1 = hspace1;
handles.hLine(n_line).hspace2 = hspace2;
handles.hLine(n_line).gbg1 = gbg1;
handles.hLine(n_line).gbg2 = gbg2;
handles.hLine(n_line).hl1 = hl1;
handles.hLine(n_line).hr1 = hr1;
handles.hLine(n_line).hl2 = hl2;
handles.hLine(n_line).hr2 = hr2;
       
end

%------------------------------------------------------
% setdata
function [x,y]=setdata(l,f,c,handles)
% l is line number
% f is figure (1 or 2)
% c is checkIndex (1,2,3 or 4)
switch handles.Options.Fig_xtype(f)
   case 0
      x=handles.EXPDATA(l).temperature;
   case 1
      x=handles.EXPDATA(l).time;
end
switch handles.EXPDATA(l).type
   case 'TGA'
      if handles.Options.rb(l,f)==0
         yy=handles.EXPDATA(l).TGA*100;
      elseif handles.Options.rb(l,f)==1
         yy=handles.EXPDATA(l).gradient*100;
      elseif handles.Options.rb(l,f)==2
         yy=handles.EXPDATA(l).opgradient*100;
      end
    case 'TGAFDS'
      if handles.Options.rb(l,f)==0
         yy=handles.EXPDATA(l).TGA*100;
      elseif handles.Options.rb(l,f)==1
         yy=handles.EXPDATA(l).gradient*100;
      elseif handles.Options.rb(l,f)==2
         yy=handles.EXPDATA(l).opgradient*100;
      end
   case 'DSC'
      yy=handles.EXPDATA(l).DSC/1000;
   case 'Cone'
      if (handles.Options.rb(l,f)==0)
         yy=handles.EXPDATA(l).ConeHRR;
      elseif (handles.Options.rb(l,f)==1)
         yy=handles.EXPDATA(l).ConeMLR;
      elseif (handles.Options.rb(l,f)==2)
         yy=handles.EXPDATA(l).ConeEHC;
      end
   case 'ConeFDS'
      if (handles.Options.rb(l,f)==0)
         yy=handles.EXPDATA(l).ConeHRR;
      elseif (handles.Options.rb(l,f)==1)
         yy=handles.EXPDATA(l).ConeMLR;
      elseif (handles.Options.rb(l,f)==2)
         yy=handles.EXPDATA(l).ConeEHC;
      end
end
if (handles.Options.NFilter(c)>1)
   y=filtNs(yy,handles.Options.NFilter(c));
else
   y = yy;
end
end

%----------------------------------------------
% updatePlot
% function plots chosen series again everytime something is called
%

function plot_data_Callback(hObject, eventdata)
handles = guidata(hObject);
%handles = updata_plot(handles);
guidata(hObject, handles);

end

function handles = updatePlot(handles)
%first collect all data to strings

%look trought all the check boxes and collect information

dataset1l = ' ';
dataset1r = ' ';
dataset2l = ' ';
dataset2r = ' ';
legend1l = ' ';
legend1r = ' ';
legend1=' ';
legend2 = ' ';
legend2l = ' ';
legend2r = ' ';
limits1=[1 0]; %min and max values of x-axis (figure1)
limits2=[1 0]; %figure2
h11 = [];
h12 = [];
h21 = [];
h22 = [];
ph_11 = [];
ph_12 = [];
ph_21 = [];
ph_22 = [];
% first check if Cone data is plotted
for i=1:handles.var.lines
   if (strcmp(handles.EXPDATA(i).type,'Cone') || strcmp(handles.EXPDATA(i).type,'ConeFDS'))
      if any(handles.EXPDATA(i).check(1:2))
         handles.Options.Fig_xtype(1)= 1;
      end
      if any(handles.EXPDATA(i).check(3:4))
         handles.Options.Fig_xtype(2) = 1;
      end      
   end
end
for i=1:handles.var.lines
   name=' ';
   nameindex=1;
   for k=1:length(handles.EXPDATA(i).name)
      if strcmp(handles.EXPDATA(i).name(k),'_')
         name(nameindex)='\';
         name(nameindex+1)='_';
         nameindex=nameindex+2;
      else
         name(nameindex) = handles.EXPDATA(i).name(k);
         nameindex=nameindex+1;
      end
   end
   %figure 1 left
   if (handles.EXPDATA(i).check(1) == 1)
      k = length(ph_11)+1;
      ph_11(k) = i;
      [x,y]=setdata(i,1,1,handles);
      if (x==0), msgbox('Problem with x-data'), return;end
      limits1(1)=min(limits1(1),min(x));
      limits1(2)=max(limits1(2),max(x));
      xx=mat2str(x);
      yy=mat2str(y);
      if length(dataset1l)==1
         dataset1l = [dataset1l xx, ',' yy];
      else
         dataset1l = [dataset1l ',' xx ',' yy];
      end
      if length(legend1l) ==1
         legend1l = [legend1l  '''' name ''''];
      else
         legend1l = [legend1l ',' '''' name ''''];
      end
   end
            
   %figure 1 right
   if (handles.EXPDATA(i).check(2) == 1)
       k = length(ph_12)+1;
      ph_12(k) = i;
      [x,y]=setdata(i,1,2,handles);
      if (x==0), msgbox('Problem with x-data'), return;end
      limits1(1)=min(limits1(1),min(x));
      limits1(2)=max(limits1(2),max(x));
      xx=mat2str(x);
      yy=mat2str(y);
      if length(dataset1r)==1
         dataset1r = [dataset1r xx, ',' yy];
      else
         dataset1r = [dataset1r ',' xx ',' yy];
      end
      if length(legend1r) ==1
         legend1r = [legend1r  '''' name ''''];
      else
         legend1r = [legend1r ',' '''' name ''''];
      end
   end
            
   %figure 2 left
   if (handles.EXPDATA(i).check(3) == 1)
       k = length(ph_21)+1;
      ph_21(k) = i;
      [x,y]=setdata(i,2,3,handles);
      if (x==0), msgbox('Problem with x-data'), return;end
      limits2(1)=min(limits2(1),min(x));
      limits2(2)=max(limits2(2),max(x));
      xx=mat2str(x);
      yy=mat2str(y);
      if length(dataset2l)==1
         dataset2l = [dataset2l xx ',' yy];
      else
         dataset2l = [dataset2l ',' xx ',' yy];
      end
      if length(legend2l) ==1
         legend2l = [legend2l  '''' name ''''];
      else
         legend2l = [legend2l ',' '''' name ''''];
      end
   end
            
   %figure 2 right
   if handles.EXPDATA(i).check(4) == 1
       k = length(ph_22)+1;
      ph_22(k) = i;
      [x,y]=setdata(i,2,4,handles);
      if (x==0), msgbox('Problem with x-data'), return;end
      limits2(1)=min(limits2(1),min(x));
      limits2(2)=max(limits2(2),max(x));
      xx=mat2str(x);
      yy=mat2str(y);
      if length(dataset2r)==1
         dataset2r = [dataset2r xx ',' yy];
      else
         dataset2r = [dataset2r ',' xx ',' yy];
      end
      if length(legend2r) ==1
         legend2r = [legend2r  '''' name ''''];
      else
         legend2r = [legend2r ',' '''' name ''''];
      end
   end
end % for 1 to lines
%
% figure1
if (length(dataset1l) ~=1 || length(dataset1r)~=1)
   if length(legend1r) == 1
      legend1 = legend1l;
   elseif length(legend1l) == 1
      legend1 = legend1r;
   elseif (length(legend1l) ~= 1) && (length(legend1r) ~=1)
      legend1 = [legend1l ',' legend1r];
   end
   if (isfield(handles,'hFigure'))
      if ishandle(handles.hFigure(1))   
          if isfield(handles,'fig_index_11')
            for i = 1:length(handles.fig_index_11(:,1))
                index = handles.fig_index_11(i,2);         
                line_handle = handles.fig_index_11(i,1);  
                
                handles.EXPDATA(index).color_l = get(line_handle, 'Color');
                style = 1;
                if strcmp(get(line_handle,'LineStyle'), '--')
                    style = 2;
                elseif strcmp(get(line_handle,'LineStyle'), ':')
                    style = 3;
                elseif strcmp(get(line_handle,'LineStyle'), '-.')
                    style = 4;
                end
                handles.EXPDATA(index).p_set_l = [style, get(line_handle,'LineWidth')];
            end
            handles.fig1_fontSize_left = get(handles.hAxes.ax11, 'FontSize');
          
          end
          
          if isfield(handles,'fig_index_12')
            for i = 1:length(handles.fig_index_12(:,1))
                index = handles.fig_index_12(i,2);         
                line_handle = handles.fig_index_12(i,1);  
                
                handles.EXPDATA(index).color_r = get(line_handle, 'Color');
                style = 1;
                if strcmp(get(line_handle,'LineStyle'), '--')
                    style = 2;
                elseif strcmp(get(line_handle,'LineStyle'), ':')
                    style = 3;
                elseif strcmp(get(line_handle,'LineStyle'), '-.')
                    style = 4;
                end
                handles.EXPDATA(index).p_set_r = [style, get(line_handle,'LineWidth')];
            end
            handles.fig1_fontSize_right = get(handles.hAxes.ax12, 'FontSize');
          end
          
         handles.fig1_position = get(handles.hFigure(1), 'Position');
         close(handles.hFigure(1));
      else
          load plot_settings.mat
          handles.fig1_fontSize_left = settings.left.FontSize;
          handles.fig1_fontSize_right = settings.right.FontSize;
          handles.fig1_position = settings.Position;
          clear settings
      end
   end
   handles.hFigure(1) = figure(...
      'Visible', 'off',...
      'Name', 'FIGURE 1',...
      'HandleVisibility','callback', ...
      'NumberTitle', 'off', ...
      'Units', 'normalized', ...
      'Position',handles.fig1_position,...
      'Color', get(0,...
      'defaultuicontrolbackgroundcolor'));
  
   set(0,'CurrentFigure',handles.hFigure(1));
   set(handles.hFigure(1), 'Visible', 'on');
   
   s1 = ['plot(' dataset1l ')'];
   if (length(dataset1l)==1)
      t=['legend(h12, ' legend1 ')'];
   elseif (length(dataset1r) == 1)
      t=['legend(h11, ' legend1 ')'];
   elseif (length(dataset1l) ~=1 && length(dataset1r)~=1)
      t = ['legend([h11;h12], ' legend1 ')'];
   end
   handles.hAxes.ax11=gca;
   % plot figure with eval
   if(length(dataset1l) ~=1)
      h11=eval(s1);
      handles.fig_index_11 = [h11, ph_11']; 
      ylabel(handles.hAxes.ax11, handles.var.yLabel1Left);
     for i=1:length(handles.fig_index_11(:,1))
            index = handles.fig_index_11(i,2);
            line_handle = handles.fig_index_11(i,1);
        if isfield(handles.EXPDATA(index), 'color_l') % if fig 1 exists
            set(line_handle, 'Color', handles.EXPDATA(index).color_l); 
        else
            %handles.EXPDATA(index).color_l = get(line_handle, 'Color');
        end
            set(line_handle, 'LineWidth', handles.EXPDATA(index).p_set_l(2));
            
            if isequal(handles.EXPDATA(index).p_set_l(1),1)
                style = '-';
            elseif isequal(handles.EXPDATA(index).p_set_l(1),2)
                style = '--';
            elseif isequal(handles.EXPDATA(index).p_set_l(1),3)
                style = ':';
            else
               style = '-.';
            end
            set(line_handle, 'LineStyle', style);
     end
   end
   set(handles.hAxes.ax11, 'FontSize', handles.fig1_fontSize_left)
   set(handles.hAxes.ax11, 'XLim', limits1);
   hold on;
   handles.hAxes.ax12=axes('Position', get(handles.hAxes.ax11, 'Position'), ...
      'YAxisLocation', 'right', ...
      'XLim', limits1, ...
      'Color', 'none', ...
      'Visible', 'off');
   s2 = ['plot('  dataset1r ', ''Parent'',handles.hAxes.ax12)'];
   hold on;
   if(length(dataset1r) ~=1)
      h12=eval(s2);
      handles.fig_index_12 = [h12, ph_12']; 
      ylabel(handles.hAxes.ax12, handles.var.yLabel1Right);
     for i=1:length(handles.fig_index_12(:,1))
            index = handles.fig_index_12(i,2);
            line_handle = handles.fig_index_12(i,1);
        if isfield(handles.EXPDATA(index), 'color_r') % if fig 1 exists
            set(line_handle, 'Color', handles.EXPDATA(index).color_r); 
        else
            %handles.EXPDATA(index).color_l = get(line_handle, 'Color');
        end
            set(line_handle, 'LineWidth', handles.EXPDATA(index).p_set_r(2));
            
            if isequal(handles.EXPDATA(index).p_set_r(1),1)
                style = '-';
            elseif isequal(handles.EXPDATA(index).p_set_r(1),2)
                style = '--';
            elseif isequal(handles.EXPDATA(index).p_set_r(1),3)
                style = ':';
            else
               style = '-.';
            end
            set(line_handle, 'LineStyle', style);
     end
   
   set(handles.hAxes.ax12, 'FontSize', handles.fig1_fontSize_right)
      colormap(hsv);
      
      ylabel(handles.hAxes.ax12, handles.var.yLabel1Right);
      set(handles.hAxes.ax12, 'Position', get(handles.hAxes.ax11, 'Position'));
      set(handles.hAxes.ax12, 'Visible', 'on');
   end
   handles.legends1 = eval(t);
   %set(handles.legends1, 'Location', 'BestOutside');
   set(handles.hAxes.ax12, 'Position', get(handles.hAxes.ax11, 'Position'));
   switch handles.Options.Fig_xtype(1)
      case 0
         xlabel(handles.hAxes.ax11,'Temperature (\circC)');
      case 1
         xlabel(handles.hAxes.ax11,'Time (s)');
   end
   guidata(handles.hFigure(1), handles);
   
   %save handle information
plot_handles.fig1.figure = handles.hFigure(1);
plot_handles.fig1.left.ax = handles.hAxes.ax11;
plot_handles.fig1.right.ax = handles.hAxes.ax11;
plot_handles.fig1.left.lines = h11;        
plot_handles.fig1.right.lines = h12;   

save plot_handles.mat plot_handles
clear plot_handles

end

%figure2
if ((length(dataset2l) ~=1) || (length(dataset2r)~=1))
   if length(legend2r) == 1
      legend2 = legend2l;
   elseif length(legend2l) == 1
      legend2 = legend2r;
   elseif (length(legend2l) ~= 1) && (length(legend2r) ~=1)
      legend2 = [legend2l ',' legend2r];
   end
   if (isfield(handles,'hFigure'))
       if ishandle(handles.hFigure(2))   
          if isfield(handles,'fig_index_21')
            for i = 1:length(handles.fig_index_21(:,1))
                index = handles.fig_index_21(i,2);         
                line_handle = handles.fig_index_21(i,1);  
                
                handles.EXPDATA(index).color_l = get(line_handle, 'Color');
                style = 1;
                if strcmp(get(line_handle,'LineStyle'), '--')
                    style = 2;
                elseif strcmp(get(line_handle,'LineStyle'), ':')
                    style = 3;
                elseif strcmp(get(line_handle,'LineStyle'), '-.')
                    style = 4;
                end
                handles.EXPDATA(index).p_set_l = [style, get(line_handle,'LineWidth')];
            end
            handles.fig1_fontSize_left = get(handles.hAxes.ax21, 'FontSize');
          
          end
          
          if isfield(handles,'fig_index_22')
            for i = 1:length(handles.fig_index_22(:,1))
                index = handles.fig_index_22(i,2);         
                line_handle = handles.fig_index_22(i,1);  
                
                handles.EXPDATA(index).color_r = get(line_handle, 'Color');
                style = 1;
                if strcmp(get(line_handle,'LineStyle'), '--')
                    style = 2;
                elseif strcmp(get(line_handle,'LineStyle'), ':')
                    style = 3;
                elseif strcmp(get(line_handle,'LineStyle'), '-.')
                    style = 4;
                end
                handles.EXPDATA(index).p_set_r = [style, get(line_handle,'LineWidth')];
            end
            handles.fig1_fontSize_right = get(handles.hAxes.ax22, 'FontSize');
          end
          
         handles.fig1_position = get(handles.hFigure(2), 'Position');
         
         close(handles.hFigure(2));
         else
          load plot_settings.mat
          handles.fig1_fontSize_left = settings.left.FontSize;
          handles.fig1_fontSize_right = settings.right.FontSize;
          handles.fig1_position = settings.Position;
          clear settings
       end
    end
   %figure2
   handles.hFigure(2) = figure(...
      'Visible', 'off',...
      'Name', 'FIGURE 2',...
      'NumberTitle', 'off', ...
      'HandleVisibility','callback', ...
      'Units', 'normalized', ...
      'Position',handles.fig1_position,...
      'Color', get(0,...
      'defaultuicontrolbackgroundcolor'));
%  handles.hSavepng2 = uicontrol(handles.hFigure(2), ...
%      'Style', 'pushbutton', ...
%      'String', 'Smaller', ...
%      'Units', 'Normalized', ...
%      'Position', [0.85 0.2 0.1 0.05], ...
%      'Callback', @savepng2_cb);
   set(0,'CurrentFigure',handles.hFigure(2));
   set(handles.hFigure(2), 'Visible', 'on');
   s1 = ['plot(' dataset2l ')'];
   if length(dataset2l)==1
      t=['legend(h22, ' legend2 ')'];
   elseif length(dataset2r) == 1
      t=['legend(h21, ' legend2 ')'];
   elseif length(dataset2l) ~=1 && length(dataset2r)~=1
      t = ['legend([h21;h22], ' legend2 ')'];
   end
   handles.hAxes.ax21=gca;
   % plot figure with eval
   if(length(dataset2l) ~=1)
      h21=eval(s1);
      handles.fig_index_21 = [h21, ph_21']; 
      ylabel(handles.hAxes.ax21, handles.var.yLabel2Left);
      for i=1:length(handles.fig_index_21(:,1))
            index = handles.fig_index_21(i,2);
            line_handle = handles.fig_index_21(i,1);
        if isfield(handles.EXPDATA(index), 'color_l') % if fig 1 exists
            set(line_handle, 'Color', handles.EXPDATA(index).color_l); 
        else
            %handles.EXPDATA(index).color_l = get(line_handle, 'Color');
        end
            set(line_handle, 'LineWidth', handles.EXPDATA(index).p_set_l(2));
            
            if isequal(handles.EXPDATA(index).p_set_l(1),1)
                style = '-';
            elseif isequal(handles.EXPDATA(index).p_set_l(1),2)
                style = '--';
            elseif isequal(handles.EXPDATA(index).p_set_l(1),3)
                style = ':';
            else
               style = '-.';
            end
            set(line_handle, 'LineStyle', style);
     end
   end
   set(handles.hAxes.ax21, 'FontSize', handles.fig1_fontSize_left)
   set(handles.hAxes.ax21, 'XLim', limits2);
   hold on;
   handles.hAxes.ax22=axes('Position', get(handles.hAxes.ax21, 'Position'), ...
      'YAxisLocation', 'right', ...
      'XLim', limits2, ...
      'Color', 'none', ...
      'Visible', 'off');
   s2 = ['plot('  dataset2r ',''Parent'',handles.hAxes.ax22)'];
   hold on;
   if(length(dataset2r) ~=1)
      h22=eval(s2);
      handles.fig_index_22 = [h22, ph_22']; 
      ylabel(handles.hAxes.ax22, handles.var.yLabel2Right);
     for i=1:length(handles.fig_index_22(:,1))
            index = handles.fig_index_22(i,2);
            line_handle = handles.fig_index_22(i,1);
        if isfield(handles.EXPDATA(index), 'color_r') % if fig 1 exists
            set(line_handle, 'Color', handles.EXPDATA(index).color_r); 
        else
            %handles.EXPDATA(index).color_l = get(line_handle, 'Color');
        end
            set(line_handle, 'LineWidth', handles.EXPDATA(index).p_set_r(2));
            
            if isequal(handles.EXPDATA(index).p_set_r(1),1)
                style = '-';
            elseif isequal(handles.EXPDATA(index).p_set_r(1),2)
                style = '--';
            elseif isequal(handles.EXPDATA(index).p_set_r(1),3)
                style = ':';
            else
               style = '-.';
            end
            set(line_handle, 'LineStyle', style);
     end
   
   set(handles.hAxes.ax22, 'FontSize', handles.fig1_fontSize_right)
   
      ylabel(handles.hAxes.ax22, handles.var.yLabel2Right);
      set(handles.hAxes.ax22, 'Position', get(handles.hAxes.ax21, 'Position'));
      set(handles.hAxes.ax22, 'Visible', 'on');
   end
   handles.legends2 = eval(t);
   
   set(handles.hAxes.ax22, 'Position', get(handles.hAxes.ax21, 'Position'));
   switch handles.Options.Fig_xtype(2)
      case 0
         xlabel(handles.hAxes.ax21,'Temperature');
      case 1
         xlabel(handles.hAxes.ax21,'Time');
   end
   guidata(handles.hFigure(2), handles);
   %save handle information
    plot_handles.fig2.figure = handles.hFigure(2);
    plot_handles.fig2.left.ax = handles.hAxes.ax21;
    plot_handles.fig2.right.ax = handles.hAxes.ax21;
    plot_handles.fig2.left.lines = h21;        
    plot_handles.fig2.right.lines = h22;   

save plot_handles.mat plot_handles
clear plot_handles
end

guidata(handles.hPyroPlot, handles);
set(handles.hPyroPlot, 'Visible', 'on');
end