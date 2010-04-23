function handles = ppread(handles,DType)
% ppREAD       PyroPlot data reading function
%
% Initialize handles
handles.heditfig = -1;
handles.hname = -1;
handles.hgaspopup = -1;
handles.hmaterial = -1;
handles.hrate = -1;
handles.hTGAUpopup = -1;
handles.hDSCUpopup = -1;
handles.hA_r = -1;
%read and save the data of selected file
dstr = DType;
if strcmp(DType,'DTA'),dstr='TGA';end
file_name = '';
dir_path = handles.var.starting_path;
if strcmp(DType, 'TGAFDS')
rdstr = ['Read TGA from FDS output'];
[file_name, dir_path]=uigetfile('*_prof_01.csv', rdstr, handles.var.starting_path);
elseif strcmp(DType, 'Cone')
    [file_name, dir_path]=uigetfile('*_red.csv', 'Read Cone', handles.var.starting_path);
elseif strcmp(DType, 'ConeFDS')
    [file_name, dir_path]=uigetfile('*_hrr.csv', 'Read Cone from output', handles.var.starting_path);
else
rdstr = ['Read ' dstr];
[file_name, dir_path]=uigetfile('*.txt', rdstr, handles.var.starting_path);
end
%
if (~isequal(file_name,0))
   % overwrite the first initialized vector
   file=fullfile(dir_path, file_name);
   handles.var.starting_path = dir_path;
   lines = handles.var.lines;
   values = [];
   time = [];
   if strcmp(DType, 'TGAFDS')
    M = readdata(file,3);
    M = M(:,[1 4 5]);
    Mass = [M(:,1),M(:,2).*M(:,3)];
    tempfile = [file_name(1:length(file_name)-12) '_devc.csv'];
    tempfile = fullfile(dir_path, tempfile);
    Temp = readdata(tempfile, 3); %2?
    
    %if time step is different in massfile and tempfile
    l=1;
     if min(Temp(:,1))>min(Mass(:,1))  %if first value smaller        
            while min(Temp(:,1))>Mass(l,1)
                l=l+1;
            end
     end
     
    Temp = interp1(Temp(:,1), Temp(:,2), Mass(l:length(Mass(:,1)),1));
    
    Temp = removeNaNM(Temp,3);
    %%%
    N = length(Mass);
    if N > length(Temp)
    Mass = Mass(1:length(Temp),:);
    elseif N < length(Temp)
    Temp = Temp(1:N);
    end
    values(:,1)=Temp;
    values(:,2)=Mass(:,2);
    time = Mass(:,1);
   %convert to mass-fraction
   values(:,2)=values(:,2)./max(values(:,2));
   elseif strcmp(DType, 'Cone')
       
       s = '%s';
   for i =1:23
       s = [s ' %f'];
   end
   title=textread(file, '%s', 'delimiter', ','); %'\n'
   titlelist=char(title);
   datalines = length(titlelist);
   fid=fopen(file);
   
   data=textscan(fid, s, datalines, 'headerlines', 2, 'delimiter', ',');
   fclose(fid);
%   data = readdata(file, 2);
   values(:,1) = cell2mat(data(:,3));
   values(:,2) = cell2mat(data(:,4)); %HRR (kw/m^2)
   values(:,3) = cell2mat(data(:,8)).*10^(-3); %MLR (in file unit is g/sm^2, now kg/sm^2)
   values(:,4) = cell2mat(data(:,5)); %EHC unit MJ/kg
   elseif strcmp(DType, 'ConeFDS')
   data1 = readdata(file, 2);
   devc_file = file(1:(length(file)-7));
   devc_file = [devc_file, 'devc.csv'];
   data2 = readdata(devc_file, 2);
   values(:,1) = data1(:,1);
   values(:,2) = data1(:,2)./handles.A; %HRR
   values(:,3) = interp1(data2(:,1),data2(:,5),data1(:,1)); %MLR (surface density)
   values(:,3) = -gradient(values(:,3),data1(:,1));
   warning off all
   values(:,4) = values(:,2)./values(:,3).*10^(-3); %EHC=HRR/MLR (MJ/kg)
   warning on
   else
   [headerlines, datalines,EXO] = getHeaderlines(file);
   fid=fopen(file);
   data=textscan(fid, '%f %f', datalines, 'headerlines', headerlines, 'delimiter', ';');
   fclose(fid);
   values = cell2mat(data);
   end
   if (length(values(:,1))<=1)
      msgbox('Too few data lines in the file')
      return
   end
   if strcmp(DType,'DTA')
      rdstr = ['Read DSC'];
      [file_name, dir_path]=uigetfile('*.txt', rdstr, handles.var.starting_path);
      file=fullfile(dir_path, file_name);
      [headerlines, datalines,EXO] = getHeaderlines(file);
      fid=fopen(file);
      data=textscan(fid, '%f %f', datalines, 'headerlines', headerlines, 'delimiter', ';');
      fclose(fid);
      values2 = cell2mat(data);
      if (length(values2(:,1))<=1)
         msgbox('Too few data lines in the file')
         return
      end
   end
   lines = lines + 1;
   handles.var.lines = lines;
   handles.EXPDATA(lines).gas = 'Air';
   handles.EXPDATA(lines).rate = 10;
   switch DType
      case('TGA')
         handles.EXPDATA(lines).temperature = values(:,1);
         handles.EXPDATA(lines).TGA = values(:,2);
         handles.EXPDATA(lines).DSC = 0;
         handles.EXPDATA(lines).EXO = 0;
         handles.EXPDATA(lines).ConeHRR=0;
         handles.EXPDATA(lines).ConeMLR=0;
         handles.EXPDATA(lines).ConeEHC=0;
         handles.EXPDATA(lines).type = 'TGA';
         handles.EXPDATA(lines).name = 'TGA ';
         handles.EXPDATA(lines).path = file;
         handles.EXPDATA(lines).check = [0 0 0 0];
         handles.EXPDATA(lines).pair = 0;
         handles.Options.rb(lines,1:2)=0;
         case('TGAFDS')
         handles.EXPDATA(lines).temperature = values(:,1);
         handles.EXPDATA(lines).TGA = values(:,2);
         handles.EXPDATA(lines).time = time;
         handles.EXPDATA(lines).DSC = 0;
         handles.EXPDATA(lines).EXO = 0;
         handles.EXPDATA(lines).ConeHRR=0;
         handles.EXPDATA(lines).ConeMLR=0;
         handles.EXPDATA(lines).ConeEHC=0;
         handles.EXPDATA(lines).type = 'TGAFDS';
         handles.EXPDATA(lines).name = 'TGA ';
         handles.EXPDATA(lines).path = file;
         handles.EXPDATA(lines).check = [0 0 0 0];
         handles.EXPDATA(lines).pair = 0;
         handles.Options.rb(lines,1:2)=0;
      case('DSC') 
         handles.EXPDATA(lines).temperature = values(:,1);
         handles.EXPDATA(lines).DSC = values(:,2);
         handles.EXPDATA(lines).EXO = EXO;
         handles.EXPDATA(lines).TGA=0;
         handles.EXPDATA(lines).gradient = 0;
         handles.EXPDATA(lines).opgradient = 0;
         handles.EXPDATA(lines).ConeHRR=0;
         handles.EXPDATA(lines).ConeMLR=0;
         handles.EXPDATA(lines).ConeEHC=0;
         handles.EXPDATA(lines).type='DSC';
         handles.EXPDATA(lines).name = 'DSC ';
         handles.EXPDATA(lines).path = file;
         handles.EXPDATA(lines).check = [0 0 0 0];
         handles.EXPDATA(lines).pair = 0;
         handles.Options.rb(lines,1:2)=0;
      case('DTA')
         handles.EXPDATA(lines).temperature = values(:,1);
         handles.EXPDATA(lines).TGA = values(:,2);
         handles.EXPDATA(lines).EXO = 0;
         handles.EXPDATA(lines).DSC = 0;
         handles.EXPDATA(lines).ConeHRR=0;
         handles.EXPDATA(lines).ConeMLR=0;
         handles.EXPDATA(lines).ConeEHC=0;
         handles.EXPDATA(lines).type = 'TGA';
         handles.EXPDATA(lines).name = 'TGA ';
         handles.EXPDATA(lines).path = file;
         handles.EXPDATA(lines).check = [0 0 0 0];
         handles.EXPDATA(lines).pair = lines+1;
         handles.Options.rb(lines,1:2)=0;
         %
         lines = lines + 1;
         handles.var.lines = lines;
         %
         handles.EXPDATA(lines).temperature = values2(:,1);
         handles.EXPDATA(lines).DSC = values2(:,2);
         handles.EXPDATA(lines).EXO = EXO;
         handles.EXPDATA(lines).TGA=0;
         handles.EXPDATA(lines).gradient = 0;
         handles.EXPDATA(lines).opgradient = 0;
         handles.EXPDATA(lines).ConeHRR=0;
         handles.EXPDATA(lines).ConeMLR=0;
         handles.EXPDATA(lines).type='DSC';
         handles.EXPDATA(lines).name = 'DSC ';
         handles.EXPDATA(lines).path = file;
         handles.EXPDATA(lines).check = [0 0 0 0];
         handles.EXPDATA(lines).pair = lines-1;
         handles.Options.rb(lines,1:2)=0;
      case('Cone') 
         [m,n] = size(values);
         time = values(:,1);
         handles.EXPDATA(lines).time=removeNan_2(time);
         HRR = values(:,2);
         handles.EXPDATA(lines).ConeHRR= removeNan_2([time,HRR]); %sample 10 cm x 10 cm
         if (n>2)
             MLR = values(:,3);
            handles.EXPDATA(lines).ConeMLR=removeNan_2([time,MLR]);
         end
         if (n>3)
             EHC = values(:,4);
             handles.EXPDATA(lines).ConeEHC=removeNan_2([time,EHC]);
         end
         handles.EXPDATA(lines).temperature=0;
         handles.EXPDATA(lines).TGA=0;
         handles.EXPDATA(lines).EXO = 0;
         handles.EXPDATA(lines).gradient = 0;
         handles.EXPDATA(lines).opgradient = 0;
         handles.EXPDATA(lines).DSC=0;
         handles.EXPDATA(lines).type='Cone';
         handles.EXPDATA(lines).rate='-';
         handles.EXPDATA(lines).name = 'Cone ';
         handles.EXPDATA(lines).path = file;
         handles.EXPDATA(lines).check = [0 0 0 0];
         handles.EXPDATA(lines).pair = 0;
         handles.Options.rb(lines,1:2)=0;
         case('ConeFDS') 
         [m,n] = size(values);
         handles.EXPDATA(lines).time=values(:,1);
         handles.EXPDATA(lines).ConeHRR=values(:,2);
         if (n>2)
            handles.EXPDATA(lines).ConeMLR=values(:,3);
         end
         if (n>3)
             handles.EXPDATA(lines).ConeEHC=values(:,4);
         end
         handles.EXPDATA(lines).temperature=0;
         handles.EXPDATA(lines).TGA=0;
         handles.EXPDATA(lines).EXO = 0;
         handles.EXPDATA(lines).gradient = 0;
         handles.EXPDATA(lines).opgradient = 0;
         handles.EXPDATA(lines).DSC=0;
         handles.EXPDATA(lines).type='ConeFDS';
         handles.EXPDATA(lines).rate='-';
         handles.EXPDATA(lines).name = 'Cone ';
         handles.EXPDATA(lines).path = file;
         handles.EXPDATA(lines).check = [0 0 0 0];
         handles.EXPDATA(lines).pair = 0;
         handles.Options.rb(lines,1:2)=0;
   end

   % GUIs for ReadDATA 
   % this should come up after user selects data to read
   handles.heditfig = figure('Visible','off',...
      'Name', DType, ...
      'Menubar', 'none', ...
      'NumberTitle', 'off', ...
      'Toolbar', 'none', ...
      'Color', get(0,'defaultuicontrolbackgroundcolor'), ...
      'Units', 'normalized', ...
      'Position',[0.36,0.5,0.15,0.4]);

   % Add items for GUI
   handles = AddGuiItems(handles,DType);
   
   % save handles
   guidata(handles.heditfig,handles)
   % Assign the GUI a name to appear in the window title.
   movegui(handles.heditfig,'center')
   % Make the GUI visible.
   set(handles.heditfig,'Visible','on');
   % Wait for ok-button
   uiwait(handles.heditfig)
   % update handles-data
   handles = guidata(handles.heditfig);
   close(handles.heditfig)
   handles = rmfield(handles,'heditfig');
   handles = rmfield(handles,'hname');
   handles = rmfield(handles,'hgaspopup');
   handles = rmfield(handles,'hmaterial');
   handles = rmfield(handles,'hrate');   
   handles = rmfield(handles,'hTGAUpopup');
   handles = rmfield(handles,'hDSCUpopup');
end
end

function handles = AddGuiItems(handles,DType)
lines = handles.var.lines;
%
% add 'display file' check button
hdf = uicontrol(handles.heditfig,'Style','checkbox',...
   'CallBack', @hdf_cb,...
   'String','Display file in new window',...
   'Units', 'normalized', ...
   'Value',0,'Position',[0.1 0.9 0.8 0.05]);
%add edit text box for material
hmaterialtitle = uicontrol(handles.heditfig,'Style','text',...
   'String','Enter material',...
   'Units', 'normalized', ...
   'Position',[0.1 0.82 0.8 0.05]);
handles.hmaterial = uicontrol(handles.heditfig,'Style','edit',...
   'String',handles.var.save_material,...
   'BackgroundColor', 'w',...
   'Units', 'normalized', ...
   'Position',[0.1 0.77 0.8 0.05]);

% add rate input
if strcmp(DType,'Cone')
   ratestr1 = 'Enter heat flux level (kW/m2)';
else
   ratestr1 = 'Enter heating rate K/min';
end
ratestr2 = num2str(handles.var.save_rate);
hratetitle = uicontrol(handles.heditfig,'Style','text',...
   'String',ratestr1,...
   'Units', 'normalized', ...
   'Position',[0.1 0.67 0.8 0.05]);
handles.hrate = uicontrol(handles.heditfig,'Style','edit',...
   'String',ratestr2,...
   'Units', 'normalized', ...
   'BackgroundColor', 'w',...
   'Position',[0.3 0.62 0.4 0.05]);

%add popup menu for gas
hgastitle = uicontrol(handles.heditfig,'Style','text',...
   'String','Purge gas',...
   'Units', 'normalized', ...
   'Position',[0.1, 0.49, 0.38, 0.05]);
handles.hgaspopup = uicontrol(handles.heditfig,'Style','popupmenu',...
   'String',{'Air','N2'},...
   'BackgroundColor', 'w',...
   'Callback', @gas_Callback, ...
   'Units', 'normalized', ...
   'Value',1,...
   'Position',[0.1, 0.44, 0.38, 0.05]);

%add edit text for surface ratio of cone
if strcmp(DType,'Cone') 
    hA_r = uicontrol(gcf,'Style','text',...
   'String','Surface ratio',...
   'Units', 'normalized', ...
   'Position',[0.1, 0.36, 0.38, 0.05]);
    handles.hA_r = uicontrol(gcf,'Style','edit',...
   'String','1',...
   'Units', 'normalized', ...
   'BackgroundColor', 'w',...
   'Position',[0.1, 0.31, 0.38, 0.05]);
end
    
%add popup menu for TGA units
if (strcmp(DType,'TGA') || strcmp(DType,'DTA')) || (strcmp(DType,'TGAFDS'))
   hTGAUtitle = uicontrol(handles.heditfig,'Style','text',...
      'String','TGA units',...
      'Units', 'normalized', ...
      'Position',[0.1, 0.36, 0.38, 0.05]);
   handles.hTGAUpopup = uicontrol(handles.heditfig,'Style','popupmenu',...
      'String',{'Percent (%)','Ratio ()','Mass (g)'},...
      'BackgroundColor', 'w',...
      'Units', 'normalized', ...
      'Value',1,...
      'Position',[0.1, 0.31, 0.38, 0.05]);
end

%add popup menu for DSC Exothermic and units
if (strcmp(DType,'DSC') || strcmp(DType,'DTA'))
   hEXOtitle = uicontrol(handles.heditfig,'Style','text',...
      'String','DSC Exothermic peak',...
      'Units', 'normalized', ...
      'Position',[0.50, 0.48, 0.41, 0.10]);
   ExoVal = 1;
   if (handles.EXPDATA(lines).EXO == -1), ExoVal = 2; end
   handles.hEXOpopup = uicontrol(handles.heditfig,'Style','popupmenu',...
      'String',{'Up (+)','Down (-)'},...
      'BackgroundColor', 'w',...
      'Units', 'normalized', ...
      'Value',ExoVal,...
      'Position',[0.55, 0.44, 0.38, 0.05]);

   hDSCUtitle = uicontrol(handles.heditfig,'Style','text',...
      'String','DSC units',...
      'Units', 'normalized', ...
      'Position',[0.55, 0.36, 0.38, 0.05]);
   handles.hDSCUpopup = uicontrol(handles.heditfig,'Style','popupmenu',...
      'String',{'mW/mg','W/g','W/kg'},...
      'BackgroundColor', 'w',...
      'Units', 'normalized', ...
      'Value',1,...
      'Position',[0.55, 0.31, 0.38, 0.05]);
end

% add the name text
hnametitle = uicontrol(gcf,'Style','text',...
   'String','Name of the data',...
   'Units', 'normalized', ...
   'Position',[0.10 0.18 0.8 0.05]);
handles.hname = uicontrol(gcf,'Style','edit',...
   'String',handles.var.save_name,...
   'CallBack', @name_cb, ...
   'Units', 'normalized', ...
   'BackgroundColor', 'w',...
   'Position',[0.10 0.13 0.8 0.05]);
% add the pushbutton 'OK'
hok = uicontrol(handles.heditfig,'Style','pushbutton',...
   'String','OK',...
   'Units', 'normalized', ...
   'Position',[0.20,0.04,0.3,0.06],...
   'Callback',@okbutton_Callback);
% add the pushbutton 'Cancel'
hcancel = uicontrol(handles.heditfig,'Style','pushbutton',...
   'String','Cancel',...
   'Units', 'normalized', ...
   'Position',[0.55,0.04,0.3,0.06],...
   'Callback',@cancelbutton_Callback);
end

function gas_Callback(hObject, eventdata)
handles = guidata(hObject);
val=get(hObject, 'Value');
str=get(hObject, 'String');
lines = handles.var.lines;
% get type
type = handles.EXPDATA(lines).type;
if strcmp(type,'TGAFDS')
    type = 'TGA (FDS)';
elseif strcmp(type,'ConeFDS')
    type = 'Cone (FDS)';
end
% get gas
gas = char(str(val));
handles.EXPDATA(lines).gas = gas;
% get material
material = get(handles.hmaterial, 'string');
% get rate
rate = num2str(get(handles.hrate, 'string'));
% set name
name = [material ' ' type ' ' rate ' ' gas];
handles.EXPDATA(lines).name = name;
set(handles.hname,'String',name)
guidata(hObject,handles)
end

function name_cb(hObject, eventdata)
handles = guidata(hObject);
name = get(hObject, 'string');
lines = handles.var.lines;
handles.EXPDATA(lines).name = name;
handles.var.save_name = name;
guidata(hObject,handles);
end

function okbutton_Callback(hObject, eventdata)
handles = guidata(hObject);
lines = handles.var.lines;
% get material
material = get(handles.hmaterial, 'string');
handles.EXPDATA(lines).material = material;
handles.var.save_material = material;
% get rate
rate_str = get(handles.hrate,'string');
rate = str2double(rate_str);
%handles.EXPDATA(lines).name = [handles.EXPDATA(lines).name ' ' rate_str];
handles.EXPDATA(lines).rate=rate;
handles.var.save_rate = rate;

% copy data name to its pair in case of DTA data.
pair = handles.EXPDATA(lines).pair;
DTAlines = lines;
if (pair > 0)
   DTAlines = [pair lines];
   handles.EXPDATA(pair).material = handles.EXPDATA(lines).material;
   handles.EXPDATA(pair).rate = handles.EXPDATA(lines).rate;
   handles.EXPDATA(pair).gas      = handles.EXPDATA(lines).gas;
   str = handles.EXPDATA(lines).name;
   nn = strfind(str,'DSC');
   if isempty(nn)
      handles.EXPDATA(pair).name = ['TGA ' handles.EXPDATA(lines).material ...
         ' ' mat2str(handles.EXPDATA(lines).rate) ' ' handles.EXPDATA(lines).gas];
   else
      if (nn>1)
         handles.EXPDATA(pair).name = [str(1:nn-1) 'TGA' str(nn+3:length(str))];
      else
         handles.EXPDATA(pair).name = ['TGA' str(4:length(str))];
      end
   end
end
for i = [DTAlines]
   % Set time
   if (~strcmp(handles.EXPDATA(i).type,'Cone') && ~strcmp(handles.EXPDATA(i).type,'TGAFDS') && ~strcmp(handles.EXPDATA(i).type,'ConeFDS'))
      MinTemp = min(handles.EXPDATA(i).temperature);
      timeMin = (handles.EXPDATA(i).temperature-MinTemp) / handles.EXPDATA(i).rate;
      handles.EXPDATA(i).time=timeMin.*60;
      
   end
   if (strcmp(handles.EXPDATA(i).type,'DSC'))
      % Set DSC direction
      ExoVal = get(handles.hEXOpopup,'value');
      switch ExoVal
         case(1)
            kEXO = +1;
         case(2)
            kEXO = -1;
      end
      % Set DSC units
      DSCUval = get(handles.hDSCUpopup,'value');
      switch DSCUval
         case {1,2} % mW/mg, W/g
            kDSC = 1000.0;
         case 3 % W/kg
            kDSC = 1.0;
      end
      handles.EXPDATA(i).DSC = kEXO * kDSC * handles.EXPDATA(i).DSC;
   end
   if (strcmp(handles.EXPDATA(i).type,'TGA'))
      % Set TGA units
      TGAUval = get(handles.hTGAUpopup,'value');
      switch TGAUval
         case 1% Percent (%)
            kTGA = 1/100;
         case 2 % Ratio
            kTGA = 1.0;
         case 3 % Mass (g)
            kTGA = 1/max(handles.EXPDATA(i).TGA);
      end
      handles.EXPDATA(i).TGA = kTGA * handles.EXPDATA(i).TGA;
      % Compute TGA gradients
      handles.EXPDATA(i).gradient = gradient(handles.EXPDATA(i).TGA, handles.EXPDATA(i).time(2)-handles.EXPDATA(i).time(1));
      handles.EXPDATA(i).opgradient = -handles.EXPDATA(i).gradient;
   end
   if (strcmp(handles.EXPDATA(i).type,'TGAFDS'))
   % Compute TGA gradients
      handles.EXPDATA(i).gradient = gradient(handles.EXPDATA(i).TGA, handles.EXPDATA(i).time(2)-handles.EXPDATA(i).time(1));
      handles.EXPDATA(i).opgradient = -handles.EXPDATA(i).gradient;
   end
   if strcmp(handles.EXPDATA(i).type,'Cone')
       A_str = get(handles.hA_r,'string');
       A_r = str2double(A_str);
      handles.EXPDATA(i).ConeHRR =  handles.EXPDATA(i).ConeHRR.*A_r;
      handles.EXPDATA(i).ConeMLR =  handles.EXPDATA(i).ConeMLR.*A_r;
   end
end
guidata(hObject,handles);
uiresume(handles.heditfig);
end

function cancelbutton_Callback(hObject,eventdata)
handles = guidata(hObject);
lines = handles.var.lines;
% copy data name to its pair in case of DTA data.
pair = handles.EXPDATA(lines).pair;
if (pair > 0)
   Nmax = pair-1;
else
   Nmax = lines-1;
end
handles.EXPDATA = handles.EXPDATA(1:Nmax);
handles.var.lines = Nmax;
guidata(hObject,handles);
uiresume(handles.heditfig);
end

function [headerlines,datalines,EXO] = getHeaderlines(filename)
[title, value]=textread(filename, '%s %s', 'delimiter', ';');
titlelist=char(title);
headerlines=0;
index = 1;
EXO = +1;
ExoFound = 0;
for i=1:length(titlelist)
   str = titlelist(i,:);
   if ~ExoFound
      n = strfind(str,'EXO');
      if (~isempty(n))
         len = length(str);
         EXO = str2double(str(n+4:len));
         ExoFound = 1;
      end
   end
   if strcmp(str(1), '#')
      index = index + 1;
   else
      headerlines = index;
      break;
   end
end
datalines = length(titlelist)-headerlines;
end

%callbacks of question box

function hdf_cb(hObject, eventdata)
handles = guidata(hObject);
lines = handles.var.lines;
open(handles.EXPDATA(lines).path);
end
