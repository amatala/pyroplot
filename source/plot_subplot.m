function handles = plot_subplot(handles);

handles.SP.RowStr = '1';
handles.SP.ColumnStr = '1';
handles.SP.Rows = 1;
handles.SP.Columns = 1;
%handles.SP.Cb = [];
handles.SP.PU_first = [];
handles.SP.PU_second = [];
handles.SP.PU_third = [];
handles.SP.fig = figure('Visible','off',...
   'Name', 'Plot subplot', ...
   'NumberTitle', 'off', ...
   'Menubar', 'none', ...
   'Toolbar', 'none', ...
   'Color', get(0,...
   'defaultuicontrolbackgroundcolor'), ...
   'Units', 'normalized', ...
   'Position',[0.36,0.3,0.4,0.7]);
hSP = handles.SP.fig;
movegui(hSP,'center')


hRow = uicontrol(hSP,'Style','text',...
   'String','Rows',...
   'Units', 'normalized', ...
   'Position',[0.25 0.955 0.1 0.03]);

hRowEdit = uicontrol(hSP,'Style','edit',...
     'BackgroundColor', 'w',...
     'String', '1', ...
     'Units', 'normalized', ...
     'Callback', {@row_cb}, ...
     'Position',[0.35 0.96 0.07 0.03]);
  
hColumn = uicontrol(hSP,'Style','text',...
   'String','Columns',...
   'Units', 'normalized', ...
   'Position',[0.5 0.955 0.1 0.03]);

hColumnEdit = uicontrol(hSP,'Style','edit',...
     'BackgroundColor', 'w',...
     'String', '1', ...
     'Units', 'normalized', ...
     'Callback', {@column_cb}, ...
     'Position',[0.6 0.96 0.07 0.03]);

 lines = handles.var.lines;
 d = 0.045;
 h = 0.03;
 for i =1:lines
     %handles.SP.check(i) = 0;
     %handles.SP.Cb(i) = uicontrol(hSP,'Style','checkbox',...
     % 'CallBack', {@Sp_cb,i},...
     % 'Units', 'normalized', ...
     % 'Value',handles.SP.check(i),...
     % 'Position',[0.1 0.95-i*d 0.03 h]);
  %guidata(handles.SP.Cb(i),handles);
  type = handles.EXPDATA(i).type;
     hName = uicontrol(hSP,'Style','text',...
    'String', handles.EXPDATA(i).name,...
    'Units', 'normalized', ...
    'Position',[0.05 0.943-i*d 0.3 h]);

    handles.SP.PU_first(i) = uicontrol(hSP, 'Style', 'popupmenu', ...
        'String','No',...
        'BackgroundColor', 'w',...
        'Units', 'normalized', ...
        'Value',1,...
        'Position',[0.45 0.95-i*d 0.1 h]);
    guidata(handles.SP.PU_first(i), handles);
   handles.SP.PU_second(i) = uicontrol(hSP, 'Style', 'popupmenu', ...
        'String','No',...
        'BackgroundColor', 'w',...
        'Units', 'normalized', ...
        'Value',1,...
        'Position',[0.65 0.95-i*d 0.1 h]);
     guidata(handles.SP.PU_second(i), handles);
     
     handles.SP.PU_third(i) = uicontrol(hSP, 'Style', 'popupmenu', ...
        'String','No',...
        'BackgroundColor', 'w',...
        'Units', 'normalized', ...
        'Value',1,...
        'Position',[0.85 0.95-i*d 0.1 h]);
     guidata(handles.SP.PU_third(i), handles);
     
        if strcmp(type,'TGA') || strcmp(type,'TGAFDS')
         hTitle1 = uicontrol(hSP,'Style','text',...
        'String', '(%)',...
        'Units', 'normalized', ...
        'Position',[0.37 0.943-i*d 0.07 h]);
        hTitle2 = uicontrol(hSP,'Style','text',...
        'String', '(%/s)',...
        'Units', 'normalized', ...
        'Position',[0.57 0.943-i*d 0.07 h]);
        hTitle3 = uicontrol(hSP,'Style','text',...
        'String', '(-%/s)',...
        'Units', 'normalized', ...
        'Position',[0.77 0.943-i*d 0.07 h]);
        elseif strcmp(type,'Cone') || strcmp(type,'ConeFDS')
        hTitle1 = uicontrol(hSP,'Style','text',...
        'String', 'HRR',...
        'Units', 'normalized', ...
        'Position',[0.37 0.943-i*d 0.07 h]);
        hTitle2 = uicontrol(hSP,'Style','text',...
        'String', 'MLR',...
        'Units', 'normalized', ...
        'Position',[0.57 0.943-i*d 0.07 h]);
        hTitle3 = uicontrol(hSP,'Style','text',...
        'String', 'EHC',...
        'Units', 'normalized', ...
        'Position',[0.77 0.943-i*d 0.07 h]); 
        else
        hTitle = uicontrol(hSP,'Style','text',...
        'String', 'DSC',...
        'Units', 'normalized', ...
        'Position',[0.37 0.943-i*d 0.07 h]);
        end
      
 end
 hPlot = uicontrol(hSP,'Style','pushbutton',...
     'String', 'Plot', ...
     'Units', 'normalized', ...
     'Callback', {@plot_cb}, ...
     'Position',[0.4 0.02 0.07 0.04]);
 hClose = uicontrol(hSP,'Style','pushbutton',...
     'String', 'Close', ...
     'Units', 'normalized', ...
     'Callback', {@close_plot_cb}, ...
     'Position',[0.53 0.02 0.07 0.04]);
 
set(hSP, 'Visible', 'on');

uiwait(hSP);
% update handles-data
handles = guidata(hSP);
guidata(handles.hPyroPlot, handles);
close(hSP);

end

function row_cb(hObject, eventdata)
handles = guidata(hObject);

str = get(hObject, 'String');
Rows = str2double(str);
Columns = handles.SP.Columns;

str2 = '''No''';
for i=1:Rows
    for j = 1:Columns
   str2 = [str2, ',' , '''', num2str(i), ', ', num2str(j), '''']; 
    end
end

handles.SP.Rows = Rows;
s1 = ['set(handles.SP.PU_first(:), ''String'', {', str2, '});'];
eval(s1);
s2 = ['set(handles.SP.PU_second(:), ''String'', {', str2, '});'];
eval(s2);
s3 = ['set(handles.SP.PU_third(:), ''String'', {', str2, '});'];
eval(s3);

guidata(hObject,handles);
end

function column_cb(hObject, eventdata)
handles = guidata(hObject);

str = get(hObject, 'String');
Columns = str2double(str);
Rows = handles.SP.Rows;
str2 = '''No''';
for i=1:Rows
    for j = 1:Columns
   str2 = [str2, ',' , '''', num2str(i), ', ', num2str(j), '''']; 
    end
end

handles.SP.Columns = Columns;

s1 = ['set(handles.SP.PU_first(:), ''String'', {', str2, '});'];
eval(s1);
s2 = ['set(handles.SP.PU_second(:), ''String'', {', str2, '});'];
eval(s2);
s3 = ['set(handles.SP.PU_third(:), ''String'', {', str2, '});'];
eval(s3);

guidata(hObject,handles);
end

function plot_cb(hObject, eventdata)
handles = guidata(hObject);
new_plot = figure('Visible','off');
hold on;
m=1;
lines = handles.var.lines;

for i=1:handles.SP.Rows
   for j=1:handles.SP.Columns
       subplot(handles.SP.Rows,handles.SP.Columns, m);
       m = m+1;
       hold on
       d = '';
      for k=1:lines
         if get(handles.SP.PU_first(k),'Value') == m
             if strcmp(handles.EXPDATA(k).type,'TGA') || strcmp(handles.EXPDATA(k).type,'TGAFDS')
                x = handles.EXPDATA(k).temperature;
                y = handles.EXPDATA(k).TGA;
             elseif strcmp(handles.EXPDATA(k).type,'Cone') || strcmp(handles.EXPDATA(k).type,'ConeFDS')
                x = handles.EXPDATA(k).time;
                y = handles.EXPDATA(k).ConeHRR;
             else
                x = handles.EXPDATA(k).temperature;
                y = handles.EXPDATA(k).DSC;
             end %if
             d = [d, mat2str(x),',', mat2str(y), ','];
             
         end %if this subplot
         if get(handles.SP.PU_second(k),'Value') == m
             if strcmp(handles.EXPDATA(k).type,'TGA') || strcmp(handles.EXPDATA(k).type,'TGAFDS')
                x = handles.EXPDATA(k).temperature;
                y = handles.EXPDATA(k).TGA;
             elseif strcmp(handles.EXPDATA(k).type,'Cone') || strcmp(handles.EXPDATA(k).type,'ConeFDS')
                x = handles.EXPDATA(k).time;
                y = handles.EXPDATA(k).ConeHRR;
             else
                x = handles.EXPDATA(k).temperature;
                y = handles.EXPDATA(k).DSC;
             end %if
             d = [d, mat2str(x),',', mat2str(y), ','];
             
         end %if this subplot
         if get(handles.SP.PU_third(k),'Value') == m
             if strcmp(handles.EXPDATA(k).type,'TGA') || strcmp(handles.EXPDATA(k).type,'TGAFDS')
                x = handles.EXPDATA(k).temperature;
                y = handles.EXPDATA(k).TGA;
             elseif strcmp(handles.EXPDATA(k).type,'Cone') || strcmp(handles.EXPDATA(k).type,'ConeFDS')
                x = handles.EXPDATA(k).time;
                y = handles.EXPDATA(k).ConeHRR;
             else
                x = handles.EXPDATA(k).temperature;
                y = handles.EXPDATA(k).DSC;
             end %if
             d = [d, mat2str(x),',', mat2str(y), ','];
             
         end %if this subplot
         end %k
         d = d(1:length(d)-1);
         if ~isempty(d)
         s = ['plot(', d, ')'];
         eval(s);
         end
         hold on;
      end %j
   end %i

set(new_plot, 'Visible','on');
guidata(hObject, handles);
end

function close_plot_cb(hObject, eventdata)
handles = guidata(hObject);
uiresume(handles.SP.fig);
end
