function varargout = mygui(varargin)
% Grad Shafranov reconstruction user interface
%
% Christian Möstl, Space Research Institute, Graz, Austria
% Based on the code by Qiang Hu, PhD thesis, Dartmouth College, 2002
% last update January 2010
%
%Reconstruction Pipeline:
%
%1. use "getace.m"
%type "mygui" to start a graphical user interface
%2. choose event
%
%use buttons to...
%3. Plot data (calls "dataplot.m")
%4. dht and MVUB analysis ("prefr.m")  
%5. GS axis  ("optcloud.m", calls several other functions)
%6. GS solver ("hu12n.m", calls "hu34n.m")
%
%
%figures and results are saved in the event directories
%e.g. ACE_303_306_05
%
%(means event is from ACE spacecraft, full data interval selected from 
%day of year 303 to 306 in year 2005)
%
%
%The files
%"ALLCASE.par" and "CASE1.par" (in Event - Directory)
%control the reduction process!
%
%ALLCASE.PAR looks like this: 
%"
%ACE_303_306_05   %event directory
%0			% 1 to write data file
%3			% porder - order of fitting polynomial
%21 131 16		% nx ny (-)dmid
%1 0.1 0.15		% get_Ab (0:yes; -1: sum for < Ab; +1: sum for > Ab) (+)dAl0 (-)dAr0
%"
%These parameters can be controlled using "mygui".


% Edit the above text to modify the response to help mygui


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mygui_OpeningFcn, ...
                   'gui_OutputFcn',  @mygui_OutputFcn, ...
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


% --- Executes just before mygui is made visible.
function mygui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mygui (see VARARGIN)

%path('C:\chris\Matlab\GS-code',path)
%cd('C:\chris\Matlab\GS-code')
handles.basedir=pwd;

% Choose default command line output for mygui
handles.output = hObject;


f0=fopen('ALLCASE.par', 'rt');
%dirstr=fscanf(f0, '%s\n')
dirstr=fscanf(f0, '%c', 14);
ctrl=fscanf(f0, '%i', 6);
write=ctrl(1);
porder=ctrl(2); nx=ctrl(3); ny=ctrl(4); dmid=ctrl(5); get_Ab=ctrl(6);
ctrl2=fscanf(f0, '%f', 2);
dAl0=ctrl2(1); dAr0=ctrl2(2);
fclose(f0);

handles.dirstring=dirstr;
handles.polyorder=porder;
handles.writefile=write;
handles.ctrl11=nx;
handles.ctrl12=ny;
handles.ctrl13=dmid;
handles.ctrl21=get_Ab;
handles.ctrl22=dAl0;
handles.ctrl23=dAr0;

pvalue=porder-1;

set(handles.edit3,'String', dirstr);
set(handles.popupmenu2,'Value', write+1);
set(handles.popupmenu3,'Value', pvalue);
set(handles.edit5,'String', num2str(nx));
set(handles.edit7,'String', num2str(ny));
set(handles.edit8,'String', num2str(dmid));
set(handles.edit9,'String', num2str(get_Ab));
set(handles.edit10,'String', num2str(dAl0));
set(handles.edit11,'String', num2str(dAr0));

cd(handles.dirstring);
fpar=fopen('CASE1.par', 'rt');
str=fscanf(fpar, '%c', 21); % str= Dateinamen der Daten, z.B.:..\data\winall291_95.mat
while feof(fpar) ~= 1
i12=fscanf(fpar,'%i %i %i\n', 3); 
end
fclose(fpar);
handles.i1=i12(1); % Anfangspunkt der Daten
handles.i2=i12(2); % Endpunkt der Daten
handles.adjp=i12(3);


if exist('zs.dat','file') ~= 0
 load zs.dat;
 [axislong,axislat]=vectodeg(zs);
 handles.axislat=axislat;
 handles.axislong=axislong;
else 
 axislat=0;
 axislong=0;
 handles.axislat=0;
 handles.axislong=0;
end

set(handles.edit12,'String',handles.i1);
set(handles.edit14,'String',handles.i2);
set(handles.popupmenu4,'Value',handles.adjp+1);

set(handles.edit23,'String', num2str(axislat));
set(handles.edit25,'String', num2str(axislong));


cd ..





% Update handles structure
guidata(hObject, handles)

% UIWAIT makes mygui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mygui_OutputFcn(hObject, eventdata, handles) 
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

prefr;
cd(handles.dirstring);
load dht.dat;
load CCoef;
load walen;
load mvaby.dat;
load zeitraum;
load la23;
s=size(zeitraum);
load d_in_AU;
[long,lat]=vectodeg(mvaby);
out1=['dHT velocity (xyz):  ', num2str(dht,'%10.2f')];
out2=['dHT Correlation:  ',  num2str(CCoef,'%10.5f')];
out3=['Walen slope:  ',  num2str(walen_slope,'%10.5f')];
out4=['MVAB orientation: longitude: ',  num2str(long,'%4.2f'),' latitude: ',num2str(lat,'%4.2f')];
out5=['Time: ',zeitraum(5:s(2))];
out6=['Diameter (AU): ',  num2str(d_in_AU,'%4.4f')];
out7=['Eigenvalue ratio 2/3: ',  num2str(Lambda2_zu_Lambda3,'%4.4f')];
outstring=strvcat(out5, out1,out2,out3,out4, out6, out7);
set(handles.text6,'String', outstring);

cd ..

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataplot;

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton3.
function pushbutton3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;
cd(handles.basedir);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

optcloud;

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hu12n;

cd(handles.dirstring);
load solvervariables;
[long, lat]=vectodeg(zs);
out0=['GS orientation: longitude: ',  num2str(long,'%4.2f'),' latitude: ',num2str(lat,'%4.2f')];
out1=['Average Beta:  ', num2str(ave_beta,'%10.3f')];
out2=['Mean speed to Alfven speed:  ',  num2str(Mta,'%10.4f')];
out3=['Fitting residue:  ',  num2str(residuef,'%2.4f')];
out4=['Bzmax (nT):  ',  num2str(Bzmax,'%4.1f'),'     Y0 (AU):  ',num2str(Y0,'%2.4f')];
out5=['Iz (MA):  ',num2str(round(current/1e6),'%6.0f')];
out6=['Axíal flux (10^21 Mx):  ',  num2str(axflux, '%6.3f')];
out7=['Poloidal flux (for L=1 AU) (10^21 Mx):  ',  num2str(polflux, '%6.3f')];
outstring=strvcat(out0,out1,out2,out3,out4, out5, out6,out7);
set(handles.text6,'String', outstring);
cd ..


% --------------------------------------------------------------------
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


d=dir;
s=size(d);

for k = 1:s(1)
    qdir{k} = getfield(d, {k,1}, 'isdir');
end
qdir=cell2mat(qdir);

for k = 1:s(1)
    namefiles{k} = getfield(d, {k,1}, 'name');
end
%namdir=cell2mat(namdir)
%namdirind=find(qdir);
namdir=namefiles(find(qdir));

s=size(namdir);
for k = 1:s(2) 
    dirchar=cell2mat(namdir(k));
    sn=size(dirchar);
    if (sn(2) == 14) | (sn(2) == 15) 
      qlength(k)=1 ; 
     else 
      qlength(k)=0;
     end
    end

namdir=namdir(find(qlength));

string_list= namdir;
set(hObject,'String',string_list);
guidata(hObject,handles);








val = get(hObject,'Value');
string_list = get(hObject,'String');
handles.dirstring = string_list{val}; % Convert from cell array
                                   % to string

set(handles.edit3,'String',handles.dirstring);
cd(handles.dirstring);

pathtolook=[handles.dirstring,'\ALLCASE.par']
fileex=exist(pathtolook,'file')

pwd

if fileex == 2
   copyfile('ALLCASE.par',handles.basedir)
   
   

%ACE_316_319_04
%1
%2
%19 221 0
%0 0.000000 0.100000
   
f0=fopen('ALLCASE.par', 'rt');
%dirstr=fscanf(f0, '%s\n')
dirstr=fscanf(f0, '%c', 14);
ctrl=fscanf(f0, '%i', 6);
write=ctrl(1);
porder=ctrl(2); nx=ctrl(3); ny=ctrl(4); dmid=ctrl(5); get_Ab=ctrl(6);
ctrl2=fscanf(f0, '%f', 2);
dAl0=ctrl2(1); dAr0=ctrl2(2);
fclose(f0);

handles.dirstring=dirstr;
handles.polyorder=porder;
handles.writefile=write;
handles.ctrl11=nx;
handles.ctrl12=ny;
handles.ctrl13=dmid;
handles.ctrl21=get_Ab;
handles.ctrl22=dAl0;
handles.ctrl23=dAr0;

pvalue=porder-1;

set(handles.edit3,'String', dirstr);
set(handles.popupmenu2,'Value', write+1);
set(handles.popupmenu3,'Value', pvalue);
set(handles.edit5,'String', num2str(nx));
set(handles.edit7,'String', num2str(ny));
set(handles.edit8,'String', num2str(dmid));
set(handles.edit9,'String', num2str(get_Ab));
set(handles.edit10,'String', num2str(dAl0));
set(handles.edit11,'String', num2str(dAr0));
     
end



fpar=fopen('CASE1.par', 'rt');
str=fscanf(fpar, '%c', 21); % str= Dateinamen der Daten, z.B.:..\data\winall291_95.mat
while feof(fpar) ~= 1
i12=fscanf(fpar,'%i %i %i\n', 3); 
end
fclose(fpar);
handles.i1=i12(1); % Anfangspunkt der Daten
handles.i2=i12(2); % Endpunkt der Daten
handles.adjp=i12(3);


set(handles.edit12,'String',handles.i1);
set(handles.edit14,'String',handles.i2);
set(handles.popupmenu4,'Value',handles.adjp+1);

cd ..

guidata(hObject,handles);
% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
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

d=dir;
s=size(d);

for k = 1:s(1)
    qdir{k} = getfield(d, {k,1}, 'isdir');
end
qdir=cell2mat(qdir);

for k = 1:s(1)
    namefiles{k} = getfield(d, {k,1}, 'name');
end
%namdir=cell2mat(namdir)
%namdirind=find(qdir);
namdir=namefiles(find(qdir));

s=size(namdir);
for k = 1:s(2) 
    dirchar=cell2mat(namdir(k));
    sn=size(dirchar);
    if (sn(2) == 14) | (sn(2) == 15) 
      qlength(k)=1 ; 
     else 
      qlength(k)=0;
     end
    end

namdir=namdir(find(qlength));

string_list= namdir;
set(hObject,'String',string_list);
guidata(hObject,handles);



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
control=[handles.ctrl11 handles.ctrl12 handles.ctrl13 handles.ctrl21];
control2=[handles.ctrl22 handles.ctrl23];
changeallcase2(handles.dirstring, handles.writefile, handles.polyorder,control,control2);


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.dirstring=get(hObject,'String');% returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

f0=fopen('ALLCASE.par', 'rt');
dirstr=deblank(fgets(f0));
fclose(f0);

handles.dirstring=dirstr;
set(hObject,'String',handles.dirstring);
guidata(hObject,handles);
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%polyorder=contents{get(hObject,'Value')} %returns selected item from popupmenu2

val = get(hObject,'Value');
switch val
case 1
% User selected the first item
handles.writefile=0;
case 2
handles.writefile=1;
end

guidata(hObject,handles);
    
% User selected the second item
% Proceed with callback...


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%set(hObject,'Value',2)
%handles.writefile=1;
guidata(hObject,handles);

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%  writefile= contents{get(hObject,'Value')} %returns selected item from popupmenu3
val = get(hObject,'Value');
switch val
case 1
handles.polyorder=2;
case 2
handles.polyorder=3;   
case 3    
handles.polyorder=4;  
case 4
handles.polyorder=5;  
case 5
handles.polyorder=6;  
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%set(hObject,'Value',2);

%handles.polyorder=3;
%guidata(hObject,handles);

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
handles.ctrl11=str2double(get(hObject,'String')); % returns contents of edit5 as a double
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%handles.ctrl11=21;

%guidata(hObject,handles);

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
handles.ctrl12=   str2double(get(hObject,'String'));% returns contents of edit7 as a double
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%handles.ctrl12=131;
%guidata(hObject,handles);
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
handles.ctrl13=  str2double(get(hObject,'String')) %returns contents of edit8 as a double
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%handles.ctrl13=0;
%guidata(hObject,handles);

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
handles.ctrl21= str2double(get(hObject,'String')); %returns contents of edit9 as a double
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%handles.ctrl21=1;
%guidata(hObject,handles);

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
handles.ctrl22= str2double(get(hObject,'String'));% returns contents of edit10 as a double
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%handles.ctrl22=0.1;
%guidata(hObject,handles);

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
handles.ctrl23= str2double(get(hObject,'String'))% returns contents of edit10 as a double
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.ctrl23=0.1;
guidata(hObject,handles);


% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
handles.i1=str2double(get(hObject,'String')); % returns contents of edit5 as a double
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
handles.i2=str2double(get(hObject,'String')); % returns contents of edit5 as a double
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(hObject,'Value');
switch val
case 1
% User selected the first item
handles.adjp=0;
case 2
handles.adjp=1;
end

guidata(hObject,handles);


% Hints: contents = get(hObject,'String') returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cd(handles.dirstring);
parid=fopen('CASE1.par', 'at');
fprintf(parid, '%i %i %i\n', [handles.i1 handles.i2 handles.adjp]);
fclose(parid);
cd('..');





function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.basedir);



% --- Executes during object creation, after setting all properties.
function text6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes during object creation, after setting all properties.
function pushbutton9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

control=[handles.ctrl11 handles.ctrl12 handles.ctrl13];
control2=[handles.ctrl21 handles.ctrl22 handles.ctrl23];
changeallcase3(handles.dirstring, handles.writefile, handles.polyorder,control,control2);


% --------------------------------------------------------------------
function newace_Callback(hObject, eventdata, handles)
% hObject    handle to newace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

acegui
%popupmenu1_CreateFcn








function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double

handles.axislat=str2double(get(hObject,'String'));
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


handles.axislong=str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cd(handles.dirstring);

zs=degtovec(handles.axislong,handles.axislat);

save zs.dat zs -ascii;
cd ..

