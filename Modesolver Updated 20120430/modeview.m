function varargout = modeview(varargin)
% MODEVIEW - GUI viewer for modesolver output data
%
% Syntax:   modeview(SZ)
%           where SZ = data structure containing mode data from solver

% -------------------------------------------------------------------------
% MODEVIEW M-file for modeview.fig
%      MODEVIEW, by itself, creates a new MODEVIEW or raises the existing
%      singleton*.
%
%      H = MODEVIEW returns the handle to a new MODEVIEW or the handle to
%      the existing singleton*.
%
%      MODEVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEVIEW.M with the given input arguments.
%
%      MODEVIEW('Property','Value',...) creates a new MODEVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before modeview_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to modeview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2005, Milos Popovic
%
% Code to do list:
% ----------------
% Dec 18, 2005 - Make mouse cursor show field amplitude at point where
%                pointing - maybe in a yellow hint bar?.  Or, on moving (or
%                click and moving) of mouse, update field amplitude and
%                maybe also x,y coordinates in some text boxes.
%              - Make field plot right-click menu enable: animate, spawn
%                static figure, spawn/export movie figure.
% Aug 29, 2005 - make one shuffle function to call for 4 places of mode
%                reshuffling (2 for erasing, 2 for mode moving fwd/back)
%              - fcn to initialize new data structure..
%              - menu item to sort modes (by ascend/descending frequency, Q, etc.)
%              - menu item to reorder modes (arbitrarily) by entering vector - similar to
%                existing mode moving
%              - option to append onto existing data structure on open/import
%              - ability to add/update notes on each mode computation.
% **           - REMOVE insertion of F.dx, F.dy prior to call of groupindex
%                function (now groupindex deals with it on its own).
%
% Code updates:
% -------------
% Nov 14, 2011 - Widened some popup boxes to make them visible on OS X
%              - Added colorbar
%              - Introduced colormaps
%              - Made a pass through to clean up code style (Matlab yellow warnings)
% Jul 20, 2006 - Commented out 'EraseMode' and 'DoubleBuffer' settings in
%                axes1_CreateFcn function because the axes don't seem to
%                have this property and cause an error at least in Matlab R14
%                (not in R13).
% Dec 17, 2005 - Corrected printing of path ("\"s) when saving data to file.
% Nov 17, 2005 - Added zoom slider for quicker field zooming.
% Oct 13, 2005 - Bug fix: in MPupdateplot, line taking real/imag/abs of field
%                changed from field(:,:,nmode) to field, since nmode
%                already chosen in previous lines.
%              - Added support for 'Fs' named field struct (besides F and Fb).
%              - Bug fix: one strdata was missing {} cell index on assignment
% Sep 11, 2005 - Added support for plotting Ste,Stm (TE/TM Poynting vectors)
%              - Added support for plotting non-uniform grid modes
%                using 'pcolorcenter'. Can choose shading type
%                (flat, faceted, interp).
% Sep  9, 2005 - v0.12
% Aug 30, 2005 - opening, saving, import/export done.
% Aug 29, 2005 - v0.11 - erasing, struct/mode moving done.


% Edit the above text to modify the response to help modeview

% Last Modified by GUIDE v2.5 14-Nov-2011 21:21:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @modeview_OpeningFcn, ...
                   'gui_OutputFcn',  @modeview_OutputFcn, ...
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


% --- Executes just before modeview is made visible.
function modeview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to modeview (see VARARGIN)

% Choose default command line output for modeview
handles.output = hObject;

handles.releasever  = 'v0.14';                                          % [MP] Versioning info
handles.releasedate = 'Nov 14, 2011';
set(handles.figure_main, 'Name', ['Mode Viewer (' handles.releasever ', ' handles.releasedate ') - report bugs to milos.popovic@colorado.edu']); % Main window title

%set(handles.axes1,'Visible','off');
%set(handles.axes1,'XTick',[]); set(handles.axes1,'YTick',[]);

if (nargin <= 3) || ischar(varargin{1})                                  % [MP] store passed-in SZ variable containing mode data
    error('modeview: you must call modeview with the SZ mode structure as parameter, modeview(SZ)');
end
handles.modedata = varargin{1};                                         % |
MPmodedata_initialize(hObject, eventdata, handles);                     % Initialize GUI with data structure
%handles = guidata(hObject);    % Retrieve handles updated by MPmodedata.. for further commands placed below

% UIWAIT makes modeview wait for user response (see UIRESUME)
% uiwait(handles.figure_main);


% --- MP Initializes GUI to a newly imported datastructure SZ (initialize
% sliders, mode plots, text boxes, popup menus, etc.)
function MPmodedata_initialize(hObject, eventdata, handles)
if(~isfield(handles.modedata,'F'))
    if(isfield(handles.modedata,'Fb'))
        for k = 1:length(handles.modedata)
            handles.modedata(k).F = handles.modedata(k).Fb;
        end
        handles.modedata = rmfield(handles.modedata,'Fb');
        fprintf('modeview: Renaming field datastructure from Fb to F\n');
    elseif(isfield(handles.modedata,'Fs'))
        for k = 1:length(handles.modedata)
            handles.modedata(k).F = handles.modedata(k).Fs;
        end
        handles.modedata = rmfield(handles.modedata,'Fs');
        fprintf('modeview: Renaming field datastructure from Fs to F\n');
    else
        delete(handles.figure_main);    % Kill the current browser instance
        error('modeview: Must have either Fb or F structure containing mode fields');
    end
end
handles.fieldcomp = 'Ex';
handles.fieldfunc = @real;
handles.veclabels = {'wc','beta','betas','betaR','betabR','gamma','dB90','dBcm','dBum','LossQ',...
                     'dB90TE','dB90TM','betabRTE','betabRTM'};   % Fields that, if present, are vectors indexed by mode number
handles.fldlabels = {'Ex','Ey','Ez','Hx','Hy','Hz'};   % Fields that, if present, are vectors indexed by mode number
% Coordinates for fields/values listed in drop box popupmenu_fieldcomp:
%handles.fldcomps  = {'Ex','Ey','Ez','Hx','Hy','Hz','Ste','Stm','Int','We','Wm'};
handles.xcoord     = {'Rr','Rz','Rz','Rz','Rr','Rr','Rr','Rz','Rr','Rr','Rr'};  % Last three are dummy settings (don't know size of Int,We,Wm yet)
handles.ycoord     = {'Zr','Zz','Zr','Zz','Zr','Zz','Zr','Zz','Zr','Zr','Zr'};

% Reset non-uniform grid plotting checkbox
handles.shading  = 'flat';
handles.colormap = 'redbluehilight';
set(handles.checkbox_grid,'Value',0);
MPpopupmenu_shading_Refresh(hObject, eventdata, handles);

fprintf('Building sliders...\n');
NSTR = length(handles.modedata);                                                              % [MP] later chg this line and one below to call to MPslider_str_Refresh
MPslider_set(handles.slider_str,  1, 1, NSTR, handles.text_strnum);
MPslider_str_Refresh(hObject, eventdata, handles);                              % Called to initialize 'erase buttons' gray-out properly
set(handles.slider_mode, 'Value', 1);
%NMODES = size(SZ(1).F.Ex,3);                                                    % [MP]
%MPslider_set(handles.slider_mode, [], 1, NMODES, handles.text_modenum);
MPslider_mode_Refresh(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

%slider_mode_Callback(hObject, eventdata, handles);                      % DOES NOT SEEM LEGIT -- Matlab doesn't change these when I rename the tags!
%slider_str_Callback(hObject, eventdata, handles);                       % DOES NOT SEEM LEGIT -- Matlab doesn't change these when I rename the tags!
popupmenu_fieldcomp_Callback(hObject, eventdata, handles);  handles = guidata(hObject); % DOES NOT SEEM LEGIT
popupmenu_real_Callback(hObject, eventdata, handles);  handles = guidata(hObject);  % DOES NOT SEEM LEGIT
popupmenu_shading_Callback(hObject, eventdata, handles);  handles = guidata(hObject);  % DOES NOT SEEM LEGIT
MPupdate_mode_plot(hObject, eventdata, handles);                                  % [MP]
MPupdate_params_textbox(hObject, eventdata, handles);
MPupdate_eigenvals_textbox(hObject, eventdata, handles);


% --- Outputs from this function are returned to the command line.
function varargout = modeview_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function slider_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function slider_mode_Callback(hObject, eventdata, handles)
% hObject    handle to slider_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
MPslider_mode_Refresh(hObject, eventdata, handles);

MPupdate_mode_plot(hObject, eventdata, handles);
MPupdate_params_textbox(hObject, eventdata, handles);
MPupdate_eigenvals_textbox(hObject, eventdata, handles);


% --- MP mode slider refresh
% - reads NMODES from data structure selected by structure slider
% - verifies/limits current mode value
% - writes current mode, min/max, and text '1/N' to mode slider and text box structures
%
function MPslider_mode_Refresh(hObject, eventdata, handles)
%fprintf('Trying to get structure %f\n', get(handles.slider_str,'Value'));
NMODES = size(handles.modedata( round( get(handles.slider_str,'Value') ) ).F.Ex, 3);   % [MP] ***shouldn't need round, but when slider2 is *dragged*, it seems to end up non-integer despite all the integerizing code in the callback fcn for slider2
if (NMODES > 1),  modeidx = round(get(handles.slider_mode,'Value'));  else  modeidx = 1;  end;
if (modeidx > NMODES), modeidx = NMODES; end;            % Relevant when the structure changes to one with fewer modes
MPslider_set(handles.slider_mode, modeidx, 1, NMODES, handles.text_modenum);
if (modeidx == 1),       st = 'off';  else  st = 'on';  end;  set(handles.pushbutton_modeshiftleft, 'Enable', st);
if (modeidx == NMODES),  st = 'off';  else  st = 'on';  end;  set(handles.pushbutton_modeshiftright, 'Enable', st);
NSTR = length(handles.modedata);
if (NSTR == 1 && NMODES == 1),  st = 'off';  else  st = 'on';  end;  set(handles.pushbutton_modeerase, 'Enable', st);
set(handles.menu_edit_delstructs, 'Enable', st);


% --- Executes during object creation, after setting all properties.
function slider_str_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function slider_str_Callback(hObject, eventdata, handles)
% hObject    handle to slider_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
MPslider_str_Refresh(hObject, eventdata, handles);
MPslider_mode_Refresh(hObject, eventdata, handles);

MPupdate_mode_plot(hObject, eventdata, handles);
MPupdate_params_textbox(hObject, eventdata, handles);
MPupdate_eigenvals_textbox(hObject, eventdata, handles);


% --- MP structure slider refresh
% - reads NSTRuctures from data structure and current struct # from slider struct
% - rounds current structure value
% - writes current structure and text '1/N' to str slider and text box structures
%
function MPslider_str_Refresh(hObject, eventdata, handles)
NSTR = length(handles.modedata);
if (NSTR > 1),  val = round(get(handles.slider_str,'Value'));  else  val = 1;  end;
MPslider_set(handles.slider_str, val, [], [], handles.text_strnum);         % [MP]
if (val == 1),     st = 'off';  else  st = 'on';  end;  set(handles.pushbutton_strshiftleft, 'Enable', st);
if (val == NSTR),  st = 'off';  else  st = 'on';  end;  set(handles.pushbutton_strshiftright, 'Enable', st);
if (NSTR == 1),    st = 'off';  else  st = 'on';  end;  set(handles.pushbutton_strerase, 'Enable', st);


% --- Executes during object creation, after setting all properties.
function popupmenu_fieldcomp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_fieldcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu_fieldcomp.
function popupmenu_fieldcomp_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_fieldcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_fieldcomp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_fieldcomp
idx = get(handles.popupmenu_fieldcomp,'Value');
fieldnames = get(handles.popupmenu_fieldcomp,'String');
handles.fieldcomp = fieldnames{idx};
guidata(hObject, handles);
%fprintf('Running popup1 callback...%d, %s\n',idx,handles.fieldcomp);
MPupdate_mode_plot(hObject, eventdata, handles);                                  % [MP]


% --- Executes during object creation, after setting all properties.
function popupmenu_real_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_real (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu_real.
function popupmenu_real_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_real (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_real contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_real
idx = get(handles.popupmenu_real,'Value');
fieldfcnnames = get(handles.popupmenu_real,'String');
handles.fieldfunc = fieldfcnnames{idx};
guidata(hObject, handles);
MPupdate_mode_plot(hObject, eventdata, handles);                                  % [MP]


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
%set(hObject,'Erasemode','xor');    % [MP] faster update (default is 'normal', but 'xor' is also fast.
%set(hObject,'Doublebuffer','off');    % [MP] faster update (default is 'on', use 'on' for animation, 'off' for faster one-time rendering.


% [MP] - Custom functions are below:
%-- Update mode plot function --
function MPupdate_mode_plot(hObject, eventdata, handles)
SZ = handles.modedata;
%%%fprintf('Trying to get structure %f\n', get(handles.slider_str,'Value'));
nstr = round( get(handles.slider_str,'Value') );       % [MP] *** shouldn't need this, but when slider_str (structures) is *dragged* to a point it seems to end up non-integer! (problem)
nmode = get(handles.slider_mode,'Value');

% Get correct spatial coordinates
xcoord = getfield( SZ(nstr).F, handles.xcoord{get(handles.popupmenu_fieldcomp,'Value')} );  % Find correct coordinate set for chosen field component
ycoord = getfield( SZ(nstr).F, handles.ycoord{get(handles.popupmenu_fieldcomp,'Value')} );  % Find correct coordinate set for chosen field component
%fprintf('MPupdate - %s\n',handles.fieldcomp);
switch handles.fieldcomp
    case handles.fldlabels                                  % regular fields: {'Ex','Ey','Ez','Hx','Hy','Hz'}
        field = getfield( SZ(nstr).F, handles.fieldcomp );  % [MP] Choose Ex,y,z/Hx,y,z/Ste,Stm,Intensity?
        field = field(:,:,nmode);                           % [MP] Choose mode number
    case 'Ste'                                              % TE Poynting vector
        [Pz,SzTE] = ecrosshdotz(SZ(nstr).F, SZ(nstr).F, nmode*[1 1], 0, SZ(nstr).N.x, SZ(nstr).N.y);
        field = SzTE;
    case 'Stm'                                              % TM Poynting vector
        [Pz,SzTE,SzTM] = ecrosshdotz(SZ(nstr).F, SZ(nstr).F, nmode*[1 1], 0, SZ(nstr).N.x, SZ(nstr).N.y);
        field = SzTM;
    otherwise
        fprintf(['MODEVIEW: Plotting of quantity ' handles.fieldcomp ' not yet supported.\n']);
        field = zeros(length(xcoord),length(ycoord));
end
%field = feval( handles.fieldfunc, field(:,:,nmode) );    % [MP] Choose real/imag/abs
field = feval( handles.fieldfunc, field );    % [MP] Choose real/imag/abs
fieldmax = max(abs(field(:)));
if(fieldmax>0), field = field/fieldmax; end             % Normalize field

axes(handles.axes1);
if( get(handles.checkbox_grid,'Value')==1 )     % If non-uniform grid plotting...
    set(gcf,'renderer','zbuffer');
%    pcolorcenter(SZ(nstr).F.Rr/1e-6, SZ(nstr).F.Zr/1e-6, field.'); set(gca,'ydir','normal');
    pcolorcenter(xcoord, ycoord, field.'); set(gca,'ydir','normal');
    shading(gca,handles.shading);
else
%    imagesc(SZ(nstr).F.Rr/1e-6, SZ(nstr).F.Zr/1e-6, field.'); set(gca,'ydir','normal');
    imagesc(xcoord, ycoord, field.'); set(gca,'ydir','normal');
end
axis image;
% [MP] note (aug 15, '05): coordinates Rr, Zr need to be fixed to use correct ones for each field component.
colorbar;
colormap(handles.colormap);

zoom = str2double(get(handles.edit_zoom,'String'));  zoom = zoom(1);                    % [MP] update field zoom
if(zoom <= 0 || isnan(zoom)), zoom = 1; end
if( find(field < 0) ),  caxis([-1 1]/zoom);  else  caxis([0 1]/zoom);  end;


%-- Square resonator mode function for icon --
% v0.1 code from square_reson_model.m (Aug 15, 2005)
function F = MPsquare_reson_model(p,q,V,xova,yova)
if (nargin <= 3),  xova = (0:0.01:1.5); yova = xova;  end
[Y,X] = meshgrid(yova,xova);    % x along 1st dimension, y along 2nd dimension
F = zeros(size(X));
phix = mod(floor(p),2) * pi/2; phiy = mod(floor(q),2) * pi/2;     % phases: 0 for p,q even; pi/2 for p,q odd

ix = find(abs(xova)<=1); iy = find(abs(yova)<=1);
F(ix,iy) = cos(X(ix,iy) * p * pi/2 - phix) .* cos(Y(ix,iy) * q * pi/2 - phiy);

alphaxa = sqrt(V^2 - (p * pi/2)^2); alphaya = sqrt(V^2 - (q * pi/2)^2);
ix2 = find(abs(xova)>1); iy2 = find(abs(yova)>1);
F(ix,iy2) = cos(X(ix,iy2) * p * pi/2 - phix) .* (cos(1 * q * pi/2 - phiy) * exp(-alphaya * (Y(ix,iy2)-1)));
F(ix2,iy) = (cos(1 * p * pi/2 - phix) * exp(-alphaxa * (X(ix2,iy)-1))) .* cos(Y(ix2,iy) * q * pi/2 - phiy);
F(ix2,iy2) = 0;


%-- Fill textbox with parameters of simulated structure
function MPupdate_params_textbox(hObject, eventdata, handles)
nstr = round( get(handles.slider_str,'Value') );       % [MP] *** shouldn't need this, but when slider_str (structures) is *dragged* to a point it seems to end up non-integer! (problem)
SZ = handles.modedata(nstr);

if ~isfield(SZ,'P')
%    set(handles.edit_params,'Enable','off');
    strdata = {'No geometry parameters stored in current data structure.'};
else
%    set(handles.edit_params,'Enable','inactive');
    pfields = fieldnames(SZ.P);
    for k = 1:length(pfields)
        valdata = getfield(SZ.P,pfields{k});
        if isnumeric(valdata)
            if (length(valdata) == 1)                                   % Store single number in string
                strdata{k} = sprintf('%s:\t%g',pfields{k},valdata);
            else                                                        % Store vectors/arrays in string
                strdata{k} = sprintf('%s:\t[ ',pfields{k});
                for mm = 1:length(valdata)
                    strdata{k} = [strdata{k} sprintf('%g ',valdata(mm))];
                end
                strdata{k} = [strdata{k} ']'];
            end
        else                                                            % Bypass for non-numerical data
            strdata{k} = sprintf('%s:\t%g',pfields{k},'<Non-numeric data>');
        end
    end
end
set(handles.edit_params,'String',strdata);


%-- Fill textbox with parameters of simulated structure
function MPupdate_eigenvals_textbox(hObject, eventdata, handles)
nstr  = round( get(handles.slider_str,'Value') );       % [MP] *** shouldn't need this, but when slider_str (structures) is *dragged* to a point it seems to end up non-integer! (problem)
SZ    = handles.modedata(nstr);
nmode = round( get(handles.slider_mode,'Value') );       % [MP] *** shouldn't need this, but when slider_str (structures) is *dragged* to a point it seems to end up non-integer! (problem)

sfields = fieldnames(SZ);
%set(handles.edit_params,'Enable','inactive');
strdata = {};
for k = 1:length(sfields)
    switch lower(sfields{k})
        case 'wc'
            wc = getfield(SZ,sfields{k});
            if (length(wc) == 1), wc(nmode) = wc; end       % If only one frequency available, then set nth mode frequency to first (note: then it must be in beta mode)
            if ~isreal(wc(nmode))
                strdata{end+1} = sprintf('Resonance:\t%g nm', 2*pi*299792458/wc(nmode)/1e-9);
                strdata{end+1} = sprintf('Qexact:\t\t%g', -real(wc(nmode))/2/imag(wc(nmode)));  % [MP] assuming exp(i*k*x) convention
            else
                strdata{end+1} = sprintf('Wavelength:\t%g nm', 2*pi*299792458/wc(nmode)/1e-9);
            end
        case {'beta','betas'}
            beta = getfield(SZ,sfields{k});
            if (length(beta) == 1 && (nmode > 1)), beta(nmode) = NaN; end
            strdata{end+1} = sprintf('neff:\t\t%g', real(beta(nmode))/(2*pi/SZ.P.l0));
            ffields = fieldnames(SZ.F);
            if (isfield(SZ.F,'Ex') && isfield(SZ.F,'Ey') && isfield(SZ.F,'Ez') && isfield(SZ.F,'Hx') && isfield(SZ.F,'Hy') && isfield(SZ.F,'Hz') && ...
                isfield(SZ,'N') && exist('groupindex','file'))                 % Compute group index *if* have fields, dielectric (N) and groupindex function
                if ~isfield(SZ.F, 'dx'),  SZ.F.dx = diff(SZ.F.Rr(1:2));  end
                if ~isfield(SZ.F, 'dy'),  SZ.F.dy = diff(SZ.F.Zr(1:2));  end
                strdata{end+1} = sprintf('ngroup*:\t\t%g', real(feval('groupindex',SZ.F,SZ.N,nmode)) );    % Compute group index from fields/dielectric
            else
                strdata{end+1} = 'ngroup*:\t\t<missing fields, diel or groupindex function>';
            end
            if ~isreal(beta(nmode))
                strdata{end+1} = sprintf('Loss*:\t\t%g dB/cm', 0.2/log(10) * imag(beta(nmode)));    % Compute loss (dB/cm) from beta
            end
        case {'betabr','betar','gamma'}
            gamma = getfield(SZ,sfields{k});
            if (length(gamma) == 1 && (nmode > 1)), gamma(nmode) = NaN; end
            strdata{end+1} = sprintf('Angular propag.:\t%g rad/rad', real(gamma(nmode)) );
            ffields = fieldnames(SZ.F);
            if (isfield(SZ.F,'Ex') && isfield(SZ.F,'Ey') && isfield(SZ.F,'Ez') && isfield(SZ.F,'Hx') && isfield(SZ.F,'Hy') && isfield(SZ.F,'Hz') && ...
                isfield(SZ,'N') && exist('groupindex','file'))                 % Compute group index *if* have fields, dielectric (N) and groupindex function
                if ~isfield(SZ.F, 'dx'),  SZ.F.dx = diff(SZ.F.Rr(1:2));  end
                if ~isfield(SZ.F, 'dy'),  SZ.F.dy = diff(SZ.F.Zr(1:2));  end
                strdata{end+1} = sprintf('ngroup* (str. appx.):\t%g', real(feval('groupindex',SZ.F,SZ.N,nmode)) );    % Compute group index from fields/dielectric
            else
                strdata{end+1} = sprintf('ngroup*:\t\t<missing fields, diel or');
                strdata{end+1} = sprintf('\t\tgroupindex function>');
            end
            if ~isreal(gamma(nmode))
                strdata{end+1} = sprintf('Loss*:\t\t%g dB/90', 10*pi/log(10) * imag(gamma(nmode)));    % Compute loss (dB/90) from gamma
                strdata{end+1} = sprintf('Qapprox* (b/2a):\t%g', real(gamma(nmode))/2/imag(gamma(nmode)));    % Compute approx Q from gamma
            end
        otherwise
%            fprintf('Unrecognized structure field %s\n.', sfields{k});
    end
end
if isempty(strdata)
%    set(handles.edit_params,'Enable','off');
    strdata = {'No eigenvalue data stored in current data structure.'};
end
set(handles.edit_eigenvals,'String',strdata);


%-- Struct/mode erase dialog box: Erase structure(s) or particular mode(s) of structure(s):
function status = MPerasedialog_Create(hObject, eventdata, handles, defstrstr, defmodestr)

% Set up erase dialog box:
hPos = get(handles.figure_main,'Position');                 % Get position of main window...
hPos = round([hPos(1) + hPos(3)/2, hPos(2) + hPos(4)/2]);   % ...center point.
CBOK  = 'set(gcbf,''UserData'',''OK''); uiresume';
CBCan = 'set(gcbf,''UserData'',''Cancel''); uiresume';
h = dialog('Units','Pixels','Name','Delete structure(s)/mode(s):', ...
           'Position', ([hPos 0 0] + [250 100]*[-0.5 0 1 0; 0 -0.5 0 1]), ...
           'KeyPressFcn',['switch uint8(get(gcbo,''CurrentCharacter'')), case 13, ' CBOK '; case 27, ' CBCan '; end']);
h2 = uicontrol(h,'Style','text','HorizontalAlignment','left', ...
                 'Units','Characters','Position',[2 5 25 1.615], ...
                 'String','Structure(s):');
h3 = uicontrol(h,'Style','edit','BackgroundColor',[1 1 1], ...
                 'HorizontalAlignment','left','Units','Characters', ...
                 'Position', [15 5.2 32 1.615], 'String', defstrstr);
h4 = uicontrol(h,'Style','text','HorizontalAlignment','left', ...
                 'Units','Characters','Position',[2 3 25 1.615], ...
                 'String','Mode(s):');
h5 = uicontrol(h,'Style','edit','BackgroundColor',[1 1 1], ...
                 'HorizontalAlignment','left','Units','Characters', ...
                 'Position', [15 3.2 32 1.615], 'String', defmodestr);
h6 = uicontrol(h,'Style','pushbutton','String','Ok', 'FontWeight', 'Bold', ...
                 'Units','Characters', 'Position',[8 0.8 15 1.5], ...
                 'Callback', CBOK);
h7 = uicontrol(h,'Style','pushbutton','String','Cancel', ...
                 'Units','Characters', 'Position',[28 0.8 15 1.5], ...
                 'Callback', CBCan);
uiwait(h);  % Stop execution here till dialog is closed

switch(char(get(h, 'UserData')))    % Doesn't work if h is replaced by gcbf or gcbo, the latter two seem to be interpreted as the main window.
    case 'OK'
        istr = str2num(get(h3,'String')); imode = str2num(get(h5,'String'));    % Get list of structures and modes to delete
        if (~isempty(istr) && ~isempty(imode))
            kk = 1;     % Start at structure #1
            while (kk <= length(istr))
                NMODES = size(getfield(handles.modedata(istr(kk)), 'F', 'Ex'), 3);  % Number of modes in structure istr(kk)
                [imodkeep] = intersect((1:NMODES), setxor((1:NMODES), imode));        % Choose non-deleted modes to keep

                % Erase structure(s)/mode(s):
                if ~isempty(imodkeep)
                    sfields = fieldnames(handles.modedata(istr(kk)));   % Delete eigenvalue-related data
                    [~,ia,ib] = intersect(lower(sfields),lower(handles.veclabels));
                    sfields = sfields(ia);
                    for mm = 1:length(sfields)
                        tmp = getfield(handles.modedata(istr(kk)), sfields{mm});
                        if (length(tmp) >= max(imodkeep))
                            handles.modedata(istr(kk)) = setfield(handles.modedata(istr(kk)), sfields{mm}, tmp(imodkeep));
                            fprintf('Structure %d, field %s, keeping modes ', istr(kk), sfields{mm});
                            fprintf('%d ',imodkeep); fprintf('\n');
                        else
                            fprintf('Can''t reduce structure %d, field %s\n', istr(kk), sfields{mm});
                        end
                    end
                    sfields = fieldnames(handles.modedata(istr(kk)).F); % Delete field patterns
                    [~,ia,ib] = intersect(lower(sfields),lower(handles.fldlabels));
                    sfields = sfields(ia);
                    for mm = 1:length(sfields)
                        tmp = getfield(handles.modedata(istr(kk)).F, sfields{mm});
                        handles.modedata(istr(kk)).F = setfield(handles.modedata(istr(kk)).F, sfields{mm}, tmp(:,:,imodkeep));
                        fprintf('Structure %d, field F.%s, keeping modes ', istr(kk), sfields{mm});
                        fprintf('%d ',imodkeep); fprintf('\n');
                    end
                    
                    % --- slider update start
                    if (get(handles.slider_str, 'Value') == istr(kk))           % If cutting modes of current structure...
                        sliderval = round(get(handles.slider_mode,'Value'));    % ...update mode slider pointer

                        icurrmodeidx = find(imodkeep >= sliderval);             % Find new index to stay on "current" mode if it was not deleted.
                        if isempty(icurrmodeidx),  icurrmodeidx = 1;  end       % ...else, go to 1st mode.

                        fprintf('modesliderval = %d, modeslidermax = %d, istr(kk) = %d\n', sliderval, length(imodkeep), istr(kk));
                        MPslider_set(handles.slider_mode, icurrmodeidx(1), [], length(imodkeep), handles.text_modenum);
                    end
                    % --- slider update end
                    
                    kk = kk+1;  % Increment kk to next structure
                else
                    NSTR = length(handles.modedata);
                    if (NSTR > 1)
                        istrkeep = intersect((1:NSTR), setxor((1:NSTR), istr(kk)));     % Find structures to keep
                        handles.modedata = handles.modedata(istrkeep);
                        fprintf('All structure %d modes deleted. Erasing structure.\n', istr(kk));

                        % --- slider update start
                        sliderval = round(get(handles.slider_str,'Value'));     % Update structure slider pointer
                        slidermax = round(get(handles.slider_str,'Max'));
%                        fprintf('sliderval = %d, slidermax = %d, istr(kk) = %d\n', sliderval, slidermax, istr(kk));
                        if ((sliderval > istr(kk)) || (sliderval == slidermax))
                            MPslider_set(handles.slider_str, sliderval-1, [], [], []);  % No 'handles.text_strnum' because it's updated 2 lines below
                        end
                        MPslider_set(handles.slider_str, [], [], slidermax-1, handles.text_strnum);

                        MPslider_mode_Refresh(hObject, eventdata, handles);
                        % --- slider update end

                        istr = istr([1:(kk-1) (((kk+1):end)-1)]);       % Resort remaining structs to be deleted
                        % Don't increment kk (since a structure was deleted, the current kk is really the next structure now)
                    else
                        fprintf('Last remaining structure. Can''t delete.\n');
                        uiwait(errordlg('Cannot delete entire last structure.','Delete Error','modal'));
                        kk = kk+1;  % Increment kk to exit loop.
                    end
                end
                disp('Done erasing structs.');
            end
            fprintf('OK\n');
            status = 1;
        else
            uiwait(errordlg('No structure(s) or mode(s) selected.','Selection Error','modal'));
            status = 0;
        end
    case 'Cancel'
        fprintf('Cancel\n');
        status = 0;
end
guidata(hObject, handles);                  % Save modified data
fprintf('Gone out of MPstrmodeerase\n');
delete(h);


%-- Slider updating function --
function MPslider_set(hObject, value, minvalue, maxvalue, hObjectText)
if ~isempty(value),     set(hObject, 'Value', value);     end
if ~isempty(minvalue),  set(hObject, 'Min',   minvalue);  end
if ~isempty(maxvalue),
    set(hObject, 'Max', max(maxvalue,2));
    minvalue = get(hObject, 'Min');
    if (maxvalue > minvalue),  stepval = 1/(maxvalue-minvalue) * [1 1];  else  stepval = [1 10000];  end
    set(hObject, 'SliderStep', stepval);
end
if ~isempty(hObjectText)        % Update text box
    if isempty(maxvalue),  maxvalue = get(hObject,'Max');  end
    tmp = get(hObject,'SliderStep'); if(tmp(1) ~= tmp(2)),  maxvalue = 1;  end
    set(hObjectText, 'String', [num2str(get(hObject,'Value')) '/' num2str(maxvalue)]);
end


% --- Executes on button press in checkbox_diel.
function checkbox_diel_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_diel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_diel


% --- Executes during object creation, after setting all properties.
function edit_zoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function edit_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zoom as text
%        str2double(get(hObject,'String')) returns contents of edit_zoom as a double
MPupdate_mode_plot(hObject, eventdata, handles);


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_open_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
       {'*.mat','MAT-files (*.mat)'; '*.*',  'All Files (*.*)'}, ...
        'Open mode data file:');
%uiload;
if ~(isequal(filename,0) || isequal(pathname,0))     % If didn't press cancel
    indata = load(fullfile(pathname,filename), '-mat');
    inflds = fieldnames(indata);
    if ~isempty(inflds)
        if (length(inflds) > 1)     % Choose if more than 1 data structure variable in file
            [sel,ok] = listdlg('PromptString','Select file variable:', ...
                               'SelectionMode','single','Name','Import data structure', ...
                               'ListString',inflds);
        else
            sel = 1; ok = 1;
        end
        if(ok)
            handles.modedata = getfield(indata, inflds{sel});   % Import new data structure
            % Call to initialize new data structure
            guidata(hObject, handles);
            MPmodedata_initialize(hObject, eventdata, handles); % Reinitialize GUI w/ new data
        end
    else
        uiwait(errordlg('File is empty.','Open Error','modal'));
    end
    clear indata;   % Remove loaded file indata from memory
end


% --------------------------------------------------------------------
function menu_file_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrel = str2num(version('-release'));            % Choose file save format
filetypes = {'*.mat','MAT-files (*.mat)'};
    flags = {''};
if (mrel > 13)      % For versions after Matlab 6.x allow saving in v6 format
    filetypes(end+1,:) = {'*.mat','MATLAB Version 6.x compatible (*.mat)'};
    flags{end+1} = '-v6';
end
filetypes(end+1,:) = {'*.mat','MATLAB Version 4 compatible (*.mat)'};
flags{end+1} = '-v4';
    
[filename, pathname, filterindex] = uiputfile( filetypes, ...
        'Save mode data as data structure SZ in file:');
if ~(isequal(filename,0) || isequal(pathname,0))     % If didn't press cancel
    SZ = handles.modedata;
    save(fullfile(pathname,filename), 'SZ', flags{filterindex}); %, '-mat');  % -mat flag forgets to add .mat extension if not explicitly added
    clear SZ;           % Remove temporary save data from memory
%    fprintf(['ModeView: Saved data structure to (' fullfile(pathname,filename) ').\n']);
%    fprintf(['ModeView: Saved data structure to (' strrep(fullfile(pathname,filename), '\', '\\') ').\n']);  % [MP] Need \\ to show backslashes.
    disp(['ModeView: Saved data structure to (' fullfile(pathname,filename) ').']);     % [MP] Easier way to print \'s (cross-platform).
end


% --------------------------------------------------------------------
function menu_file_import_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vars = evalin('base','who');
if ~isempty(vars)
    [sel,ok] = listdlg('PromptString','Select workspace variable:', ...
                       'SelectionMode','single','Name','Import data structure', ...
                       'ListString',vars);
    if(ok)
        handles.modedata = evalin('base',vars{sel});        % Import new data structure
        % Call to initialize new data structure
        guidata(hObject, handles);
        MPmodedata_initialize(hObject, eventdata, handles); % Reinitialize GUI w/ new data
    end
else
    uiwait(errordlg('Base workspace is empty.','Import Error','modal'));
end


% --------------------------------------------------------------------
function menu_file_export_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
outvar = inputdlg('Output variable name:','Export',1);
if isempty(outvar{1})
    fprintf('Export: blank export variable - no action taken.\n');
else
    if ~isvarname(outvar{1})
        uiwait(errordlg({'Invalid MATLAB variable name.  A valid variable name is a ', ...
                         'character string of letters, digits and underscores, with ', ...
                         'length <= namelengthmax and the first character a letter.'},'Export error','modal'));
    else    % Export variable
        flag_var_exists = evalin('base',['exist(''' outvar{1} ''',''var'')']);  % Does variable already exist?
        flag_overwrite = 1;
        if(flag_var_exists)
            if( strcmp( questdlg('Variable exists in base workspace. Overwrite it?','Export','Yes','No','No'), 'No' ) )
                flag_overwrite = 0;
            end
        end
        if(flag_overwrite)
            assignin('base', outvar{1}, handles.modedata);      % Export the variable
            fprintf(['ModeView: structure exported to base workspace as ' outvar{1} '.\n']);
        end
    end
end


% --------------------------------------------------------------------
function menu_file_exit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure_main);    % [MP] close the gui


% --------------------------------------------------------------------
function menu_edit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_edit_delstructs_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_delstructs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton_modeerase_Callback(hObject, eventdata, handles);         % [MP] **Dirty call type - Matlab won't rename this line if I rename the pushbutton.


% --------------------------------------------------------------------
function menu_help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_help_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgtext = {['Modeview ' handles.releasever ' (' handles.releasedate ')'], ...
           '(C) 2005-2011 Milos Popovic', '', ...
           'Electromagnetic structure mode viewer for MATLAB for use with mode solver.', '', ...
           'Report bugs to <milos@mit.edu>.'};
%icondata=1:64; icondata=(icondata'*icondata)/64;
%load('squaremode'); icondata = real(F3);
F3 = MPsquare_reson_model(5.7, 3.8, 2*pi/1.55*1.8*sqrt(4^2-1^2), (0:1/32:1-1/32)+0.5/32, (0:1/32:1-1/32)+0.5/32); % Get icon
F3 = F3-F3';  icondata = [flipud([fliplr(-F3) F3]); fliplr(F3) -F3];                       % Make diagonal mode
icondata = round(64*(icondata - min(icondata(:)))/(max(icondata(:)) - min(icondata(:))));  % Normalize
msgbox(msgtext,'About .\\odeView','custom',icondata,jet(64),struct('WindowStyle','modal','Interpreter','tex'));
% TO DO: put in webpage for updates, free-software license note.


% --- Executes on button press in checkbox_pml.
function checkbox_pml_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_pml (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_pml


% --- Executes on button press in pushbutton_strshiftleft.
function pushbutton_strshiftleft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_strshiftleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
istr = round(get(handles.slider_str, 'Value'));                         % Get struct #; should be fine for 1 struct long cell-array, since callback always resets # to 1 then
handles.modedata = handles.modedata([1:istr-2 istr istr-1 istr+1:end]); % Right shift current structure
guidata(hObject,handles);                                               % Store new data
MPslider_set(handles.slider_str, istr-1, [], [], []); %, handles.text_strnum);
MPslider_str_Refresh(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_strshiftright.
function pushbutton_strshiftright_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_strshiftright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
istr = round(get(handles.slider_str, 'Value'));                         % Get struct #; should be fine for 1 struct long cell-array, since callback always resets # to 1 then
handles.modedata = handles.modedata([1:istr-1 istr+1 istr istr+2:end]); % Right shift current structure
guidata(hObject,handles);                                               % Store new data
MPslider_set(handles.slider_str, istr+1, [], [], []); %, handles.text_strnum);
MPslider_str_Refresh(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_modeshiftleft.
function pushbutton_modeshiftleft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_modeshiftleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
istr  = round(get(handles.slider_str,  'Value'));                       % Get struct #; should be fine for 1 struct long cell-array, since callback always resets # to 1 then
imode = round(get(handles.slider_mode, 'Value'));                       % Get struct #; should be fine for 1 mode long cell-array, since callback always resets # to 1 then
NMODES= round(get(handles.slider_mode, 'Max'));                         % Get struct #; should be fine for 1 mode long cell-array, since callback always resets # to 1 then
iord = [1:imode-2 imode imode-1 imode+1:NMODES];                        % Index vector for right shifting current mode

sfields = fieldnames(handles.modedata(istr));                           % Reorder eigenvalue-related data
[tmp,ia,ib] = intersect(lower(sfields),lower(handles.veclabels));  sfields = sfields(ia);   % Find relevant datfields present in structure
for mm = 1:length(sfields)
    tmp = getfield(handles.modedata(istr), sfields{mm});
    if (length(tmp) >= max(iord))
        handles.modedata(istr) = setfield(handles.modedata(istr), sfields{mm}, tmp(iord));
        fprintf('Structure %d, field %s, keeping modes ', istr, sfields{mm});
        fprintf('%d ',iord); fprintf('\n');
    else    % e.g. when wc has length 1, it can't be reordered..:
        fprintf('Can''t reorder structure %d, field %s\n', istr, sfields{mm});
    end
end
sfields = fieldnames(handles.modedata(istr).F);                         % Reorder field patterns
[tmp,ia,ib] = intersect(lower(sfields),lower(handles.fldlabels));  sfields = sfields(ia);
for mm = 1:length(sfields)
    tmp = getfield(handles.modedata(istr).F, sfields{mm});
    handles.modedata(istr).F = setfield(handles.modedata(istr).F, sfields{mm}, tmp(:,:,iord));
    fprintf('Structure %d, field F.%s, keeping modes ', istr, sfields{mm});
    fprintf('%d ',iord); fprintf('\n');
end

guidata(hObject,handles);                                               % Store new data
MPslider_set(handles.slider_mode, imode-1, [], [], []); %, handles.text_modenum);
MPslider_mode_Refresh(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_modeshiftright.
function pushbutton_modeshiftright_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_modeshiftright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
istr  = round(get(handles.slider_str,  'Value'));                       % Get struct #; should be fine for 1 struct long cell-array, since callback always resets # to 1 then
imode = round(get(handles.slider_mode, 'Value'));                       % Get struct #; should be fine for 1 mode long cell-array, since callback always resets # to 1 then
NMODES= round(get(handles.slider_mode, 'Max'));                         % Get struct #; should be fine for 1 mode long cell-array, since callback always resets # to 1 then
iord = [1:imode-1 imode+1 imode imode+2:NMODES];                        % Index vector for right shifting current mode

sfields = fieldnames(handles.modedata(istr));                           % Reorder eigenvalue-related data
[tmp,ia,ib] = intersect(lower(sfields),lower(handles.veclabels));  sfields = sfields(ia);   % Find relevant datfields present in structure
for mm = 1:length(sfields)
    tmp = getfield(handles.modedata(istr), sfields{mm});
    if (length(tmp) >= max(iord))
        handles.modedata(istr) = setfield(handles.modedata(istr), sfields{mm}, tmp(iord));
        fprintf('Structure %d, field %s, keeping modes ', istr, sfields{mm});
        fprintf('%d ',iord); fprintf('\n');
    else    % e.g. when wc has length 1, it can't be reordered..:
        fprintf('Can''t reorder structure %d, field %s\n', istr, sfields{mm});
    end
end
sfields = fieldnames(handles.modedata(istr).F);                         % Reorder field patterns
[tmp,ia,ib] = intersect(lower(sfields),lower(handles.fldlabels));  sfields = sfields(ia);
for mm = 1:length(sfields)
    tmp = getfield(handles.modedata(istr).F, sfields{mm});
    handles.modedata(istr).F = setfield(handles.modedata(istr).F, sfields{mm}, tmp(:,:,iord));
    fprintf('Structure %d, field F.%s, keeping modes ', istr, sfields{mm});
    fprintf('%d ',iord); fprintf('\n');
end

guidata(hObject,handles);                                               % Store new data
MPslider_set(handles.slider_mode, imode+1, [], [], []); %, handles.text_modenum);
MPslider_mode_Refresh(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_strerase.
function pushbutton_strerase_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_strerase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nstr   = round( get(handles.slider_str, 'Value') );     % [MP] *** shouldn't need this, but when slider_str (structures) is *dragged* to a point it seems to end up non-integer! (problem)
NMODES = round( get(handles.slider_mode,'Max')   );     % [MP] *** shouldn't need this, but when slider_str (structures) is *dragged* to a point it seems to end up non-integer! (problem)

if( MPerasedialog_Create(handles.figure_main, eventdata, handles, num2str(nstr), ['[1:' num2str(NMODES) ']']) )
    handles = guidata(hObject);     % Refresh handles, modified by MPerase.. for calls below
    fprintf('Done erasing struct parts.. current str %d mode %d\n', get(handles.slider_str,'Value'), get(handles.slider_mode,'Value'));

    % Refresh sliders and plots
    MPslider_str_Refresh(hObject, eventdata, handles);          % Second refresh needed (although sliders already updated in MPerase.., to update erase buttons 'enabled' property)
    MPslider_mode_Refresh(hObject, eventdata, handles);
    MPupdate_mode_plot(hObject, eventdata, handles);                                  % [MP]
    MPupdate_params_textbox(hObject, eventdata, handles);
    MPupdate_eigenvals_textbox(hObject, eventdata, handles);
end


% --- Executes on button press in pushbutton_modeerase.
function pushbutton_modeerase_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_modeerase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
istr  = round( get(handles.slider_str, 'Value') );     % [MP] *** shouldn't need this, but when slider_str (structures) is *dragged* to a point it seems to end up non-integer! (problem)
imode = round( get(handles.slider_mode,'Value') );     % [MP] *** shouldn't need this, but when slider_str (structures) is *dragged* to a point it seems to end up non-integer! (problem)

if( MPerasedialog_Create(handles.figure_main, eventdata, handles, num2str(istr), num2str(imode)) )
    handles = guidata(hObject);     % Refresh handles, modified by MPerase.. for calls below
    fprintf('Done erasing modes..\n');

    % Refresh sliders and plots
    MPslider_str_Refresh(hObject, eventdata, handles);          % Second refresh needed (although sliders already updated in MPerase.., to update erase buttons 'enabled' property)
    MPslider_mode_Refresh(hObject, eventdata, handles);
    MPupdate_mode_plot(hObject, eventdata, handles);                                  % [MP]
    MPupdate_params_textbox(hObject, eventdata, handles);
    MPupdate_eigenvals_textbox(hObject, eventdata, handles);
end


% --- Executes during object creation, after setting all properties.
function edit_params_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function edit_eigenvals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_eigenvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function edit_eigenvals_Callback(hObject, eventdata, handles)
% hObject    handle to edit_eigenvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_eigenvals as text
%        str2double(get(hObject,'String')) returns contents of edit_eigenvals as a double



% --- Executes on button press in checkbox_grid.
function checkbox_grid_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_grid
MPpopupmenu_shading_Refresh(hObject, eventdata, handles);  handles = guidata(hObject);
MPupdate_mode_plot(hObject, eventdata, handles);                                  % [MP]


% --- Executes during object creation, after setting all properties.
function popupmenu_shading_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_shading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu_shading.
function popupmenu_shading_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_shading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_shading contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_shading
MPpopupmenu_shading_Refresh(hObject, eventdata, handles);  handles = guidata(hObject);
MPupdate_mode_plot(hObject, eventdata, handles);                                  % [MP]


% MP: refresh status/value of shading box
function MPpopupmenu_shading_Refresh(hObject, eventdata, handles)
contents = get(handles.popupmenu_shading,'String');                   % Store shading choice in variable
handles.shading = contents{get(handles.popupmenu_shading,'Value')};
contents = get(handles.popupmenu_colormap,'String');                  % Store colormap choice in variable
handles.colormap = eval(contents{get(handles.popupmenu_colormap,'Value')});
guidata(hObject, handles);
if( get(handles.checkbox_grid,'Value') == 0 ),  stat='off';  else  stat='on';  end;     % Do or don't gray out shading type popupmenu
set(handles.popupmenu_shading,'Enable',stat); set(handles.text_shading,'Enable',stat);  % Set gray-out of shading menu and title text


% --- Executes during object creation, after setting all properties.
function slider_zoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function slider_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to slider_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Check zoom direction
if (get(handles.slider_zoom,'Value') == 1),  zmdir = -1;  else  zmdir = +1;  end

% Update zoom factor
%fprintf('Zoom pushed, %d\n', get(hObject,'Value'));
zoomval = str2num(get(handles.edit_zoom, 'String'));
%zoomval = 10.^(round(3*log10([1 2 5 10]))/3)
zoomval = 10.^( ( round(3*log10(zoomval))+zmdir )/3 );   % Make it round to 1, 2, 5, 10, 20, 50, 100...
dec = floor(log10(zoomval));  zoomval = round(zoomval/10.^dec)*10.^dec;   % Round to first sig fig.
zoomtxt = sprintf('%1.5g',zoomval);
if ispc, zoomtxt = strrep(zoomtxt, 'e+0', 'e+');   end
set(handles.edit_zoom, 'String', zoomtxt);
set(hObject, 'Value', 2);

% Reset zoom slider to 2
set(handles.slider_zoom,'Value',2);   % Set back to middle for next click

MPupdate_mode_plot(hObject, eventdata, handles);


% --- Executes on selection change in popupmenu_colormap.
function popupmenu_colormap_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_colormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_colormap
MPpopupmenu_shading_Refresh(hObject, eventdata, handles);  handles = guidata(hObject);
MPupdate_mode_plot(hObject, eventdata, handles);                                  % [MP]

% --- Executes during object creation, after setting all properties.
function popupmenu_colormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
