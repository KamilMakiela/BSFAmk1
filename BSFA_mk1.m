

% char(cell_object) powoduje opuczszenie cudzyslowow stringa
% cellstr(string_object) powoduje dodanie cudzyslowow do stringa
function varargout = BSFA_mk1(varargin)
% BSFA_MK1 M-file for BSFA_mk1.fig
%      BSFA_MK1, by itself, creates a new BSFA_MK1 or raises the existing
%      singleton*.
%
%      H = BSFA_MK1 returns the handle to a new BSFA_MK1 or the handle to
%      the existing singleton*.
%
%      BSFA_MK1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BSFA_MK1.M with the given input arguments.
%
%      BSFA_MK1('Property','Value',...) creates a new BSFA_MK1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BSFA_mk1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BSFA_mk1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BSFA_mk1

% Last Modified by GUIDE v2.5 15-Jul-2017 15:52:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @BSFA_mk1_OpeningFcn, ...
    'gui_OutputFcn',  @BSFA_mk1_OutputFcn, ...
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

% --- Executes just before BSFA_mk1 is made visible.
function BSFA_mk1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BSFA_mk1 (see VARARGIN)
% Choose default command line output for BSFA_mk1

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes BSFA_mk1 wait for user response (see UIRESUME)
% uiwait(handles.BSFA_mk1);

% --- Outputs from this function are returned to the command line.
function varargout = BSFA_mk1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Import_Y.
function Import_Y_Callback(hObject, eventdata, handles)
% hObject    handle to Import_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.user.n ~= 0 && handles.user.t ~= 0)
    [handles.user.ylabel, handles.user.ydata, handles.user.dmulabel, handles.user.timelabel] = pobierz(1, handles.user.n, handles.user.t);
    handles.user.y = log(handles.user.ydata);
    handles.user.yimport = 1;
else
    msgbox('You did not provide the number of DMUs and/or the number of time periods', 'Well this is embarrassing... :)', 'error');
end
guidata(hObject,handles);

% --- Executes on button press in Import_X.
function Import_X_Callback(hObject, eventdata, handles)
% hObject    handle to Import_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.user.n ~= 0 && handles.user.t ~= 0 && handles.user.factors ~=0)
    [handles.user.xlabel, handles.user.xdata, dmu_label, time_label] = pobierz(handles.user.factors, handles.user.n, handles.user.t);
    handles.user.xRaw = log(handles.user.xdata);
    handles.user.ximport = 1;
else
    msgbox('You did not provide the number of DMUs and/or the number of time periods', 'Well this is embarrassing... :)', 'error');
end
guidata(hObject,handles);


% --- Executes on button press in Policz.
function Policz_Callback(hObject, eventdata, handles)
% hObject    handle to Policz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global km_handles
global W2
handles.user.W = W2;
if handles.user.ModelType == 7 || handles.user.ModelType == 8
    h = msgbox('This model type has been removed from this version', 'Well this is embarrassing...', 'warn');
    handles.user.OutputFolder = 'lipa';
    waitfor(h);
elseif ~(handles.user.n == 0 || handles.user.t == 0 || handles.user.factors == 0 || handles.user.TotalCycles == 0 || handles.user.BurninCycles == 0 || handles.user.yimport == 0 || handles.user.ximport == 0 || handles.user.TotalCycles <= handles.user.BurninCycles)
    %appropriate data preparation
    [handles.user.COLS.elast, handles.user.COLS.beta, handles.user.COLS.u, handles.user.COLS.sse, handles.user.x, handles.user.xlabel, handles.user.vecs.all, handles.user.vecs.obj, handles.user.vecs.time, handles.user.restrictions, handles.user.OLSresiduals, handles.user.priors] ...
        = PrepareData(handles.user.y, handles.user.xRaw, handles.user.n, handles.user.t, handles.user.xlabel, handles.user.ModelFunction, handles.user.RestrOptions);
    
    %SIMULATION
    %%%%%%%%%%%%%%
    %%%% For PRODUCTION models
    [handles.user.OutputFolder, handles.user.OutputFiles] = Simulation(handles.user.y, handles.user.x, handles.user.n, handles.user.t, handles.user.COLS.u, handles.user.COLS.beta, handles.user.TotalCycles, handles.user.ModelFunction, km_handles.H, handles.user.MaxChainSize, handles.user.BurninCycles, handles.user.restrictions, handles.user.ModelType, handles.user.med_apr, handles.user.priors);
    %%%% For COST model
    %handles.user.y = - handles.user.y; handles.user.x = - handles.user.x; 
    %[handles.user.COLS.elast, handles.user.COLS.beta, handles.user.COLS.u, handles.user.COLS.sse, handles.user.OLSresiduals] = licz_COLS(handles.user.y, handles.user.x, 0);
    %[handles.user.OutputFolder, handles.user.OutputFiles] = Simulation(handles.user.y, handles.user.x, handles.user.n, handles.user.t, handles.user.COLS.u, handles.user.COLS.beta, handles.user.TotalCycles, handles.user.ModelFunction, km_handles.H, handles.user.MaxChainSize, handles.user.BurninCycles, handles.user.restrictions, handles.user.ModelType, handles.user.med_apr, handles.user.priors);
    %%%%%%%%%%%%%%
    
    handles.user.ready2decom = 1;
    %appropriate analysis
    if ~strcmp(handles.user.OutputFolder,'lipa')
        [handles.user.BSFA.beta, handles.user.BSFA.elast, handles.user.BSFA.ElastScale, handles.user.BSFA.eff, handles.user.BSFA.LambdaOrOmega2, handles.user.BSFA.s, handles.user.BSFA.mdd, handles.user.BSFA.SSE, handles.user.CyclesUsedPerFile] ...
            = stats(handles.user.OutputFolder, handles.user.BurninCycles, handles.user.ModelFunction, km_handles.H, handles.user.x, handles.user.y, handles.user.vecs, handles.user.OutputFiles, handles.user.n, handles.user.t);
    end
else
    h = msgbox('Everything is input correctly and try again', 'Well this is embarrassing...', 'warn');
    handles.user.OutputFolder = 'lipa';
    waitfor(h);
end
fprintf('folder with results is: %s\n', handles.user.OutputFolder);
assignin('base', 'OutputResults', handles.user);
user = handles.user;
save(strcat(handles.user.OutputFolder,'\','results'),'user');
guidata(hObject,handles);

% --- Executes on selection change in TypModelu.
function TypModelu_Callback(hObject, eventdata, handles)
% hObject    handle to TypModelu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TypModelu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        TypModelu
handles.user.ModelFunction = get(hObject,'Value');
%disp(handles.user.ModelFunction);

% do aktualizacji
switch handles.user.ModelFunction
    case 3
        msgbox('This model is not available in this version'); 
        handles.user.ModelFunction = 1;
    case {2, 4, 7, 8, 9, 10, 11, 12}
        if handles.user.t == 1
            msgbox('You chose a model with trend but the number of observations over time is 1; select a proper model', 'Well this is embarrassing', 'error');
            handles.user.ModelFunction = 1;
        end
    otherwise
        disp('Model choice Ok');
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function TypModelu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TypModelu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dmu_Callback(hObject, eventdata, handles)
% hObject    handle to dmu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dmu as text
%        str2double(get(hObject,'String')) returns contents of dmu as a double

handles.user.n = str2double(get(hObject,'String'));
%disp(handles.user.n);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function dmu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dmu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function time_Callback(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time as text
%        str2double(get(hObject,'String')) returns contents of time as a double
handles.user.t = str2double(get(hObject,'String'));
%disp(handles.user.t);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function Projekt_PSk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BSFA_mk1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global km_handles
km_handles.H = gca; %axis handler in GUI!
%disp(km_handles.H);

handles.output = hObject;
handles.user.ModelType = 1;
handles.user.ModelFunction = 1;
handles.user.n = 0;
handles.user.t = 0;
handles.user.factors = 0;
handles.user.med_apr = 0.75;
handles.user.TotalCycles = 50000;
handles.user.BurninCycles = 10000;
handles.user.MaxChainSize = 100000;
handles.user.OLSresiduals = 0;
handles.user.RestrOptions.case = 1;
handles.user.RestrOptions.m = 0;
handles.user.RestrOptions.zm = 0;
handles.user.yimport = 0;
handles.user.ximport = 0;
handles.user.ready2decom = 0;
%handles.user.shift = 0;
assignin('base', 'InitialParameters', handles.user);

guidata(hObject, handles);

function n_var_Callback(hObject, eventdata, handles)
% hObject    handle to n_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_var as text
%        str2double(get(hObject,'String')) returns contents of n_var as a double

handles.user.factors = str2double(get(hObject,'String'));
%disp(handles.user.factors);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function n_var_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function iteracji_Callback(hObject, eventdata, handles)
% hObject    handle to iteracji (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iteracji as text
%        str2double(get(hObject,'String')) returns contents of iteracji as
%        a double

handles.user.TotalCycles = str2double(get(hObject,'String'));
%disp(handles.user.TotalCycles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function iteracji_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iteracji (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spalonych_Callback(hObject, eventdata, handles)
% hObject    handle to spalonych (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spalonych as text
%        str2double(get(hObject,'String')) returns contents of spalonych as a double
handles.user.BurninCycles = str2double(get(hObject,'String'));
%disp(handles.user.BurninCycles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function spalonych_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spalonych (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function objpliku_Callback(hObject, eventdata, handles)
% hObject    handle to objpliku (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of objpliku as text
%        str2double(get(hObject,'String')) returns contents of objpliku as a double
handles.user.MaxChainSize = str2double(get(hObject,'String'));
%disp(handles.user.MaxChainSize);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function objpliku_CreateFcn(hObject, eventdata, handles)
% hObject    handle to objpliku (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in gen_reset.
function gen_reset_Callback(hObject, eventdata, handles)
% hObject    handle to gen_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rng(0);
disp(rng);

% --- Executes on button press in gen_shuffle.
function gen_shuffle_Callback(hObject, eventdata, handles)
% hObject    handle to gen_shuffle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rng('shuffle');
disp(rng);

% --- Executes on selection change in typ_generatora.
function typ_generatora_Callback(hObject, eventdata, handles)
% hObject    handle to typ_generatora (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns typ_generatora contents as cell array
%        contents{get(hObject,'Value')} returns selected item from typ_generatora
% !!UWAGA funkcja domyslnie ustawie SEED na ZERO!!
opcja = get(hObject,'Value');
switch opcja
    case 1
        handles.user.generator = 'twister';
    case 2
        handles.user.generator = 'combRecursive';
    case 3
        handles.user.generator = 'multFibonacci';
    case 4
        handles.user.generator = 'v5uniform';
    case 5
        handles.user.generator = 'v5normal';
    case 6
        handles.user.generator = 'v4';
    otherwise
        disp('Something is wrong with the generator choice');
end
rng(0,handles.user.generator);
disp(rng);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function typ_generatora_CreateFcn(hObject, eventdata, handles)
% hObject    handle to typ_generatora (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Sposob_restr.
function Sposob_restr_Callback(hObject, eventdata, handles)
% hObject    handle to Sposob_restr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Sposob_restr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Sposob_restr
handles.user.RestrOptions.case = get(hObject,'Value');
disp(handles.user.RestrOptions);
guidata(hObject, handles);
switch handles.user.RestrOptions.case
    case {2, 6, 8}
        msgbox('Provide the number of DMUs to be accounted in restrictions. Enter this value below.','Remember','help');
    case 3
        msgbox('Enter the variable number you wish to impose restrictions on (the number is based on the order you selected inputs). Enter this value below. ','Remember','help');
    case 4
        msgbox('Provide the number of DMUs to be accounted in restrictions and the variable number you wish to impose restrictions on (the number is based on the order you selected inputs). Enter this value below. ','Remember','help');
    otherwise
        
end

% --- Executes during object creation, after setting all properties.
function Sposob_restr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sposob_restr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Liczba_m_Callback(hObject, eventdata, handles)
% hObject    handle to Liczba_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Liczba_m as text
%        str2double(get(hObject,'String')) returns contents of Liczba_m as a double
handles.user.RestrOptions.m = str2double(get(hObject,'String'));
disp(handles.user.RestrOptions);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Liczba_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Liczba_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function res_zmienna_Callback(hObject, eventdata, handles)
% hObject    handle to res_zmienna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of res_zmienna as text
%        str2double(get(hObject,'String')) returns contents of res_zmienna as a double
handles.user.RestrOptions.zm = str2double(get(hObject,'String'));
disp(handles.user.RestrOptions);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function res_zmienna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to res_zmienna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dekompozycja.
function dekompozycja_Callback(hObject, eventdata, handles)
% hObject    handle to dekompozycja (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.user.ready2decom == 1
    choice = questdlg(...
        {'There are two options for decomposing output growth rates' ...
        'Makiela (2014) - decomposes productivity change into technical and efficiency change' ...
        'Makiela, K., Marzec J., Pisulewski, A. (2016) - decomposes productivity change into technical, scale and efficiency change' ...
        }, ...
        'Decomposition menu', ...
        'Makiela (2014)', 'Makiela et al. (2015)', 'Makiela (2014)');
    
    switch choice
        case 'Makiela (2014)'
            handles.user.decomp = decompose(handles.user.x,handles.user.n,handles.user.t,handles.user.OutputFolder,handles.user.OutputFiles,handles.user.CyclesUsedPerFile,handles.user.ModelFunction,handles.user.factors,handles.user.vecs,handles.user.xRaw);
        case 'Makiela et al. (2017)'
            handles.user.decomp = decompose2(handles.user.x,handles.user.n,handles.user.t,handles.user.OutputFolder,handles.user.OutputFiles,handles.user.CyclesUsedPerFile,handles.user.ModelFunction,handles.user.factors,handles.user.vecs,handles.user.xRaw);
    end
    assignin('base', 'OutputResults', handles.user);
    disp('Decomposition results have been added to OutputResults structure. See "decomp" field');
else
    msgbox('No model has been estimated. If you wish to decompose previous a estimation use the function in decompose.m directly from the command line. ','No data for decomposition','error');
end

guidata(hObject, handles);


% --- Executes on button press in Test_reszt.
function Test_reszt_Callback(hObject, eventdata, handles)
% hObject    handle to Test_reszt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~(handles.user.n == 0 || handles.user.t == 0 || handles.user.factors == 0 || handles.user.yimport == 0 || handles.user.ximport == 0)
    [handles.user.COLS.elast, handles.user.COLS.beta, handles.user.COLS.u, handles.user.COLS.sse, handles.user.x, handles.user.xlabel, handles.user.vecs.all, handles.user.vecs.obj, handles.user.vecs.time, handles.user.restrictions, handles.user.OLSresiduals] ...
        = PrepareData(handles.user.y, handles.user.xRaw, handles.user.n, handles.user.t, handles.user.xlabel, handles.user.ModelFunction, handles.user.RestrOptions);
    skwns = skewness(handles.user.OLSresiduals);
    
    %assignin('base', 'reszty', handles.user.OLSresiduals);
    assignin('base', 'OutputResults', handles.user);
    if skwns < 0
        msg = sprintf('Distribution of OLS residuals has a negative asymmetry: %4.4f.\nThis reassures us about the existence of inefficiency in the data\nEverything seems ok.',skwns);
        msgbox(msg);
    else
        msg = sprintf('Warning! Distribution of OLS residuals has a positive asymmetry: %4.4f. \nSomething may be wrong with the data or the model you chose. ',skwns);
        msgbox(msg);
        %     data.stdev = std(handles.user.OLSresiduals);
        %     data.mom2 = moment(handles.user.OLSresiduals,2);
        %     data.mom3 = moment(handles.user.OLSresiduals,3);
        %assignin('base', 'wyszlo', skwns);
    end
else
    msgbox('To run this test first you need to import data and choose your model; see steps 1-6 in the manual', 'Well this is embarrassing...', 'warn');
end


% --- Executes on selection change in model_kl.
function model_kl_Callback(hObject, eventdata, handles)
% hObject    handle to model_kl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_kl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_kl
handles.user.ModelType = get(hObject,'Value');
disp(handles.user.ModelType);


if handles.user.ModelType == 3 || handles.user.ModelType == 4
    ile_pobran = handles.user.ModelType - 2;
    global W2;
    [~, W2, ~, ~] = pobierz(ile_pobran, handles.user.n, handles.user.t);
    
    %disp(W2);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function model_kl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_kl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function mediana_apr_Callback(hObject, eventdata, handles)
% hObject    handle to mediana_apr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mediana_apr as text
%        str2double(get(hObject,'String')) returns contents of mediana_apr as a double
mediana = str2double(get(hObject,'String'));
if mediana < 0.5 || mediana > 0.98
    msgbox('The value should be set between 0.5 and 0.98','Wrong value','error');
else
    disp(mediana);
    handles.user.med_apr = mediana;
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function mediana_apr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mediana_apr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
