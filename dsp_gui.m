function varargout = dsp_gui(varargin)
% DSP_GUI MATLAB code for dsp_gui.fig
%      DSP_GUI, by itself, creates a new DSP_GUI or raises the existing
%      singleton*.
%
%      H = DSP_GUI returns the handle to a new DSP_GUI or the handle to
%      the existing singleton*.
%
%      DSP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSP_GUI.M with the given input arguments.
%
%      DSP_GUI('Property','Value',...) creates a new DSP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dsp_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dsp_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dsp_gui

% Last Modified by GUIDE v2.5 16-Jun-2020 11:16:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dsp_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @dsp_gui_OutputFcn, ...
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


% --- Executes just before dsp_gui is made visible.
function dsp_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no out_graph args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dsp_gui (see VARARGIN)

set(handles.record_panel,'visible','off');          %Hide the record panel

% Choose default command line out_graph for dsp_gui
handles.op = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dsp_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dsp_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning out_graph args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line out_graph from handles structure
varargout{1} = handles.op;



% --- Executes on button press in all_graph_back.
function all_graph_back_Callback(hObject, eventdata, handles)
% hObject    handle to all_graph_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.output_panel,'visible','on');       %Make the output panel visible

% --- Executes on button press in sound_play.
function sound_play_Callback(hObject, eventdata, handles)
% hObject    handle to sound_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sound(handles.data, 48000);         %play the original sound
guidata(hObject, handles);

% --- Executes on button press in sound_view.
function sound_view_Callback(hObject, eventdata, handles)
% hObject    handle to sound_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(4);
plot(handles.data,'b','LineWidth',2);       %plot the original sound
title('Original Sound');
xlabel('time');
ylabel('Amplitude');

% --- Executes on button press in voiced_play.
function voiced_play_Callback(hObject, eventdata, handles)
% hObject    handle to voiced_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sound(handles.voiced, 48000);           %play the voiced signal
guidata(hObject, handles);

% --- Executes on button press in voiced_view.
function voiced_view_Callback(hObject, eventdata, handles)
% hObject    handle to voiced_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(5);
plot(handles.voiced,'r','LineWidth',2);     %plot the voiced signal
title('Voiced before filtering');
xlabel('time');
ylabel('Amplitude');

% --- Executes on button press in out_play.
function out_play_Callback(hObject, eventdata, handles)
% hObject    handle to out_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sound(handles.op, 48000);                   %play the output signal
guidata(hObject, handles);


% --- Executes on button press in out_view.
function out_view_Callback(hObject, eventdata, handles)
% hObject    handle to out_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(6);
plot(handles.op,'y','LineWidth',2);             %plot the output signal
title('Output');
xlabel('time');
ylabel('Amplitude');


% --- Executes on button press in zcr_play.
function zcr_play_Callback(hObject, eventdata, handles)
% hObject    handle to zcr_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sound(handles.zcr, 48000);                  %play the small ZCR signal
guidata(hObject, handles);

% --- Executes on button press in zcr_view.
function zcr_view_Callback(hObject, eventdata, handles)
% hObject    handle to zcr_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(7);
plot(handles.zcr,'g','LineWidth',2);         %plot the small ZCR signal
title('Small ZCR Signal');
xlabel('time');
ylabel('Amplitude');

% --- Executes on button press in energy_play.
function energy_play_Callback(hObject, eventdata, handles)
% hObject    handle to energy_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sound(handles.energy, 48000);       %play the energy signal
guidata(hObject, handles);

% --- Executes on button press in energy_view.
function energy_view_Callback(hObject, eventdata, handles)
% hObject    handle to energy_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(8);
plot(handles.energy,'g','LineWidth',2);     %plot the energy signal
title('High Energy Signal');
xlabel('time');
ylabel('Amplitude');

% --- Executes on button press in unvoiced_play.
function unvoiced_play_Callback(hObject, eventdata, handles)
% hObject    handle to unvoiced_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sound(handles.unvoiced, 48000);         %play the unvoiced signal
guidata(hObject, handles);

% --- Executes on button press in unvoiced_view.
function unvoiced_view_Callback(hObject, eventdata, handles)
% hObject    handle to unvoiced_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(9);
plot(handles.unvoiced,'r','LineWidth',2);  %plot the unvoiced signal
title('Unvoiced');
xlabel('time');
ylabel('Amplitude');


% --- Executes on button press in output_back.
function output_back_Callback(hObject, eventdata, handles)
% hObject    handle to output_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.start_panel,'visible','on');         %Make the start panel visible

% --- Executes on button press in original_play.
function original_play_Callback(hObject, eventdata, handles)
% hObject    handle to original_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sound(handles.data, 48000);         %play the original signal
guidata(hObject, handles);


% --- Executes on button press in original_view.
function original_view_Callback(hObject, eventdata, handles)
% hObject    handle to original_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(2);
plot(handles.data,'b','LineWidth',2);       %plot the original signal
title('Original Sound');
xlabel('time');
ylabel('Amplitude');

% --- Executes on button press in output_play.
function output_play_Callback(hObject, eventdata, handles)
% hObject    handle to output_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sound(handles.op, 48000);                   %play the output signal
guidata(hObject, handles);

% --- Executes on button press in output_view.
function output_view_Callback(hObject, eventdata, handles)
% hObject    handle to output_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(3);
plot(handles.op,'y','LineWidth',2);     %plot the output signal
title('Output');
xlabel('time');
ylabel('Amplitude');

% --- Executes on button press in view_all.
function view_all_Callback(hObject, eventdata, handles)
% hObject    handle to view_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.output_panel,'visible','off');        %Hide the output panel
guidata(hObject,handles);

% --- Executes on button press in load_sound.
function load_sound_Callback(hObject, eventdata, handles)
% hObject    handle to load_sound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 [filename,pathname]=uigetfile('*.wav','Select an image File');     %Select the file from the system
 [input_song,freq]=audioread(fullfile(pathname,filename));          %Get the data from the selected file
 audiowrite('input.wav', double(input_song), 48000);                %Store the data
 pro(hObject, eventdata, handles)                                   %Go to the processing function
 
 set(handles.start_panel,'visible','off');   %Hide the start panel
 
 
% --- Executes on button press in record.
function record_Callback(hObject, eventdata, handles)
% hObject    handle to record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.record_panel,'visible','on');       %Make the record panel visible
guidata(hObject,handles);

% --- Executes on button press in start_record.
function start_record_Callback(hObject, eventdata, handles)
% hObject    handle to start_record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = audiorecorder(48000, 8, 1);             %creates an 48000 Hz, 8-bit, 1 channel audiorecorder object
fprintf('Recording the audio\n');
record(a, 5);                               %record the audio for 5 seconds

f = waitbar(0, 'Speak Now...');             %Initialize the waitbar
for i=1:10
    waitbar(i/10, f,'Speak Now...');       
    pause(0.5)

end
close(f)                                    %Close the waitbar

b = getaudiodata(a);                        %get the recorded data
audiowrite('input.wav', b, 48000);           %store the recorded audio in the system
pro(hObject, eventdata, handles)            %Go to the processing function

set(handles.record_panel,'visible','off');  %Hide the record panel   
set(handles.start_panel,'visible','off');   %Hide the start panel



% --- Executes on button press in record_back.
function record_back_Callback(hObject, eventdata, handles)
% hObject    handle to record_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.record_panel,'visible','off');      %Hide the record panel
guidata(hObject,handles);

% --- Processing Function ---%
function pro(hObject, eventdata, handles)

%reading the speech file
[data, fs] = audioread('input.wav');            %get the audio data from the file named audio.wav
data = awgn(data, 10, 'measured');
data_Orig = data;                               %store the audio data in the variable
fs_Orig = fs;                                   %store the sampling frequency in the variable

handles.data = data_Orig;                       %store the data in the handles       
handles.fs = fs_Orig;                           %store the sampling rate in the handles  


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the energy_graph of the speech signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = data / abs(max(data));                   %normalize the entire data

f_d = 0.025;                                    %set frame duration as 25ms
data_frames = framing(data, fs, f_d);           %Split the signal into frames. It is like 0% overlap with rectangular window
frames = pre_process(data_frames);              %pre-process the frames obtained. 
                                                %This includes passing the frames through Hamming window 
                                                %and pre-emphasis of the signal
% calculate frames energy_graph
[r,c] = size(frames);                           %Get the number of rows and columns in the frames                                
ste = 0;                                        %initialize the energy to zero
for i = 1 : r
    ste(i) = sum(frames(i,:).^2);               %calculating energy of each frame
end                                            

ste = ste./max(ste);                            %normalize the energy data obtained

%Energy wave
f_size = round(f_d * fs);                       
ste_wave = 0;                                   %Initialize the energy wave to zero
for j = 1 : length(ste)
    l = length(ste_wave);
    ste_wave(l : l + f_size) = ste(j);          %Get the energy wave 
end
handles.ste_wave = ste_wave;                    %Store the energy wave

high_energy_signal_frames = data_frames;        %initialize the high energy signal to frames of original signal
high_to_low = 0;                                %to check if the wave changes from high to low
index = 0;                                      %index where the wave changes from high to low
for i = 1:length(ste)
     if ste(i) < 0.3                            %if the energy is less than 0.3
          if i >= 4 && index > 4                %if i is greater than 4 and index is greater than 4
              if high_to_low == 1               %if there if a change from high to low
                  if i == index+1 || i == index+2 || i == index+3 || i == index+4
                      continue;                 %keep the next four frames as it is
                  end
                  high_to_low = 0;              %again initialize the high to low indicator to zero
              end
          end
          high_energy_signal_frames([i,1],:) = 0;                                   %if the energy is less than 0.3, reduce it to zero
     else
         if i > 4
              high_energy_signal_frames([i - 1,1],:) = data_frames([i-1,1],:);      %restore the data of the previous four frames
              high_energy_signal_frames([i - 2,1],:) = data_frames([i-2,1],:);
              high_energy_signal_frames([i - 3,1],:) = data_frames([i-3,1],:);
              high_energy_signal_frames([i - 4,1],:) = data_frames([i-4,1],:);
             
              high_to_low = 1;                                                      %indicate that wave will change from high to low
              index = i;                                                            %index where the wave will change from high to low
         end
     end
end
        
high_energy_signal_frames_2 = high_energy_signal_frames;                            %backup the high energy signal

%High energy_graph Signal
high_energy_signal = (reshape(high_energy_signal_frames',1,[]))';                   %reconstruct signal i.e reshape it into original form

handles.energy = high_energy_signal;                                                %Store the data in handles


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the ZCR_GRAPH of the speech signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = data_Orig;                               %Get the original data
fs = fs_Orig;                                   %Get the original sampling frequency



data = data / abs(max(data));                   %normalize the entire data
f_d = 0.025;                                    %set frame duration as 25ms
frames = framing(data, fs, f_d);                %Split the signal into frames. It is like 0% overlap with rectangular window

% finding ZCR_GRAPH of all frames
[r,c] = size(frames);                           %Get the number of rows and columns in the frames 

for i = 1 : r
    x = frames(i, :);                           %get the indivisual frames
    ZCRf1(i) = 0;                               %Initialize the zero-crossing count of each frame to zero
    for k = 1:length(x) - 1
        if ((x(k) < 0) && (x(k + 1) > 0 ))      %If the sign changes from negative to positive
            ZCRf1(i) = ZCRf1(i) + 1;            %Increment the count by 1
   
        elseif ((x(k) > 0) && (x(k + 1) < 0))   %If the sign changes from positive to negative
            ZCRf1(i) = ZCRf1(i) + 1;            %Increment the count by 1
        end
    end
end

ZCRr1 = ZCRf1/length(x);                        %Calculate the Zero-Crossing Rate
ZCRr1 = ZCRr1/max(ZCRr1);                       %Normalise the ZCR

%ZCR wave
f_size = round(f_d * fs);                       %Calculate the frame size
zcr_wave = 0;                                   %Initialize the ZCR wave to zero
for j = 1 : length(ZCRr1)
    l = length(zcr_wave);
    zcr_wave(l : l + f_size) = ZCRr1(j);        %Get the ZCR wave
end
handles.zcr_wave = zcr_wave;                    %Store the ZCR wave


small_zcr_signal_frames = frames;               %initialize the small ZCR signal to frames of original signal
high_to_low = 0;                                %to check if the wave changes from high to low
index = 0;                                      %index where the wave changes from high to low
for i = 1:length(ZCRr1)
     if ZCRr1(i) >= 0.15                        %if the ZCR is more than 0.15
          if i >= 4 && index > 4                %if i is greater than 4 and index is greater than 4
              if high_to_low == 1               %if there if a change from high to low
                  if i == index+1 || i == index+2 || i == index+3 || i == index+4
                      continue;                 %keep the next four frames as it is
                  end
                  high_to_low = 0;              %again initialize the high to low indicator to zero
              end
          end
          small_zcr_signal_frames([i,1],:) = 0;                                 %if the ZCR is more than 0.3, reduce it to zero
     else
         if i > 4
              small_zcr_signal_frames([i - 1,1],:) = frames([i-1,1],:);         %restore the data of the previous four frames
              small_zcr_signal_frames([i - 2,1],:) = frames([i-2,1],:);
              small_zcr_signal_frames([i - 3,1],:) = frames([i-3,1],:);
              small_zcr_signal_frames([i - 4,1],:) = frames([i-4,1],:);
             
              high_to_low = 1;                                                  %indicate that wave will change from high to low
              index = i;                                                        %index where the wave will change from high to low
         end
     end
end
        
small_zcr_signal_frames_2 = small_zcr_signal_frames;                            %backup the small ZCR signal

%High energy_graph Signal
small_zcr_signal = (reshape(small_zcr_signal_frames',1,[]))';                   %reconstruct signal i.e reshape it into original form

handles.zcr = small_zcr_signal;                                                 %Store the data in the handles



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Comparing the voiced_graph signals with the original signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unvoiced = data_Orig;                                           %Initialize the unvoiced signal to the Original Signal
for i = 1:length(data_Orig)
   if high_energy_signal(i) ~= 0 && small_zcr_signal(i)~= 0     %If the signal does has both high Energy and small ZCR
       voiced(i) = data_Orig(i);                                %Add that part to voiced segment
       unvoiced(i) = 0;                                         %Make the unvoiced segment equal to zero
   else
       voiced(i) = 0;                                           %Make the voiced segment equal to zero
   end
end

handles.voiced = voiced;                                        %Store the data in handles
handles.unvoiced= unvoiced;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Passing the signal through filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 7;                                                          %Order to the filter
beginFreq = 1000/(fs/2);                                        
endFreq = (3000 / (fs/2));         
[b, a] = butter(n , [beginFreq, endFreq], 'bandpass');          %filter coefficients

fout = filter(b, a, voiced);                                    %Applying filter on the voiced signal

handles.op = fout;                                              %Store the data in handles


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Storing the out_graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

audiowrite('output.wav', fout, 48000);
audiowrite('unvoiced.wav', unvoiced, 48000);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 axes(handles.original_sound);                          %Graph of original signal on view all graph screen
 t = [0 : 1/fs : length(handles.data)/fs]; % time in sec
 t = t(1:end - 1);
 plot(t,handles.data); hold on;
 title('Original Sound');
 xlabel('time');
 ylabel('Amplitude');
 hold off;
 
 axes(handles.original_sound_graph);                    %Graph of original signal on output screen
 t = [0 : 1/fs : length(handles.data)/fs]; % time in sec
 t = t(1:end - 1);
 plot(t,handles.data); hold on;
 title('Original Sound');
 xlabel('time');
 ylabel('Amplitude');
 hold off;
 
 axes(handles.energy_graph);                            %Graph of high energy signal
 t = [0 : 1/fs : length(handles.energy)/fs]; % time in sec
 t = t(1:end - 1);
 plot(t, handles.energy,'g','LineWidth',2);
 title('High Energy Signal');
 xlabel('time');
 ylabel('Amplitude');
 
 axes(handles.zcr_graph);                               %Graph of small ZCR signal
 t = [0 : 1/fs : length(handles.zcr)/fs]; % time in sec
 t = t(1:end - 1);
 plot(t, handles.zcr,'g');
 title('Small ZCR Signal');
 xlabel('time');
 ylabel('Amplitude');
 
 axes(handles.unvoiced_graph);                          %Graph of unvoiced signal
 t = [0 : 1/fs : length(handles.unvoiced)/fs]; % time in sec
 t = t(1:end - 1);
 plot(t, handles.unvoiced,'r','LineWidth',2);
 title('Unvoiced');
 xlabel('time');
 ylabel('Amplitude');
 
 axes(handles.voiced_graph);                            %Graph of voiced signal
 t = [0 : 1/fs : length(handles.voiced)/fs]; % time in sec
 t = t(1:end - 1);
 plot(t, handles.voiced, 'r', 'LineWidth',2);
 title('Voiced Before filtering');
 xlabel('time');
 ylabel('Amplitude');
 
 axes(handles.out_graph);                               %Graph of output signal on view all screen
 t = [0 : 1/fs : length(handles.op)/fs]; % time in sec
 t = t(1:end - 1);
 plot(t, handles.op,'y','LineWidth',2);
 title('Output');
 xlabel('time');
 ylabel('Amplitude');
 
 axes(handles.output_graph);                            %Graph of output signal on output screen
 t = [0 : 1/fs : length(handles.op)/fs]; % time in sec
 t = t(1:end - 1);
 plot(t, handles.op,'y','LineWidth',2);
 title('Output');
 xlabel('time');
 ylabel('Amplitude');
 
 guidata(hObject,handles);


% --- Executes on button press in energy_wave.
function energy_wave_Callback(hObject, eventdata, handles)
% hObject    handle to energy_wave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% plot the STE with Signal
t = [0 : 1/handles.fs : length(handles.data)/handles.fs]; % time in sec
t = t(1:end - 1);
t1 = [0 : 1/handles.fs : length(handles.ste_wave)/handles.fs];
t1 = t1(1:end - 1);
figure(20)
plot(t,handles.data); hold on;
plot(t1,handles.ste_wave,'r','LineWidth',2);
xlabel('time');
ylabel('Amplitude');
legend('Speech Signal','Short Term Energy (Frame Energy)');
hold off;

% --- Executes on button press in zcr_wave.
function zcr_wave_Callback(hObject, eventdata, handles)
% hObject    handle to zcr_wave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% plot the ZCR with Signal
t = [0 : 1/handles.fs : length(handles.data)/handles.fs]; % time in sec
t = t(1:end - 1);
t1 = [0 : 1/handles.fs : length(handles.zcr_wave)/handles.fs];
t1 = t1(1:end - 1);
figure(13)
plot(t,handles.data); hold on;
plot(t1,handles.zcr_wave,'r','LineWidth',2);
xlabel('time');
ylabel('Amplitude');
legend('Speech Signal','Zero-Crossing Rate');
hold off;