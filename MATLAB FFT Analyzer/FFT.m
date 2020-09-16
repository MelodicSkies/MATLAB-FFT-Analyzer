function varargout = FFT(varargin)
% FFT MATLAB code for FFT.fig
%      FFT, by itself, creates a new FFT or raises the existing
%      singleton*.
%
%      H = FFT returns the handle to a new FFT or the handle to
%      the existing singleton*.
%
%      FFT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FFT.M with the given input arguments.
%
%      FFT('Property','Value',...) creates a new FFT or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FFT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FFT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FFT

% Last Modified by GUIDE v2.5 22-Apr-2019 06:49:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FFT_OpeningFcn, ...
                   'gui_OutputFcn',  @FFT_OutputFcn, ...
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
end
% End initialization code - DO NOT EDIT

% --- Executes just before FFT is made visible.
function FFT_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

set(handles.graphView,'Value',true);
set(handles.tableView,'Value',false);
set(handles.graphView2,'Value',false);

set(handles.graph,'Visible',true);
set(handles.table,'Visible',false);
set(handles.graph2,'Visible',false);

set(handles.filePath,'String',"");
set(handles.calibrationTest,'Value',false);
set(handles.differenceThreshold,'enable','off');
end

function varargout = FFT_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
end

function filePath_Callback(~, ~, ~)
%audio source filepath
end

function browse_Callback(~, ~, handles)
[filename, filepath] = uigetfile({'*.*';'*.mp3';'*.wav'},'Search sound file to be analyzed');
fullname = [filepath filename];
set(handles.filePath,'string',fullname);
axes(handles.graph);
cla reset;
axes(handles.graph2);
cla reset;
end

function graphView_Callback(~, ~, handles)
set(handles.graph,'visible','on');
set(handles.table,'visible','off');
set(handles.graph2,'visible','off')
set(handles.tableView,'Value',false);
set(handles.graphView2,'Value',false);

axes(handles.graph);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
end

function tableView_Callback(~, ~, handles)
set(handles.graph,'visible','off');
set(handles.table,'visible','on');
set(handles.graph2,'visible','off');
set(handles.graphView,'Value',false);
set(handles.graphView2,'Value',false);
end

function graphView2_Callback(~, ~, handles)
set(handles.graph2,'visible','on');
set(handles.graph,'visible','off');
set(handles.table,'visible','off');
set(handles.tableView,'Value',false);
set(handles.graphView,'Value',false);

axes(handles.graph2);
xlabel('Time (s)');
ylabel('Magnitude (dB)');
end

function ampFloor_Callback(~, ~, ~)
%minimum amplitude for a tone to be included
end

function calibrationTest_Callback(~, ~, handles)
if(get(handles.calibrationTest,'Value') == true)
    set(handles.ampFloor,'enable','off');
    set(handles.calibration,'enable','off');
    set(handles.customAmpFloor,'enable','off');
else
    set(handles.ampFloor,'enable','on');
    set(handles.calibration,'enable','on');
    set(handles.customAmpFloor,'enable','on');
end
end


function differenceThreshold_Callback(~, ~, ~)
%dB difference for two tones to be considered different after FFT
%not yet implemented
end

function customAmpFloor_Callback(~, ~, handles)
if(get(handles.customAmpFloor,'Value') == true)
    set(handles.ampFloor, 'enable', 'on');
else
    set(handles.ampFloor, 'enable', 'off');
    set(handles.ampFloor,'String',0.5);
end
end

function calibration_Callback(~, ~, handles)
%offset for calibrating 0 dB
end

function analyze_Callback(~, ~, handles)
filepath = get(handles.filePath,'String');
pattern = [".wav",".mp3"];

if endsWith(filepath,pattern,'IgnoreCase',true)
    [y,Fs] = audioread(filepath);
    channels = size(y,2);

    if channels == 2
        mono = sum(y, 2) / size(y, 2);
       
        L = (length(y));
        time_length = (L/Fs)%*1000; %ms
        
        NFFT = L;
        Y = fft(mono,NFFT)/Fs;
        
    else
        L = (length(y));
        time_length = (L/Fs)%*1000; %ms
        
        NFFT = L;
        Y = fft(y,NFFT)/Fs;
    end
    
     f = Fs/2*linspace(0,1,NFFT/2+1);
     calibrationDB = str2double(get(handles.calibration,'String'));
     ydb = mag2db(2*abs(Y(1:NFFT/2+1))) - calibrationDB;
     phase = angle(y(1:NFFT/2+1));
    
     axes(handles.graph);
     xlabel('Frequency (Hz)');
     ylabel('Magnitude (dB)')
     
     plot(f,ydb,'Parent', handles.graph);
     grid on;
     
     maxPeaks = {};
     
     if(get(handles.calibrationTest,'Value') == true)
         
         [M,I] = max(ydb);
         maxFreq = f(I); %index
         
         maxPeaks{end + 1, 1} = M;
         maxPeaks{end + 1, 2} = maxFreq;
         peakPhase = zeros(1,1);
         peakPhase(1,1) = phase(I);
         
         set(handles.calibration,'String',sprintf('%0.6f',M));
         set(handles.table,'Data',maxPeaks);
        
     else
         minAmp = str2double(get(handles.ampFloor,'String'));
         [peaks,freq] = findpeaks(ydb,'MinPeakHeight',minAmp,'MinPeakDistance',100); 
         %change 100 to custom separation distance if needed
         peakPhase = zeros(1,length(peaks));
         
         for i=1:length(freq)
             peakPhase(1,i) = phase((freq(i)));
             freq(i) = f(freq(i));
             maxPeaks{end + 1, 1} = peaks(i);
             maxPeaks{end, 2} = freq(i);
         end
         
         set(handles.table,'Data',maxPeaks);
         
     end
       
       axes(handles.graph2); 
       
       t = 0:1e-4:1e-2; %change the time parameters if graph looks off
       tempFreq = freq;
       tempAmp = peaks;
       tempPhase = peakPhase;
       
       y = zeros(length(tempFreq),length(t));
       for i=1:length(tempFreq)
           y(i,:) = tempAmp(i)*sin(2*pi*tempFreq(i)*t + (tempPhase(i)/180));
       end
       plot(t,y);
       grid on;
       %for i=1:length(maxPeaks{1})
           %y(i) = maxPeaks{i,1}*sin(2*pi*maxPeaks{i,2}*t + (peakPhase(i)/180));
       %end
       
       %plot(t,y);

           %plot(t,maxPeaks{i,1}*sin(2*pi*maxPeaks{i,2}*t + (peakPhase(i)/180)),'Parent',handles.graph2);
           %tempFreq = maxPeaks{i,2};
           %tempAmp = maxPeaks{i,1};
           %tempPhase = peakPhase(i);
           
           %f = tempFreq;
           %a = tempAmp;
           %y = a*sin(2*pi*f*t + (tempPhase/180));
           %plot(t,y,'Parent',handles.graph2);

           %t = 0:.01:time_length;
           %y = tempAmp*cos(2*pi*tempFreq*t*9.545 + (tempPhase/180));
           %y = tempAmp*sin(2*pi*t + (tempPhase/180));
           %y = tempAmp*sin(2*pi*tempFreq*t + (tempPhase/180));
           %plot(t,y);
       %end
       
       Amp = string(peaks) + 'dB';
      
       fid = fopen('AmplitudeFFT.txt','wt');
       for ii = 1:size(Amp)
       fprintf(fid,Amp(ii,:));
       fprintf(fid,'\n');
       end
       fclose(fid);
       
       Freq = freq;
       
       save('FreqlistFFT.txt', 'Freq','-ascii');
        
        %Test code IGNORE---------------------------------------------
        %ydb = mag2db(2*abs(Y(1:NFFT/2+1)));
        %plot(f,ydb);
        
        %ydb = 2*abs(Y(1:NFFT/2+1));
        %bx = NFFT*f/Fs + 1;
        %As = 20*log10(ydb(bx));
        %ydb(bx) = 0;
        %An = 10*log10(sum(ydb.^2));
        %SNR = As - An;
        
        %ydb2 = 20*log10(abs(Y(1:NFFT/2+1)));
        %plot(f,mag2db(abs(fft(Y))));
        %ydb = 20*log10((2*abs(Y(1:NFFT/2+1)))/NFFT);
        %ydb = 20*log10((2*abs(Y(1:NFFT/2+1)))/32768);
        %ydb = 20*log10(2*abs(Y(1:NFFT/2+1)));
        %plot(f,ydb);
        
        %[PKAmp,PKTime]=findpeaks(f,ydb);
        %[pks,frqs] = findpeaks(ydb,f);
                
        %bin_vals = [0 : L-1];
        %L_2 = ceil(L/2);
        %plot(f, 10*log10(Y(1:L_2)));
        
        
        %plot(abs(fft(y,L)));
        %f = Fs/2*linspace(0,1,NFFT/2+1);
        %plot(f,abs(Y(1:NFFT/2+1)));
        
        %f = Fs/2*linspace(0,1,NFFT/2+1);
        %axes(handles.ToneNoisePlot);
        %plot(f,abs(Y(1:NFFT/2+1)));
        
        %[B,IX] = sort(2*abs(Y(1:NFFT/2+1)));
        %BFloor = 0.1; %BFloor is the minimum amplitude value (ignore small values)
        %Amplitudes = B(B >= BFloor); %find all amplitudes above the BFloor
        %Frequencies = f(IX(1+end-numel(Amplitudes):end)); %frequency of the peaks
        
        %Amp = reshape(Amplitudes,[],1);
        %Freq = reshape(Frequencies,[],1);
        
        %save('FreqlistTest.txt', 'Freq','-ascii');
        %save('AmplitudeTest.txt', 'Amp','-ascii');
        
        %end
        %plot(Amp, Freq);
        %title('Single-Sided Amplitude Spectrum of S(t)')


else
    warndlg('Only .mp3 and .wav files supported.');
end
set(handles.graph2,'visible','off');
set(handles.graph,'visible','on');
set(handles.table,'visible','off');
set(handles.tableView,'Value',false);
set(handles.graphView,'Value',true);
set(handles.graphView2,'Value',false);
end
%coded with 16 bit audio files in mind, mag2dB equation assumes this,
%other bits per sample will have to be converted differently according to
%their dynamic range
