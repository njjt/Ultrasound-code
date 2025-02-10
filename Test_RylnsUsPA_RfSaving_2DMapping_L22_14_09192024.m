% log: 2016/06/01 modifide the time interval between US Foci shift from 372
% microsecs to 100  microsecs, which could save 29 ms in total for 128 foci
% US imaging.

%0919 test changed line 1142 1098

clear all
close all

%% intialize Zaber Stage
delete(instrfind);
comX = serial('com3');
comY = serial('com4');

fopen(comX); fopen(comY);
zaberHome(comX); zaberHome(comY);
% zaberReset(comX);
% zaberReset(comY);
% zaberSetSpeed(comX, 12);
% zaberSetSpeed(comY, 12);
% 
% % zaberSetZeroHoldCurrent(comX); 
% zaberSetHoldCurrent(25, comX)
% 
% xHoldCurrent = zaberReadHoldCurrent(comX);
% disp(['xHoldCurrent is :' , num2str(xHoldCurrent)]);
% 
% 
% % zaberSetZeroHoldCurrent(comY); 
% zaberSetHoldCurrent(25, comY)
% 
% yHoldCurrent = zaberReadHoldCurrent(comY);
% disp(['yHoldCurrent is :' , num2str(yHoldCurrent )]);

xPos = 0;
yPos = 0;

%% linear scan parameters 
scanFlag = 0;
savCurrentFrameFlag = 0;

stageCurPosX = 0;     % mm
scanEdPosX = 0;      % mm
scanStepX = 0.2;      % mm
scanAxisX = stageCurPosX: scanStepX: scanEdPosX;

stageCurPosY = 0;     % mm
scanEdPosY = 0;      % mm
scanStepY = 9.6;      % mm
scanAxisY = stageCurPosY: scanStepY: scanEdPosY;

ScanUsImg = [];
ScanPaImg = [];
frameNo = 0;
tempFrame = [];
tempFrameUs = [];
tempFramePa = [];
dest_offset = 0;
UsFileName = [];
PaFileName = [];
currentFrameImagP = [];
FilePathForCurFrame = [];
flagSav = 0;
NumFrames2Sav = 0;
%%
paramFname = [];
UsFname = [];
PaFname = [];
UsIQFname = [];
PaIQFname = [];


%% IQ data saving

tempFrameIQ = [];
tempFrameUsIQdata = [];
tempFramePaIQdata = [];


%% Set PA parameters
oneway = 1;     % (logical) oneway=1 turns off the transmitters by setting TX.Apod to zero for all transmitters
flash2Qdelay = 0; % microseconds between trigger input and start of acquisition (which outputs a trigger pulse at time=0)
ne = 1;         % ne = number of acquisitions in PA ensemble for coherent addition in I/Q buffer.
PA_Angle = 0;   % angle of transmit plane wave for use in testing the PA mode in simulation
PA_PRF = 10;   % PA PRF in Hz. To be set in accordance with laser rep rate, when not using the input trigger mode. 
                % When using the input trigger mode, remove the TTNA (SeqControl(4)) from the PA events to avoid
                % "missed TTNA" messages.

    if oneway==1
        disp(' *** PhotoAcoustic mode: Using one-way receive-only reconstruction ***')
    else
        disp(' *** Ultrasound Transmit mode: Using conventional T/R reconstruction ***')
    end
% ===============================================================================


%% P parameter started to use from V 3.0
P.startDepth = 5; %5
% changed to 5 for testing new US probe - Lu 20161027, roginal is 150
P.endDepth = 300;     %300       % Acquisition depth in wavelengths
P.txFocus = 80;  % Initial transmit focus.
P.numRays = 128;              % no. of Rays

nr = P.numRays;     % num of rays
%%
% Define system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L22-14v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
% Trans.frequency = 6.25;  % frequency in megaHertz
% note nominal center frequency in computeTrans is 6.25 MHz
Trans = computeTrans(Trans);  % L11-4v transducer is 'known' transducer so we can use computeTrans.

RcvProfile.LnaZinSel = 31;  %% new set up for L22-4v

nElem = Trans.numelements;
Trans.maxHighVoltage = 30;  % set maximum high voltage limit for pulser supply.

%% Specify SFormat structure array.
% US



%% Specify PData structure array.
% - US
% PData(1).sFormat = 1;      % use first SFormat structure.
PData(1).pdeltaX = 1;
PData(1).pdeltaZ = 0.5;
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).pdeltaZ);
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).pdeltaX);
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr

% - specify 128 Region structures.
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','Rectangle',...
                    'Position',[0,0,P.startDepth],...
                    'width',Trans.spacing,...
                    'height',P.endDepth-P.startDepth)),1,128);
% - set position of regions to correspond to beam spacing.
for i = 1:128
    PData(1).Region(i).Shape.Position(1) = (-63.5 + (i-1))*Trans.spacing;
end
PData(1).Region = computeRegions(PData(1));


% - PA
% PData(2).sFormat = 2;      % use first SFormat structure.
PData(2).pdeltaX = 1;
PData(2).pdeltaZ = 0.5;
PData(2).Size(1) = ceil((P.endDepth-P.startDepth)/PData(2).pdeltaZ); % startDepth, endDepth and pdelta set PData.Size.
PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).pdeltaX);
PData(2).Size(3) = 1;      % single image page
PData(2).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.


%% Specify Resources.
Resource.RcvBuffer.datatype = 'int16';
Resource.RcvBuffer.rowsPerFrame = 983040+2048; % max depth of 384 * 2 for round trip * 4 smpls/wave * 2*nr acquisitions 
Resource.RcvBuffer.colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer.numFrames = 10;

% 1 IQ buffer * 1 Frame for US
Resource.InterBuffer(1).datatype = 'complex';
% Resource.InterBuffer(1).rowsPerFrame = 1024; % this is for greatest depth
% Resource.InterBuffer(1).colsPerFrame = PData(1).Size(2);
Resource.InterBuffer(1).numFrames = 10;  % one intermediate buffer needed.

% 1 IQ buffer * 1 Frame for PA
Resource.InterBuffer(2).datatype = 'complex';
% Resource.InterBuffer(2).rowsPerFrame = 1024; % this is for greatest depth
% Resource.InterBuffer(2).colsPerFrame = PData(2).Size(2);
Resource.InterBuffer(2).numFrames = 10;  % one intermediate buffer needed.

% 1 Image buffer * 10 Frames for US
Resource.ImageBuffer(1).datatype = 'double';
% Resource.ImageBuffer(1).rowsPerFrame = 1024;
% Resource.ImageBuffer(1).colsPerFrame = PData(1).Size(2);
Resource.ImageBuffer(1).numFrames = 10;

% 1 Image buffer * 10 Frames for PA
Resource.ImageBuffer(2).datatype = 'double';
% Resource.ImageBuffer(2).rowsPerFrame = 1024;
% Resource.ImageBuffer(2).colsPerFrame = PData(2).Size(2);
Resource.ImageBuffer(2).numFrames = 10;

% % 1 Image buffer * 10 Frames for PA + US
% Resource.ImageBuffer(3).datatype = 'double';
% Resource.ImageBuffer(3).rowsPerFrame = 1024;
% Resource.ImageBuffer(3).colsPerFrame = PData(2).Size(2);
% Resource.ImageBuffer(3).numFrames = 10;

%% Display Window
Resource.DisplayWindow(1).Title = 'RyLns-US';
Resource.DisplayWindow(1).pdelta = 1;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).pdeltaX/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).pdeltaZ/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);


Resource.DisplayWindow(2).Title = 'PA';
Resource.DisplayWindow(2).pdelta = 1;
ScrnSize = get(0,'ScreenSize');
DwWidthP = ceil(PData(2).Size(2)*PData(2).pdeltaX/Resource.DisplayWindow(2).pdelta);
DwHeightP = ceil(PData(2).Size(1)*PData(2).pdeltaZ/Resource.DisplayWindow(2).pdelta);
Resource.DisplayWindow(2).Position = [250,(ScrnSize(4)-(DwHeightP+150))/2, ...  % lower left corner position
                                      DwWidthP, DwHeightP];
Resource.DisplayWindow(2).ReferencePt = [PData(2).Origin(1),PData(2).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).Colormap = gray(256);


% Specify TW structure array.
% - TW for US
TW(1).type = 'parametric';
% TW(1).Parameters = [Trans.frequency,.67,2,1];   % A, B, C, D
TW(1).Parameters = [18,.67,2,1];   % A, B, C, D

% - TW for PA, if relevant
TW(2).type = 'parametric';
% TW(2).Parameters = [Trans.frequency,.67,2,1];   % A, B, C, D
TW(2).Parameters = [18,.67,2,1];   % A, B, C, D

% Specify nr TX structure arrays. Transmit centered on element n in the array for event n.
% txFocus = 320;  % Initial transmit focus.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, P.numRays+1);

% Determine TX aperture based on focal point and desired f number.
txFNum = 20;  % set to desired f-number value for transmit (range: 1.0 - 20)
P.numTx = round((P.txFocus/txFNum)/Trans.spacing); % no. of elements in 1/2 aperture.
txNumEl = floor(P.numTx/2);

if txNumEl > (Trans.numelements/2 - 1), txNumEl = floor(Trans.numelements/2 - 1); end   
% txNumEl is the number of elements to include on each side of the
% center element, for the specified focus and sensitivity cutoff.
% Thus the full transmit aperture will be 2*txNumEl + 1 elements.
%display('Number of elements in transmit aperture:');
%disp(2*txNumEl+1);
           
% - Set event specific TX attributes.
for n = 1:nr   % 128 transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = [(-63.5 + (n-1))*Trans.spacing, 0.0, 0.0];
    % Set transmit Apodization so (1 + 2*TXnumel) transmitters are active.
    lft = n - txNumEl;
    if lft < 1, lft = 1; end;
    rt = n + txNumEl;
    if rt > Trans.numelements, rt = Trans.numelements; end;
    TX(n).Apod(lft:rt) = 1.0;
%     TX(n).Apod(lft:rt) = 0.0;
    TX(n).Delay = computeTXDelays(TX(n));
end


% PA mode, TX should be shut down
TX(nr+1).waveform = 1;            % use 1st TW structure.
TX(nr+1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(nr+1).focus = 0;
TX(nr+1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(nr+1).Apod = zeros(1,Trans.numelements);
TX(nr+1).Delay = computeTXDelays(TX(nr+1));



%% Specify TPC structures.
TPC(1).name = 'US';
TPC(1).maxHighVoltage = 30;

% This allows one to use different transmit profile for PA ... only relevant if transmitters are active
% --- currently TPC(2) is not used ---
TPC(2).name = 'PA';     
TPC(2).maxHighVoltage = 30;

%% Specify TGC Waveform structure.

% TGC 1 - US
TGC(1).CntrlPts = [300,511,716,920,1023,1023,1023,1023]; %[0,138,260,287,385,593,674,810];
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));

% TGC 2 - PA, 1st zero is to suppress the reflection
TGC(2).CntrlPts = [0 662 662 662 662 662 662 662]; % TO BE MODIFIED HERE AFTER COLLECTING PA DATA ;
TGC(2).rangeMax = P.endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));


%% Set LNA gain to be 24 for PA mode ---Lu
% default LNA gain for US
RcvProfile(1).LnaGain = 18;     % 12, 18, or 24 dB  (18=default)
RcvProfile(1).condition = 'immediate';

% Larger LNA gain for PA
RcvProfile(2).LnaGain = 24;     % 12, 18, or 24 dB  (18=default)
RcvProfile(2).condition = 'immediate';


%% Specify Receive structure arrays. 
% - We need 2*nr+1 Receives for every frame (synthetic aperture).


% sampling center frequency is 15.625, but we want the bandpass filter
% centered on the actual transducer center frequency of 18 MHz with 67%
% bandwidth, or 12 to 24 MHz.  Coefficients below were set using
% "G3_BPFdevelopment" with normalized cf=1.15 (18 MHz), bw=0.85, 
% xsn wdth=0.41 resulting in -3 dB 0.71 to 1.6 (11.1 to 25 MHz), and
% -20 dB 0.57 to 1.74 (8.9 to 27.2 MHz)
% 
BPF1 = [ -0.00009 -0.00128 +0.00104 +0.00085 +0.00159 +0.00244 -0.00955 ...
         +0.00079 -0.00476 +0.01108 +0.02103 -0.01892 +0.00281 -0.05206 ...
         +0.01358 +0.06165 +0.00735 +0.09698 -0.27612 -0.10144 +0.48608 ];


maxAcqLengthUS = sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2) - P.startDepth;
maxAcqLengthPA = sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2) - P.startDepth;


wlsPer128 = 128/(2*2); % wavelengths in 128 samples for 2 samplesPerWave
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLengthUS/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'samplesPerWave', 2,...
                        'InputFilter', BPF1, ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (nr+1)*Resource.RcvBuffer(1).numFrames); % +1 for PA at last Acq
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (nr+1)*(i-1);
    % Acq for US
    for j = 1:P.numRays
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
    
    % PA acquisitions
    
    Receive(k+j+1).framenum = i;
    Receive(k+j+1).acqNum = j+1;        % PA acqNums continue after 2D
    Receive(k+j+1).startDepth = P.startDepth;
    Receive(k+j+1).samplesPerWave = 2;
    Receive(k+j+1).endDepth = P.startDepth + wlsPer128*ceil(maxAcqLengthPA/wlsPer128);
    Receive(k+j+1).TGC = 2;           % TGC(1) is tied to the GUI sliders
    
end

%% Specify Recon structure array. 
Recon(1) = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', zeros(1,nr));
           
% - Set RINums values for frame.
Recon(1).RINums(:) = (1:nr)';  % 2*nr ReconInfos needed for nr rays, 2-1 synthetic aperture

Recon(2) = struct('senscutoff', 0.6, ...
               'pdatanum', 2, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [2,1], ...
               'ImgBufDest', [2,-1], ...
               'RINums', nr+1);


% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 0, ...  % replace IQ data first half SA.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, nr+1);
               
               
% - Set specific ReconInfo attributes -- US.
for j = 1:nr 
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end

               
ReconInfo(nr+1).txnum = nr+1;
ReconInfo(nr+1).rcvnum = nr+1;
ReconInfo(nr+1).regionnum = 0;
               

%% Specify Process structure arrays.
cpt = 22;       % define here so we can use in UIControl below
cpers = 10;     % define here so we can use in UIControl below
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...          % number of ImageBuffer to process.
                         'framenum',-1,...          % frame number in source buffer (-1 => lastFrame)
                         'pdatanum',1,... 
                         'norm',1,...               % normalization method(1=fixed)
                         'pgain',20.0,...            % pgain is image processing gain
                         'persistMethod','simple',...
                         'persistLevel',5,...       % persistence averaging of image ( 0 < pers < 100)
                         'interp',1,...             % method of interpolation (1 = 4pt interp)
                         'compression',0.3,...      % X^0.5 is Amplitude = sqrt(Intensity)
                         'reject',10,...
                         'mappingMode','full',...
                         'displayWindow',1,...
                         'display',1};              % don't display image after processing

Process(2).classname = 'Image';                     % image display for color data.
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',2,...          % number of buffer to process.
                         'framenum',-1,...          % frame number in source buffer (-1 => lastFrame)
                         'pdatanum',2,... 
                         'norm',1,...               % normalization method(2=none)
                         'pgain',60.0,...            % pgain is image processing gain
                         'persistMethod','simple',...
                         'persistLevel',cpers,...
                         'interp',1,...             % method of interpolation (1=4pt interp)
                         'mappingMode','full',...
                         'reject',10,...
                         'threshold', 22,...
                         'displayWindow',2,...
                         'display',1};              % do display image after processing

% External Processing to control Zaber stage for scanning
Process(3).classname = 'External';
Process(3).method = 'stageMov';
Process(3).Parameters = {'srcbuffer','none'};            

% External Processing - US Display and svae Frame - Lu.
Process(4).classname = 'External';
Process(4).method = 'UsDisplaySave';
Process(4).Parameters = {'srcbuffer','image',... % name of buffer to process.
    'srcbufnum',1,...
    'srcframenum',-1,... % process the most recent frame.
    'dstbuffer','none'};

% External Processing - US Display and svae Frame - Lu.
Process(5).classname = 'External';
Process(5).method = 'PaDisplaySave';
Process(5).Parameters = {'srcbuffer','image',... % name of buffer to process.
    'srcbufnum',2,...
    'srcframenum',-1,... % process the most recent frame.
    'dstbuffer','none'};

% External Processing - RfData of US sav - Lu.
Process(6).classname = 'External';
Process(6).method = 'UsIQSave';
Process(6).Parameters = {'srcbuffer','inter',... % name of buffer to process.
    'srcbufnum',1,...
    'srcframenum',1,... % process the most recent frame.
    'dstbuffer','none'};

% External Processing - RfData of US sav - Lu.
Process(7).classname = 'External';
Process(7).method = 'PaIQSave';
Process(7).Parameters = {'srcbuffer','inter',... % name of buffer to process.
    'srcbufnum',2,...
    'srcframenum',1,... % process the most recent frame.
    'dstbuffer','none'};


% Specify SeqControl structure arrays.
%  - Time between acquisitions in usec
% t1 = 4*384*.2 + 20; % 44 msec at max Receive.endDepth.
t1 = 500; % 44 msec at max Receive.endDepth.
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = t1;
%  - Time between frames at 12 fps at max Receive.endDepth.
SeqControl(2).command = 'timeToNextAcq';
% SeqControl(2).argument = round((1e+06-256*t1*12)/12);
SeqControl(2).argument = 1000;

% -- Change to Profile 2 (PA)
SeqControl(3).command = 'setTPCProfile';
SeqControl(3).condition = 'next';
SeqControl(3).argument = 2;
% -- Time between 2D acquisition and PA ensemble. Set to allow time for profile change.
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = 100; % time in usec
% -- PRF for PA ensemble
SeqControl(5).command = 'timeToNextAcq';
SeqControl(5).argument = round(1/(PA_PRF*1e-06)); % (10 msecs for PA_PRF=100 Hz) 
% -- Change to Profile 1 (2D)
SeqControl(6).command = 'setTPCProfile';
SeqControl(6).condition = 'next';
SeqControl(6).argument = 1;
% -- Time between PA and next 2D acquisition. Set to allow time for profile change.
SeqControl(7).command = 'timeToNextAcq';
SeqControl(7).argument = 700; % time in usec
% -- Jump back to start.
SeqControl(8).command = 'jump';
SeqControl(8).argument = 1;
% set receive profile
SeqControl(9).command = 'setRcvProfile';        
SeqControl(9).argument = 1;  % profile 1 has default gain
SeqControl(10).command = 'setRcvProfile';        
SeqControl(10).argument = 2;  % profile 2 has +6 dB gain
% output trigger
SeqControl(11).command = 'triggerOut';
% input trigger 
SeqControl(12).command = 'pause';
SeqControl(12).condition = 'extTrigger';
SeqControl(12).argument = 17; % Trigger input 1, enable with rising edge 
% noop delay between trigger in and start of acquisition
SeqControl(13).command = 'noop';
SeqControl(13).argument = fix(flash2Qdelay)*5; % noop counts are in 0.2 microsec increments
% SeqControl(13).argument = 3*5; % noop counts are in 0.2 microsec increments
% SeqControl(13).argument = 1e3; % noop counts are in 0.2 microsec increments
% sync command
SeqControl(14).command = 'sync';

% noop delay between IQ buffer reading 
SeqControl(15).command = 'noop';
SeqControl(15).argument = 0*5; % noop counts are in 0.2 microsec increments
% SeqControl(13).argument = 1e3; % noop counts are in 0.2 microsec increments
% sync command
SeqControl(16).command = 'sync';

nsc = 17;  % next SeqControl number

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 2;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:nr                      % Acquire frame
        Event(n).info = 'Aqcuisition.';
        Event(n).tx = j;   % use next TX structure.
        Event(n).rcv = (nr+1)*(i-1)+j;   
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 1; % seqCntrl
        n = n+1;
    end
    % Replace last events SeqControl for inter-frame timeToNextAcq.
   Event(n-1).seqControl = [2,10];
   
   
   % Wait for input trigger from flash lamp firing
    Event(n).info = 'Wait for Trigger IN';
    Event(n).tx = 0;        % no TX
    Event(n).rcv = 0;       % no Rcv
    Event(n).recon = 0;     % no Recon
    Event(n).process = 0;   % no Process
    Event(n).seqControl = 12; % Trigger in
    n = n+1;

    % Pause for optical buildup
    Event(n).info = 'noop and sync';
    Event(n).tx = 0;        % no TX
    Event(n).rcv = 0;       % no Rcv
    Event(n).recon = 0;     % no Recon
    Event(n).process = 0;   % no Process
    Event(n).seqControl = [13,14]; % pause and sync
    n = n+1;

    % send trigger output at start of every PA acquisition to fire Q-switch
    Event(n).info = 'Acquire PA event';
    Event(n).tx = nr+1;      % use next TX structure after 2D (need tx event even if TX.Apod all zero)
    Event(n).rcv = (nr+1)*(i-1)+j+1;   % rcv PA data
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = 10; % no TTNA to prevent "missed time" message, and trigger out every time
    %         Event(n).seqControl = [4,10]; % TTNA for PA when not using trigger input, and trigger out every time
    n = n+1;
   
   
    Event(n-1).seqControl = [11,7,9]; % replace last PA acquisition Event's seqControl with longer TTNA and RCV profile change
   

    TransferToHost_Reference = nsc;    
    Event(n).info = 'Transfer Data';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % display processing
    Event(n).seqControl = nsc; % transferToHost
    n = n+1;
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer (each needs a different value of nsc)
      nsc = nsc+1;   
    
%     % Pause for IQ buffer reading
%     Event(n).info = 'noop';
%     Event(n).tx = 0;        % no TX
%     Event(n).rcv = 0;       % no Rcv
%     Event(n).recon = 0;     % no Recon
%     Event(n).process = 0;   % no Process
%     Event(n).seqControl = 15; % pause and sync
%     n = n+1;     

    Event(n).info = 'recons and US process'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;  % reconstruction for US
    Event(n).process = 1;    % process US
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1;

    Event(n).info = 'recons and PA process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 2;      % no reconstruction
    Event(n).process = 2;    % processing for PA, and display overlay of PA and 2D
    Event(n).seqControl = 0; % no seqCntrl
    n=n+1; 
    
    
    
    %  Display & Sav   
    Event(n).info = 'US display and sav'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;  % reconstruction for both 2D and PA
    Event(n).process = 4;    % process 2D
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1;

    Event(n).info = 'PA image display and sav';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 5;    % processing for PA, and display overlay of PA and 2D
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1; 
    
    Event(n).info = 'US IQ sav';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 6;    % processing for PA, and display overlay of PA and 2D
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1;
%     
    Event(n).info = 'PA IQ sav';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 7;    % processing for PA, and display overlay of PA and 2D
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1;
    
%     % sync after IQ buffer reading
%     Event(n).info = 'sync';
%     Event(n).tx = 0;        % no TX
%     Event(n).rcv = 0;       % no Rcv
%     Event(n).recon = 0;     % no Recon
%     Event(n).process = 0;   % no Process
%     Event(n).seqControl = 16; % pause and sync
%     n = n+1; 

    Event(n).info = 'Stage Move'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;  % reconstruction for US
    Event(n).process = 3;    % process US
    Event(n).seqControl = 0; % no seqCntrl      

    if floor(i/frameRateFactor) == i/frameRateFactor     % Exit to Matlab only every 3rd frame to prevent slowdown
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;   % no Processing
Event(n).seqControl = 8;

%% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Current Stage X position Slider - Lu
UI(2).Control = {'UserB6','Style','VsSlider','Label','X Pos.', 'FontUnits',4,'SliderMinMaxVal',[0, 150, 0],...
                 'SliderStep',[0.002, 0.02],'ValueFormat','%6.3f'};
UI(2).Callback = text2cell('%-UI#2Callback');

% - Scan X End pos. Slider - Lu
UI(3).Control = {'UserB5','Style','VsSlider','Label','X End Pos.','SliderMinMaxVal',[0, 150, 0],...
                 'SliderStep',[0.002, 0.02],'ValueFormat','%6.3f'};
UI(3).Callback = text2cell('%-UI#3Callback');

% - Xaxis Scan Step size Slider - Lu
UI(4).Control = {'UserB4','Style','VsSlider','Label','X step','SliderMinMaxVal',[0, 5, 0.2],...
                 'SliderStep',[0.01, 0.1],'ValueFormat','%6.3f'};
UI(4).Callback = text2cell('%-UI#4Callback');


% - Current Stage Y position Slider - Lu
UI(5).Control = {'UserB3','Style','VsSlider','Label','Y Pos.', 'FontUnits',4,'SliderMinMaxVal',[0, 75, 0],...
                 'SliderStep',[0.002, 0.02],'ValueFormat','%6.3f'};
UI(5).Callback = text2cell('%-UI#5Callback');

% - Scan Y End pos. Slider - Lu
UI(6).Control = {'UserB2','Style','VsSlider','Label','Y End Pos.','SliderMinMaxVal',[0, 75, 0],...
                 'SliderStep',[0.002, 0.02],'ValueFormat','%6.3f'};
UI(6).Callback = text2cell('%-UI#6Callback');

% - Yaxis Scan Step size Slider - Lu
UI(7).Control = {'UserB1','Style','VsSlider','Label','Y step','SliderMinMaxVal',[0, 20, 12.8],...
                 'SliderStep',[0.01, 0.1],'ValueFormat','%6.3f'};
UI(7).Callback = text2cell('%-UI#7Callback');



% - Scan start Button - Lu
UI(8).Control = {'UserC1','Style','VsPushButton','Label','StartScan?'};
UI(8).Callback = text2cell('%-UI#8Callback');


% - UI save push button for current frame - Lu
UI(9).Control = {'UserC3','Style','VsPushButton','Label','Save CurrentFrame'};
UI(9).Callback = text2cell('%-UI#9Callback');

% - Transmit focus change
UI(10).Control = {'UserA1','Style','VsSlider','Label','TX Focus',...
                 'SliderMinMaxVal',[20,320,100],'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(10).Callback = text2cell('%-UI#10Callback');


%% PA Recon Parameter Cal and Set
NoElem = 64;
[ImgDepthAxis, xAxis, zAxis, ApertureAxis] = myPA_ReconParameterCal_v307(NoElem);
% calculate US Display Axis
[Us_zAxis, Us_xAxis, Us_zAxisNum, Us_xAxisNum] = myUS_DisplayParameterCal_v307();

%% Call External Control or display function
EF(1).Function = text2cell('%EF#1');        % stage scanning control
EF(2).Function = text2cell('%EF#2');        % Us display & save
EF(3).Function = text2cell('%EF#3');        % PA display & save

EF(4).Function = text2cell('%EF#4');        % US IQ save
EF(5).Function = text2cell('%EF#5');        % PA IQ save

%% Save all the structures to a .mat file, and run VSX automatically
filename = ('L22-14vFlashAngleUS_PA');   % define variable 'filename' to permit VSX to skip user query for matfile
save (filename)                 % save the structures to a matfile
VSX                             % invoke VSX automatically when running this Setup script

return


%% **** Callback routines to be encoded by text2cell function. ****
% ---------------------------------------------------------------------
%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%-UI#2Callback -  X Scan start pos. update
    flag = evalin('base', 'scanFlag');
    if ~flag   % only work when scanning not going on
        newPosX = UIValue;
        pX = evalin('base', 'comX');
%         fclose(pX);
%         fopen(pX);
%         keyboard
        stageCurPosX = evalin('base', 'stageCurPosX');
        zaberMove_rel_wthrespnd(newPosX - stageCurPosX, 0 , pX);
        pause(0.1)
        m=zaberPos(pX)
        assignin('base','stageCurPosX', m);
        evalin('base', 'scanAxisX = stageCurPosX: scanStepX: scanEdPosX;'); %update scanAxis
    end
%-UI#2Callbacko

%-UI#3Callback - update X scan ed position
    flag = evalin('base', 'scanFlag');
    if ~flag   % only work when scanning not going on
        assignin('base','scanEdPosX', UIValue);
        evalin('base', 'scanAxisX = stageCurPosX: scanStepX: scanEdPosX;'); %update scanAxis
    end
%-UI#3Callback

%-UI#4Callback - update X scan step size
    flag = evalin('base', 'scanFlag');
    if ~flag   % only work when scanning not going on       
        assignin('base','scanStepX', UIValue);
        evalin('base', 'scanAxisX = stageCurPosX: scanStepX: scanEdPosX;'); %update scanAxis 
    end
%-UI#4Callback

%-UI#5Callback -  Y Scan start pos. update
    flag = evalin('base', 'scanFlag');
    if ~flag   % only work when scanning not going on
        newPosY = UIValue;
        pY = evalin('base', 'comY');
%         fclose(pX);
%         fopen(pX);
        stageCurPosY = evalin('base', 'stageCurPosY');
        zaberMove_rel_wthrespnd(newPosY - stageCurPosY, 0 , pY);
        pause(0.1)
        m=zaberPos(pY)
        assignin('base','stageCurPosY', m);
        evalin('base', 'scanAxisY = stageCurPosY: scanStepY: scanEdPosY;'); %update scanAxis
    end
%-UI#5Callback

%-UI#6Callback - update Y scan ed position
    flag = evalin('base', 'scanFlag');
    if ~flag   % only work when scanning not going on
        assignin('base','scanEdPosY', UIValue);
        evalin('base', 'scanAxisY = stageCurPosY: scanStepY: scanEdPosY;'); %update scanAxis
    end
%-UI#6Callback

%-UI#7Callback update Y scan step size
    flag = evalin('base', 'scanFlag');
    if ~flag   % only work when scanning not going on       
        assignin('base','scanStepY', UIValue);
        evalin('base', 'scanAxisY = stageCurPosY: scanStepY: scanEdPosY;'); %update scanAxis 
    end
%-UI#7Callback

%-UI#8Callback - set scanFlag 1 to start scanning
    flag = evalin('base', 'scanFlag');
    if ~flag
        comX = evalin('base', 'comX');
        comY = evalin('base', 'comY');
        
        fclose(comX); fopen(comX);
        fclose(comY); fopen(comY);
        
        cPosX = zaberPos(comX);
        cPosY = zaberPos(comY);
        
        disp([cPosX, cPosY]);
        
        scanAxisX = evalin('base', 'scanAxisX');
        scanAxisY = evalin('base', 'scanAxisY');
        
        mappingMatrix = createMappingMatrix(scanAxisX, scanAxisY);
        
        assignin('base', 'mappingMatrix', mappingMatrix);
        
        r = evalin('base', 'Us_zAxisNum');
        c = evalin('base', 'Us_xAxisNum');
        
        M = r; 
        N = c * length(mappingMatrix(:,1));
        
        % initiate temp frames for scanning images 
       evalin('base', ['tempFrameUs=zeros(', num2str(r),',',num2str(N),',''single'');']);
       evalin('base', ['tempFramePa=zeros(', num2str(r),',',num2str(N),',''single'');']);
       evalin('base', ['tempFrame=zeros(', num2str(r),',',num2str(c),',''single'');']);
       % initiate temp IQ frames
       evalin('base', ['tempFrameIQ=complex(zeros(', num2str(r),',',num2str(c),'));']);
       evalin('base', ['tempFrameUsIQdata=complex(zeros(', num2str(r),',',num2str(N),'));']);
       evalin('base', ['tempFramePaIQdata=complex(zeros(', num2str(r),',',num2str(N),'));']);
       
       % intitial US+PA data matrix to store scan data
       assignin('base','frameNo', 1);
       % flag to save frames
       assignin('base', 'flagSav',1);
       assignin('base', 'NumFrames2Sav', length(mappingMatrix(:,1)));
       
       disp([ scanAxisX(1), scanAxisX(1) - cPosX])
       disp([ scanAxisY(1), scanAxisY(1) - cPosY])
       
       zaberMove_rel_wthrespnd((scanAxisX(1) - cPosX),0 , comX); 
       pause(0.2);
       zaberMove_rel_wthrespnd((scanAxisY(1) - cPosY),0 , comY);
       pause(0.2);
       
       assignin('base','scanFlag', 1);
       assignin('base','stageCurPosX', scanAxisX(1));
       assignin('base','stageCurPosY', scanAxisY(1));

    end
    
%-UI#8Callback

%-UI#9Callback - save current frame
    flagSav = evalin('base', 'flagSav'); 
    if ~flagSav
        assignin('base','flagSav', 1);
        % UI to get # of frames to save
%         keyboard
        prompt = {'# of Frames to save:'};
        dlg_title = 'Enter # of Frames to save';
        num_lines = 1;
        value = inputdlg(prompt,dlg_title,num_lines);
        a = str2double(value{1});
        assignin('base', 'NumFrames2Sav', a);
        assignin('base','frameNo',1);
        % initiate temp matrix
        r = evalin('base', 'Us_zAxisNum');
        c = evalin('base', 'Us_xAxisNum');
%         keyboard
        M = r;
        N = c*a;
        
        evalin('base', ['tempFrameUs=zeros(', num2str(r),',',num2str(N),',''single'');']);
        evalin('base', ['tempFramePa=zeros(', num2str(r),',',num2str(N),',''single'');']);
        evalin('base', ['tempFrame=zeros(', num2str(r),',',num2str(c),',''single'');']);
        
        
        % initiate temp IQ frames
       evalin('base', ['tempFrameIQ=complex(zeros(', num2str(r),',',num2str(c),'));']);
       evalin('base', ['tempFrameUsIQdata=complex(zeros(', num2str(r),',',num2str(N),'));']);
       evalin('base', ['tempFramePaIQdata=complex(zeros(', num2str(r),',',num2str(N),'));']);
       
%         keyboard
    end
%-UI#9Callback

%-UI#10Callback - TX focus changel
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.txFocus'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P = evalin('base','P');
P.txFocus = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.txFocus = UIValue*scaleToWvl;  
    end
end
assignin('base','P',P);
TX = evalin('base', 'TX');
for n = 1:128   % 128 transmit events
    TX(n).focus = P.txFocus;
    TX(n).Delay = computeTXDelays(TX(n));
end
assignin('base','TX', TX);
txFNum = evalin('base', 'txFNum'); % get f-number value for transmit (range: 1.0 - 20)
% txNumEl is the number of elements to include on each side of the center element, for the specified
%    focus and sensitivity cutoff.  Thus the full transmit aperture will be 2*txNumEl + 1 elements.
txNumEl = round((P.txFocus/txFNum)/Trans.spacing/2); % no. of elements in 1/2 aperture.
if txNumEl > (Trans.numelements/2 - 1), txNumEl = floor(Trans.numelements/2 - 1); end  
assignin('base','txNumEl',txNumEl);
% - Redefine event specific TX attributes for the new focus.
TX = evalin('base', 'TX');
nr = evalin('base',  'nr');
% SFormat = evalin('base', 'SFormat');
% - Set event specific TX attributes.
for n = 1:nr   % 128 transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin =[(-63.5 + (n-1))*Trans.spacing, 0.0, 0.0];
    % Set Apod vector back to all zeros
    TX(n).Apod(:) = 0;
    % write new focus value to TX
    TX(n).focus = P.txFocus;
    % Set transmit Apodization so (1 + 2*TXnumel) transmitters are active.
    nCenter = (n-1) * round(128/nr) + round(128/nr/2);
    lft = max(nCenter - txNumEl, 1);
    rt = min(nCenter + txNumEl, 128);
    TX(n).Apod(lft:rt) = 1.0;
    TX(n).Delay = computeTXDelays(TX(n));
end
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%-UI#10Callback

%EF#1
stageMov()
% check if flagSav is on, and all fames neeeded is acquired
% if yes, save data using UIsave 
flagSav = evalin('base', 'flagSav');
frameNo = evalin('base','frameNo');
NumFrames2Sav = evalin('base','NumFrames2Sav');
if flagSav &&  frameNo > NumFrames2Sav
   timestamp=(datestr(clock,0));
   ind = regexp(timestamp, ':');
   timestamp(ind) = '-'; 
%    keyboard
   dname = uigetdir('C:\Users\verasonics\Desktop');
   
   paramFname = [dname, '\', 'param.mat'];
   UsFname = [dname,'\', 'Us.mat'];
   PaFname = [dname, '\', 'Pa.mat'];
   UsIQFname = [dname, '\', 'UsIQ.mat'];
   PaIQFname = [dname, '\', 'PaIQ.mat'];
   
   assignin('base', 'paramFname', paramFname);
   assignin('base', 'UsFname', UsFname);
   assignin('base', 'PaFname', PaFname);
   assignin('base', 'UsIQFname', UsIQFname);
   assignin('base', 'PaIQFname', PaIQFname);
   
%    keyboard
   evalin('base', 'param = struct(''Trans'', Trans, ''Receive'', Receive,''TX'', TX, ''TW'', TW,  ''scanAxisX'', scanAxisX, ''scanAxisY'', scanAxisY, ''Us_zAxis'', Us_zAxis, ''Us_xAxis'', Us_xAxis, ''NumFrames2Sav'', NumFrames2Sav);');
   evalin('base', 'save(paramFname, ''param'')');
   disp('param saved.');
   
%    evalin('base', 'save(UsFname, ''tempFrameUs'')');
   tic; evalin('base', 'savefast(UsFname, ''tempFrameUs'')'); toc
   disp('US frames data saved.');
   
%    evalin('base', 'save(PaFname, ''tempFramePa'')');
   tic; evalin('base', 'savefast(PaFname, ''tempFramePa'')'); toc
   disp('PA frames data saved.');
   
%    evalin('base', 'save(UsIQFname, ''tempFrameUsIQdata'', ''-v7.3'')');
   tic; evalin('base', 'savefast(UsIQFname, ''tempFrameUsIQdata'')'); toc
   disp('US frames IQ data saved.');
   
%    evalin('base', 'save(PaIQFname, ''tempFramePaIQdata'', ''-v7.3'')');
   tic; evalin('base', 'savefast(PaIQFname, ''tempFramePaIQdata'')'); toc
   disp('PA frames IQ data saved.');
   
%    evalin('base', 'param = struct(''Trans'', Trans, ''Receive'', Receive, ''SFormat'', SFormat,''TX'', TX, ''TW'', TW);');
%    evalin('base', 'uisave({''param'',''NumFrames2Sav'',''tempFrameUs'', ''tempFramePa'', ''tempFrameUsIQdata'', ''tempFramePaIQdata'', ''scanAxisX'', ''scanAxisY'', ''Us_zAxis'',''Us_xAxis''}, ''FramesTosav'')');
   % set flagSav and NumFrames2Sav to 0
   assignin('base', 'flagSav', 0);
   assignin('base', 'NumFrames2Sav', 0);
   
   scanflag = evalin('base', 'scanFlag');
   if scanflag
       comX = evalin('base', 'comX');
       comY = evalin('base', 'comY');
       
       fclose(comX); fclose(comY);
       fopen(comX); fopen(comY);
       cPosX = zaberPos(comX);  cPosY = zaberPos(comY);
%        keyboard
       xSt = evalin('base', 'scanAxisX(1)');
       ySt = evalin('base', 'scanAxisY(1)');
       zaberMove_rel_wthrespnd(xSt-cPosX, 0 , comX);
       pause(0.2);
       zaberMove_rel_wthrespnd(ySt-cPosY, 0 , comY);
       pause(0.2);
       assignin('base', 'scanFlag', 0);
%        keyboard
   end
   
end

% check if needed to move stage; if yes, move
scanflag = evalin('base', 'scanFlag');
% keyboard
if scanflag 
%    keyboard 
    mappingMatrix = evalin('base', 'mappingMatrix');
    comX = evalin('base', 'comX');
    comY = evalin('base', 'comY');
    
    [finish, success, comX, comY, posX, posY] = stage2D_Mapping(mappingMatrix, frameNo, comX, comY, 0, 0);
%     keyboard
    pause(1e-3)
%     disp('stage moving...')
%     disp(['nextFrame', num2str(frameNo)])
%     disp([posX, posY])
    assignin('base', 'stageCurPosX', posX);
    assignin('base', 'stageCurPosX', posY);
end
% toc
return
%EF#1

%EF#2
UsDisplaySave(RData)
% keyboard
myHandleUS = [];
% persistent myHandleUS
% disp('US-display --- ')
% tic
r = evalin('base', 'Us_zAxisNum');
c = evalin('base', 'Us_xAxisNum');
rAxis = evalin('base', 'Us_zAxis');
cAxis = evalin('base', 'Us_xAxis');
UsMat = RData(1:r,1:c);
% keyboard
flagSav = evalin('base', 'flagSav');
if ~flagSav    
    if isempty(myHandleUS)||~ishandle(myHandleUS)
    figure('Name', 'US');
    myHandleUS = imagesc(cAxis, rAxis, (log(UsMat/max(max(UsMat)))), 'EraseMode', 'none');
    colormap gray(256);
    colorbar;
    axis equal;
    axis tight;
    end
    set(myHandleUS, 'CData', (sqrt(UsMat/max(max(UsMat)))));
    drawnow;     
end


frameNo = evalin('base', 'frameNo');
NumFrames2Sav = evalin('base', 'NumFrames2Sav');
% tic
if flagSav && frameNo <= NumFrames2Sav
%     keyboard
%     disp('US sav')
    temp = single(UsMat);
%     stageP = evalin('base', 'stageCurPos'); 
%     disp('US sav')
    assignin('base', 'tempFrame', temp );
    dest_offset = r*c*(frameNo-1);
    evalin('base', ['MatrixOffsetAppend(tempFrameUs, tempFrame,', num2str(dest_offset),' );']);
end
% toc
return
%EF#2

%EF#3
PaDisplaySave(RData)
myHandlePA = [];
% persistent myHandlePA
% disp('PA-display --- ')
% tic
r = evalin('base', 'Us_zAxisNum');
c = evalin('base', 'Us_xAxisNum');
rAxis = evalin('base', 'Us_zAxis');
cAxis = evalin('base', 'Us_xAxis');
PAMat = RData(1:r,1:c);
% keyboard
flagSav = evalin('base', 'flagSav');
if ~ flagSav
    
    if isempty(myHandlePA)||~ishandle(myHandlePA)
        figure('Name', 'PA');
        myHandlePA = imagesc(cAxis, rAxis, log(PAMat/max(max(PAMat))), 'EraseMode', 'none');
        colormap Hot;
        colorbar;
        axis equal;
        axis tight;
    end
    set(myHandlePA, 'CData', sqrt(PAMat/max(max(PAMat))));
    drawnow;
    
end

% tic
% flagSav = evalin('base', 'flagSav');
frameNo = evalin('base', 'frameNo');
NumFrames2Sav = evalin('base', 'NumFrames2Sav');
if flagSav && frameNo <= NumFrames2Sav
%     keyboard
    temp = single(PAMat); 
%     stageP = evalin('base', 'stageCurPos');
%     disp([stageP, frameNo])
%     disp('PA sav')
    assignin('base', 'tempFrame', temp );
    dest_offset = r*c*(frameNo-1);
    evalin('base', ['MatrixOffsetAppend(tempFramePa, tempFrame,', num2str(dest_offset),' );']);
%     assignin('base','frameNo', frameNo+1);    
end
% toc
return
%EF#3

%EF#4
UsIQSave(RData)
% tic
% keyboard
r = evalin('base', 'Us_zAxisNum');
c = evalin('base', 'Us_xAxisNum');
rAxis = evalin('base', 'Us_zAxis');
cAxis = evalin('base', 'Us_xAxis');
% tic
UsMat = RData(1:r,1:c);
% toc
flagSav = evalin('base', 'flagSav');
frameNo = evalin('base', 'frameNo');
NumFrames2Sav = evalin('base', 'NumFrames2Sav');
% tic
if flagSav && frameNo <= NumFrames2Sav
%     keyboard
    temp = (UsMat);
%     stageP = evalin('base', 'stageCurPos'); 
%     disp([stageP , frameNo])
%     disp('US IQ sav')
    assignin('base', 'tempFrameIQ', temp );
    dest_offset = r*c*(frameNo-1);
    evalin('base', ['CmplxMatrixOffsetAppend(tempFrameUsIQdata, tempFrameIQ,', num2str(dest_offset),' );']);
end
% toc
return
%EF#4

%EF#5
PaIQSave(RData)
% tic
r = evalin('base', 'Us_zAxisNum');
c = evalin('base', 'Us_xAxisNum');
rAxis = evalin('base', 'Us_zAxis');
cAxis = evalin('base', 'Us_xAxis');
PAMat = RData(1:r,1:c);

flagSav = evalin('base', 'flagSav');
frameNo = evalin('base', 'frameNo');
NumFrames2Sav = evalin('base', 'NumFrames2Sav');
if flagSav && frameNo <= NumFrames2Sav
%     keyboard
    temp = (PAMat); 
%     stageP = evalin('base', 'stageCurPos');
%     disp([stageP, frameNo])
%     disp('PA IQ sav')
    assignin('base', 'tempFrameIQ', temp );
    dest_offset = r*c*(frameNo-1);
%     keyboard
    evalin('base', ['CmplxMatrixOffsetAppend(tempFramePaIQdata, tempFrameIQ,', num2str(dest_offset),' );']);
    assignin('base','frameNo', frameNo+1);    
end
% toc
return
%EF#5
