clear all
close all

%% initialization
frameNo = 0;
tempFrame = [];
tempFrameUs = [];
tempFramePa = [];
UsFileName = [];
PaFileName = [];

savCurrentFrameFlag = 0;

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

%%
P.startDepth = 5;
P.endDepth = 210;
P.txFocus = 75;  % Initial transmit focus.
P.numTx = 32;  % number of transmit elements in TX aperture (where possible).

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1480;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L7-4';
% Trans.name = 'L22-14v';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L7-4 transducer is 'known' transducer so we can use computeTrans.
% note nominal center frequency from computeTrans is 6.25 MHz
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
% Trans.maxHighVoltage = 30;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];  % x, y, z pdeltas
PData(1).pdeltaX = Trans.spacing;
PData(1).pdeltaZ = 0.5;
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
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

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*128;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;       % 20 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L7-4_128RyLns';
Resource.DisplayWindow(1).pdelta = 0.7;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify TW structure array.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% Specify nr TX structure arrays. Transmit centered on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 128);          
% - Set event specific TX attributes.
for n = 1:128   % 128 transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = [(-63.5 + (n-1))*Trans.spacing, 0.0, 0.0];
    % Set transmit Apodization.
    lft = n - floor(P.numTx/2);
    if lft < 1, lft = 1; end;
    rt = n + floor(P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end;
    TX(n).Apod(lft:rt) = 1.0;
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify Receive structure arrays. 
% - We need 128 Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wlsPer128 = 128/(2*2); % wavelengths in 128 samples for 2 samplesPerWave
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW',...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 128*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 128*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:128
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,138,260,287,385,593,674,810];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure array. 
Recon = struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:128);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, 128);
% - Set specific ReconInfo attributes.
for j = 1:128 
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
                     
                     

% External Processing - US Display and svae Frame - Lu.
Process(2).classname = 'External';
Process(2).method = 'ImgSave';
Process(2).Parameters = {'srcbuffer','image',... % name of buffer to process.
    'srcbufnum',1,...
    'srcframenum',-1,... % process the most recent frame.
    'dstbuffer','none'};


% Specify SeqControl structure arrays.
% add ext trigger control here
% input trigger 
SeqControl(1).command = 'noop';
SeqControl(1).argument = 0; % Trigger input 1, enable with rising edge 
% noop delay between trigger in and start of acquisition
SeqControl(2).command = 'noop';
SeqControl(2).argument = fix(0)*5; % noop counts are in 0.2 microsec increments
% SeqControl(13).argument = 3*5; % noop counts are in 0.2 microsec increments
% SeqControl(13).argument = 1e3; % noop counts are in 0.2 microsec increments
% sync command
SeqControl(3).command = 'sync';
SeqControl(3).argument = 2e9;   % timeout for waiting ext trigger

%  - Time between acquisitions in usec
t1 = round(2*384*(1/Trans.frequency)); % acq. time in usec for max depth
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = t1;
%  - Time between frames at 20 fps at max endDepth.
SeqControl(5).command = 'timeToNextAcq';
SeqControl(5).argument = 50000; 
%  - Return to Matlab
SeqControl(6).command = 'returnToMatlab';
%  - Jump back to start.
SeqControl(7).command = 'jump';
SeqControl(7).argument = 1;
nsc = 8; % next SeqControl number

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    
    % Wait for input trigger from flash lamp firing
    Event(n).info = 'Start Acquisition';
    Event(n).tx = 1;        % no TX
    Event(n).rcv = 1;       % no Rcv
    Event(n).recon = 0;     % no Recon
    Event(n).process = 0;   % no Process
    Event(n).seqControl = 1; % Trigger in
    n = n+1;
    
    % Pause for optical buildup
    Event(n).info = 'noop and sync';
    Event(n).tx = 0;        % no TX
    Event(n).rcv = 0;       % no Rcv
    Event(n).recon = 0;     % no Recon
    Event(n).process = 0;   % no Process
    Event(n).seqControl = [2,3]; % pause and sync
    n = n+1;

    for j = 1:128                      % Acquire frame
        Event(n).info = 'Aqcuisition.';
        Event(n).tx = j;   % use next TX structure.
        Event(n).rcv = 128*(i-1)+j;   
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 4; % seqCntrl
        n = n+1;
    end
    % Replace last events SeqControl for inter-frame timeToNextAcq.
   Event(n-1).seqControl = 5;
    
    Event(n).info = 'Transfer frame to host.';
    Event(n).tx = 0;        % no TX
    Event(n).rcv = 0;       % no Rcv
    Event(n).recon = 0;     % no Recon
    Event(n).process = 0; 
    Event(n).seqControl = nsc; 
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;
    n = n+1;

    Event(n).info = 'recon and process'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 1;    % process
    Event(n).seqControl = 0;
    n = n + 1;
    
    %  Display & Sav   
    Event(n).info = 'US sav'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;  % reconstruction for both 2D and PA
    Event(n).process = 2;    % process 2D
    Event(n).seqControl = 0; % no seqCntrl
    
    if (floor(i/5) == i/5)&&(i ~= Resource.RcvBuffer(1).numFrames)  % Exit to Matlab every 5th frame
        Event(n).seqControl = 6; % return to Matlab
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = 7;


% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');
             
% - Transmit focus change
UI(3).Control = {'UserB4','Style','VsSlider','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[20,320,P.txFocus]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%TxFocusCallback');
             
% - F number change
UI(4).Control = {'UserB3','Style','VsSlider','Label','F Number',...
                 'SliderMinMaxVal',[1,20,round(P.txFocus/(P.numTx*Trans.spacing))],'SliderStep',[0.05,0.1],'ValueFormat','%2.0f'};
UI(4).Callback = text2cell('%FNumCallback');

% - UI save push button for current frame - Lu
UI(5).Control = {'UserC3','Style','VsPushButton','Label','Save CurrentFrame'};
UI(5).Callback = text2cell('%-UI#5Callback');

%% PA Recon Parameter Cal and Set
NoElem = 64;
% [ImgDepthAxis, xAxis, zAxis, ApertureAxis] = myPA_ReconParameterCal_v307(NoElem);
% calculate US Display Axis
[Us_zAxis, Us_xAxis, Us_zAxisNum, Us_xAxisNum] = myUS_DisplayParameterCal_v307();

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;

%% Call External Control or display function
EF(1).Function = text2cell('%EF#1');        % stage scanning control

%% Save all the structures to a .mat file.

filename = ('L7-4vFlashExtTrigUS0206');   % define variable 'filename' to permit VSX to skip user query for matfile
save (filename)                 % save the structures to a matfile
VSX                             % invoke VSX automatically when running this Setup script


return

% **** Callback routines to be converted by text2cell function. ****
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

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;    
    end
end
assignin('base','P',P);

PData = evalin('base','PData');
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','Rectangle',...
                    'Position',[0,0,P.startDepth],...
                    'width',Trans.spacing,...
                    'height',P.endDepth-P.startDepth)),1,128);
% - set position of regions to correspond to beam spacing.
for i = 1:128
    PData(1).Region(i).Shape.Position(1) = (-63.5 + (i-1))*Trans.spacing;
end
assignin('base','PData',PData);
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
evalin('base','if VDAS==1, Result = loadTgcWaveform(1); end');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%TxFocusCallback - TX focus changel
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
% Update Fnumber based on new P.txFocus
evalin('base','set(UI(4).handle(2),''Value'',round(P.txFocus/(P.numTx*Trans.spacing)));');
evalin('base','set(UI(4).handle(3),''String'',num2str(round(P.txFocus/(P.numTx*Trans.spacing))));');
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%TxFocusCallback

%FNumCallback - F number change
simMode = evalin('base','Resource.Parameters.simulateMode');
P = evalin('base','P');
Trans = evalin('base','Trans');
% No F number change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',round(P.txFocus/(P.numTx*Trans.spacing)));
    return
end
P.txFNum = UIValue;
P.numTx = round(P.txFocus/(P.txFNum*Trans.spacing));
assignin('base','P',P);
% - Redefine event specific TX attributes for the new P.numTx.
TX = evalin('base', 'TX');
for n = 1:128   % 128 transmit events
    % Set transmit Apodization.
    lft = n - floor(P.numTx/2);
    if lft < 1, lft = 1; end;
    rt = n + floor(P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end;
    TX(n).Apod = zeros(1,Trans.numelements);
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
%FNumCallback

%-UI#5Callback - save current frame
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
%-UI#5Callback

%EF#1
ImgSave(RData)
% keyboard
% persistent myHandleUS
% disp('US-display --- ')
% tic
% keyboard

flagSav = evalin('base', 'flagSav');
frameNo = evalin('base','frameNo');
NumFrames2Sav = evalin('base','NumFrames2Sav');

% tic
if flagSav
%     keyboard
%     disp('US sav')
    r = evalin('base', 'Us_zAxisNum');
    c = evalin('base', 'Us_xAxisNum');
    rAxis = evalin('base', 'Us_zAxis');
    cAxis = evalin('base', 'Us_xAxis');
    UsMat = RData(1:r,1:c);

    if frameNo <= NumFrames2Sav
        temp = single(UsMat);
        assignin('base', 'tempFrame', temp );
        dest_offset = r*c*(frameNo-1);
        evalin('base', ['MatrixOffsetAppend(tempFrameUs, tempFrame,', num2str(dest_offset),' );']);
        assignin('base','frameNo', frameNo+1);   
    else
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
       evalin('base', 'param = struct(''Trans'', Trans, ''Receive'', Receive,''TX'', TX, ''TW'', TW, ''Us_zAxis'', Us_zAxis, ''Us_xAxis'', Us_xAxis, ''NumFrames2Sav'', NumFrames2Sav);');
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

       % set flagSav and NumFrames2Sav to 0
       assignin('base', 'flagSav', 0);
       assignin('base', 'NumFrames2Sav', 0);
    end
    
end
% toc
return
%EF#1
