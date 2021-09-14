% File name: Sequence_LiveView_and_save_uAM.m - Ultrafast multi-angle plane wave imaging 
%                                               with and amplitude modulation paradigm
%                                          
% Description: 
%   Sequence programming file for L22-4v Linear array, using plane wave transmits
%   with multiple steering angles combined in a single pulse code for each channel. 
%   The sign is changed on pulse within each angle, so that by combining all the angle
%   acquisitions before reconstruction, the individual angle tranmit information can be
%   recovered. 
%   The mutliplane wave transmssion is combined with an amplitude
%   modulation paradigm to obtain the nonlinear signal from the medium.
%   Processing is asynchronous with respect to acquisition.
%
% Rabut et al., Appl. Phys. Lett. 118, 244102 (2021); https://doi.org/10.1063/5.0050807


clear all

P.pathName = 'D:\YourPath' ; 
P.numberImagesToSave = 10 ;

P.startDepth_mm = 2.5;   
P.endDepth_mm = 10;   


% Aperture
P.imaging_aperture = 90 ; % Determine imaging aperture
Aperture(1) = (128-P.imaging_aperture)/2 + 1 ;
Aperture(2) = 128 - (128-P.imaging_aperture)/2 ;

P.Apod = zeros(1,128);   % From aperture, calculate Apodization
for x = Aperture(1): Aperture(2)
P.Apod(1,x) = 1;
end

DisplayScaling = 'power' ; % if ‘power’, 20=linear, 40 = sqrt - if ‘log’, 40 = 40dB dyn. range
Display_output = 20 ; %Display mode only: 20 gives a power of 1.0 for a linear output, 40 gives a power of 0.5 for a square root compression


na = 8;      % Set na = number of angles (2, 4 or 8).
if (na > 1), dtheta = (14*pi/180)/(na-1); P.startAngle = -14*pi/180/2; else dtheta = 0; P.startAngle=0; end % set dtheta to range over +/- 7 degrees.


% Define system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.simulateMode = 0;
% Resource.Parameters.simulateMode = 1 ;% forces simulate mode, even if hardware is present.

% If simulation 
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';


% Specify Trans structure array
Trans.name = 'L22-14v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans.maxHighVoltage = 25;  % set maximum high voltage limit for pulser supply.
Trans = computeTransMG(Trans);  % L11-4v transducer is 'known' transducer so we can use computeTrans.


% Convert mm to wavelength
demodFreq = Trans.frequency; % demodulation frequency
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.startDepth = P.startDepth_mm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepth_mm*scaleToWvl;      % Acquisition depth in wavelengths



% Specify PData structure array.
PData(1).PDelta =  [Trans.spacing/2, 0, Trans.frequency/Resource.Parameters.speedOfSound];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((P.imaging_aperture*Trans.spacing)/PData.PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(P.imaging_aperture-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.
P.lateral_voxel_size = PData(1).PDelta(1)*100e-6 ;
P.axial_voxel_size = PData(1).PDelta(3)*100e-6 ;

% - specify Region structures.
PData(1).Region = struct('Shape',struct( ...
    'Name','Rectangle',...
    'Position',[0,0,P.startDepth],...
    'width',Trans.spacing*P.imaging_aperture,...
    'height',P.endDepth-P.startDepth )) ;
PData(1).Region = computeRegions(PData(1));


%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = na*4096; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = 128;
Resource.RcvBuffer(1).numFrames = 2; 


% - Define a 2nd RcvBuffer to receive the process RF data from RcvBuffer 1.
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = na*4096; % this size allows for maximum range
Resource.RcvBuffer(2).colsPerFrame = 128;
Resource.RcvBuffer(2).numFrames = 1;    % 1 frame RcvBuffer for processed acquisitions.


Resource.InterBuffer(1).numFrames = 1; % one intermediate buffer needed.

Resource.ImageBuffer(1).numFrames = 5;

Resource.DisplayWindow(1).Title = 'uAM';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).numFrames = 5;
Resource.DisplayWindow(1).Colormap = hot(256);

% Specify Transmit waveform structure. 
% - Define a normal and inverted parametric TW to get the pulse code values. 
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];   % A, B, C, D
% - Inverse polarity.
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,.67,2,-1];   % A, B, C, D
% Call computeTWWaveform here to compute pulse codes for TW 1&2.
[~, ~, ~, rc, TW] = computeTWWaveform(TW);
TW(1).PulseCode= [1 -3 4 5 1 ; 3 -5 4 3 1 ];
TW(2).PulseCode= [1 3 4 -5 1 ; 3  5  4 -3 1] ;
pcRows = size(TW(1).PulseCode,1); % no. of rows in pulse code for waveform

% - Define na additional TWs which will receive the calculated pulse codes.
for j=1:na
    TW(2+j).type = 'pulseCode';
    TW(2+j).PulseCode = [];
end

% Specify TX structure array.
%    - Define na TXs that can be used for specifying the single transmit delays for each angle.
%      These will be used for reconstruction of the process RF data --> The delay compensation is thus directly take into account!
%    - Define a second set of 3*na TXs that will be used with the multiwave TWs for acquisition
%       - TX(1) --> TX(na) : definition of the transmit only for the reconstruction
%       - TX(na+1)   -->   TX(2*na) : na transmit of na successive plane waves of half-amplitude (pulse 1 of AM-imaging)
%       - TX(2*na+1) -->   TX(3*na) : na transmit of na successive plane waves of half-amplitude (pulse 2 of AM-imaging)
%       - TX(3*na+1) -->   TX(4*na) : na transmit of na successive plane waves of full-amplitude (pulse 3 of AM-imaging)
%    ==> So in total 4*na TX structures

TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', ones(1,128), ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, 4*na); 
               
% - Set TX attributes for normal flash angle scan in first na TXs.
if fix(na/2) == na/2       % if na evenI
    P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(na/2)*dtheta;
end

% - Calculate TX.Delays so that there is no overlap in the delay times. Add 3 wls to separate waveforms
for n = 1:na   % na transmit events
    TX(n).Steer = [(P.startAngle+(n-1)*dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
    if n~=1, TX(n).Delay = TX(n).Delay + (TX(n-1).Delay(1)-TX(n).Delay(1)) + 3; end
end


% Specify the TW waveforms for the acquisition TXs

k = na;
for p = 1:3
    for n = 1:na        
        TX(k+na*(p-1)+n).waveform = 2+n;      % Same Waveform combination for the three set of pulses
    end
end
% Different apodistion for the three pulses
for n = 1:na
    TX(k+n).Apod = P.Apod ;
    TX(k+n).Apod(1:2:128) = 0;   % Pulse 1
    TX(k+na+n).Apod(2:2:128) = 0;   % Pulse 2
end


% Calculate the pulse codes for the acquisition TWs
% - Define Hadamard matrix (n=8)
H = hadamard(na) ; 
  
  
% - Compose TW pulseCode definitions for each angle that include pulses for all angles.  
%   Use the Hadamard row values to determine the sign of the transmit waveform.
wls2clks = (1/Trans.frequency)/(1/250); % no. of clks in one wavelength
for i= 1:na   % compute the pulse codes for all angles
    TW(2+i).PulseCode = zeros(na*pcRows,5,128); % Initialize pulseCode array.
    SumClks = zeros(1,128);  % SumClks has the total clks for each channel
    for j=1:na
        rowstart = (j-1)*pcRows + 1; % starting row for pulse code for angle
        rowend = rowstart + pcRows - 1;
        % Get sign of pulse from Hadamard array
        if H(i,j)==1 
            PCode = TW(1).PulseCode; % use positive code
        else
            PCode = TW(2).PulseCode;  % use negative code
        end
        for k = Aperture(1) :Aperture(2)
            % Calculate pulse code for channel, based on delay.
            TW(2+i).PulseCode(rowstart:rowend,:,k) = PCode;
            TW(2+i).PulseCode(rowstart,1,k) = PCode(1,1) + round(TX(j).Delay(k)*wls2clks) - SumClks(k);
            % Add to SumClks the clks used in the added pulse code
            SumCols = sum(abs(TW(2+i).PulseCode(rowstart:rowend,1:4,k)));
            SumClks(k) = SumClks(k) + sum(SumCols);
        end
    end
end
        
%% Specify TGC Waveform structure.
TGC.CntrlPts =  1023*ones(1,8); 
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays. 
% - We need na Receives for the acquisition frames and one frame for RcvBuffer 2, which will
%   be used for reconstruction.
maxAcqLength = ceil(sqrt(P(1).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', P.Apod, ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 3*na*Resource.RcvBuffer(1).numFrames+na);
                    
                    
% - Set event specific Receive attributes for each frame

Apod_pulse_1_2 = -1 ;
if Apod_pulse_1_2 == 1
    display('!!!!!!!'),display('WARNING: Linear mode!!')
    display('Change value of [Apod_pulse_1_2] to -1 for AM mode'),display('!!!!!!!')
end

for i = 1:Resource.RcvBuffer(1).numFrames
    
    Receive(na*(i-1)+1).callMediaFunc = 1;
    
    for p = 1:3
            
            
        for j = 1:na
            
            if p==1
                Receive(3*na*(i-1) + na*(p-1) + j).Apod(1:128) = P.Apod*Apod_pulse_1_2;   % Pulses 1 (Apod in reception = -1 )
            end
            
            if p==2
                Receive(3*na*(i-1) + na*(p-1) + j).mode = 1;          % Summation AM in the hardware
                Receive(3*na*(i-1) + na*(p-1) + j).Apod(1:128) = P.Apod*Apod_pulse_1_2;  % Pulses 2 (Apod in reception = -1 )
            end
            
            if p==3
                Receive(3*na*(i-1) + na*(p-1) + j).mode = 1;          % Summation AM in the hardware, Pulses 3 (Apod in reception = +1 )
            end
            
            Receive(3*na*(i-1) + na*(p-1) + j).framenum = i;
            Receive(3*na*(i-1) + na*(p-1) + j).acqNum = j;
            
        end
        
    end
end

% Define the last frame of na Receives for reconstruction from RcvBuffer 2.
k = 3*na*Resource.RcvBuffer(1).numFrames;
for j = 1:na
    Receive(k+j).bufnum = 2;
    Receive(k+j).acqNum = j;
end

%% Specify Recon structure arrays.
% - We need one Recon structures which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
    'pdatanum', 1, ...
    'IntBufDest', [1,1], ...
    'ImgBufDest', [1,-1], ...
    'RINums', 1:na);

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % Compounding done in hardware
    'txnum', 1, ...
    'rcvnum', 1, ...
    'regionnum', 1), 1, na);
% - Set specific ReconInfo attributes.
k = 3*na*Resource.RcvBuffer(1).numFrames;

ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
ReconInfo(2).mode = 'replaceIQ'; % replace IQ data: selected compounding, do not take the first hadamard summation
ReconInfo(3).mode = 'replaceIQ'; % replace IQ data: selected compounding, do not take the second hadamard summation as well


for j = 1:na  % For each angle
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = k+j;
end
ReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect


%% Specify Process structure arrays.
pers = 20;                   
                     
Process(1).classname = 'External';
Process(1).method = 'processRF';
Process(1).Parameters = {'srcbuffer','receive',...  % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...
                         'dstbuffer','receive',...
                         'dstbufnum',2,...
                         'dstframenum',1};
EF(1).Function = text2cell('%#EF-1');

                     
% Save Image Data
Process(2).classname = 'External';     % process structure for 1st Doppler ensemble
Process(2).method = 'saveImageData';
Process(2).Parameters = {'srcbuffer','image',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...            % process the most recent frame.
                         'dstbuffer','none'};                     
EF(2).Function = text2cell('%#EF-2');
                     

Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',2,...            % pgain is image processing gain
                         'reject',2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod',DisplayScaling,...
                         'compressFactor',Display_output,...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
                     
                     
                     
% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 160;  % 160 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000 - (na-1)*160; % in usec % 20 msec
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'loopCnt'; % - Set loop count. for looping on blocs
SeqControl(5).argument = P.numberImagesToSave-1;  %
SeqControl(5).condition = 'counter1';
SeqControl(6).command = 'loopTst';  % loop test
SeqControl(6).argument = [];    % set apres
SeqControl(6).condition = 'counter1';
SeqControl(7).command = 'timeToNextAcq';  % time between frames
SeqControl(7).argument = 2e6; % 
nsc = 8; % nsc is count of SeqControl objects

%% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for p = 1:3
        for j = 1:na                      % Acquire frame
            Event(n).info = 'Acquisition';
            Event(n).tx = na+na*(p-1)+j;   % use acquisition TX structures.
            Event(n).rcv = 3*na*(i-1) + na*(p-1) + j;
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 2;
            n = n+1;
        end
    end
    
    Event(n-1).seqControl = 3; % modify last acquisition Event's seqControl
    
    Event(n).info = 'Transfer To Host';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    n = n+1;
    
    Event(n).info = 'Wait for transfer complete';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'waitForTransferComplete'; % transfer frame to host buffer
    SeqControl(nsc).argument = nsc-1; % reference previous transferToHost
    nsc = nsc+1;
    n = n+1;
    
    Event(n).info = 'Hadamard summing';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'markTransferProcessed';
    SeqControl(nsc).argument = nsc-2; % reference previous transferToHost
    nsc = nsc+1;
    n = n+1;
    
    Event(n).info = 'recon and display image';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 3;    % process
    Event(n).seqControl = 4; % return to Matlab
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = 1;
n = n+1;

% Acquire Image Data
nStartAcquire = n;

% set loop count
Event(n).info = 'start counter';
Event(n).tx = 0;   % use next TX structure.
Event(n).rcv = 0;
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = 5;
n = n+1;
SeqControl(6).argument = n;

for i = 1:Resource.RcvBuffer(1).numFrames
    for p = 1:3
        for j = 1:na                      % Acquire frame
            Event(n).info = 'Acquisition';
            Event(n).tx = na+na*(p-1)+j;   % use acquisition TX structures.
            Event(n).rcv = 3*na*(i-1) + na*(p-1) + j;
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 2;
            n = n+1;
        end
    end
    
    Event(n-1).seqControl = 7; % modify last acquisition Event's seqControl
    
    Event(n).info = 'Transfer To Host';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    n = n+1;
    
    Event(n).info = 'Wait for transfer complete';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'waitForTransferComplete'; % transfer frame to host buffer
    SeqControl(nsc).argument = nsc-1; % reference previous transferToHost
    nsc = nsc+1;
    n = n+1;
    
    Event(n).info = 'Hadamard summing';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'markTransferProcessed';
    SeqControl(nsc).argument = nsc-2; % reference previous transferToHost
    nsc = nsc+1;
    n = n+1;
    
    Event(n).info = 'recon and display image';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 3;    % process
    Event(n).seqControl = 4; % return to Matlab
    n = n+1;
    
end


Event(n).info = 'save Image Data'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 2;    % process
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Dummy transmit to impose a break';
Event(n).tx = 1;  
Event(n).rcv = 0;
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = 7;



Event(n).info = 'Loop back until all images are acquired';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % process
Event(n).seqControl = 6 ;
n = n+1;



% UI controls
% - Acquire RF data
UI(1).Control = {'UserC1','Style','VsPushButton','Label','Save Image'};
UI(1).Callback = text2cell('%StartSaving');


% Save all the structures to a .mat file.
save('MatFiles/Sequence_LiveView_and_save_uAM');
filename = 'Sequence_LiveView_and_save_uAM';
VSX
return



%StartSaving
display('Start Saving')
Control = repmat(struct('Command','set&Run','Parameters',[]),1,1);
nStartAcquire = evalin('base','nStartAcquire');
Control(1).Parameters = {'Parameters',1,'startEvent',nStartAcquire};
evalin('base','Resource.Parameters.startEvent = nStartAcquire;');
assignin('base','Control', Control);
%StartSaving


%% External Functions definition

%#EF-1
RFOut = processRF(RData)
% For speed, this function could be implemented as a mex function written in C.
persistent HA
if isempty(HA), HA = evalin('base','H'); end

RFOut = zeros(size(RData,1),size(RData,2));  % make zeroed copy of RData
RFOut = int16(RFOut) ;

Receive = evalin('base', 'Receive');
na = evalin('base','na');

% Determine startSample and endSample for the various acqNums.
ss = zeros(1,na);
se = zeros(1,na);
for j=1:na
    ss(j) = Receive(j).startSample;
    es(j) = Receive(j).endSample;
end
for i=1:na
    for j=1:na
        RFOut(ss(i):es(i),:) = RFOut(ss(i):es(i),:) + HA(j,i)*RData(ss(j):es(j),:);
    end
end

return
%#EF-1


%#EF-2
saveImageData(ImgData)

persistent bloc_count;
   if isempty(bloc_count)
      bloc_count = 1;
   end
P = evalin('base','P');
Path = P.pathName;

if bloc_count == 1 
   mkdir(Path)
   save([Path '\P.mat'],'P');
end
eval(['uAM_',num2str(bloc_count),' = squeeze(ImgData);']);
save([Path '\uAM_' num2str(bloc_count) '.mat'],['uAM_',num2str(bloc_count)]);
display(['Saved image #',num2str(bloc_count) '!'])
bloc_count = bloc_count+1;
if bloc_count == P.numberImagesToSave + 1
    display('End of the acquisition!')
end

%#EF-2
    
 
