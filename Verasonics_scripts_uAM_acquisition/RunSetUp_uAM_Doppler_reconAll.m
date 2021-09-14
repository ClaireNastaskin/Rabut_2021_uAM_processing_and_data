% Sequence Combined AM and Doppler (Multiplane Wave)
% File name: RunSetUp_ultrafast_AM_MPW_Doppler_reconAll
% This second script has been written to reconstruct and displays all acquisitions. To be compatible with the acquired data,
%   it must use the same RcvBuffer size and is intended to be run in simulateMode=2 using the previously collected data. It works
%   in real-time as well, but the processing cannot keep up with acquisition, and some superframes are skipped but the script
%   displays all sub-frames within a super frame.

% Description:
%   Sequence programming for L22-14v Linear array, using multiplane plane wave
%   transmits with multiple steering angles. For each angle, three transmits
%   are sent: pulse 1 -> half-amplitude
%             pulse 2 -> half-amplitude
%             pulse 3 -> full amplitude
%
% Apod = -1 for pulse 1 and pulse 2 and Apod = 1 for pulse 3
%
% Hadamard summing: external process
% Beamforming and coherent summing of the angles during reconstruction process.
%----last update ---- C.Rabut july 2020
% All the IQs are saved

clear all

global UF na dir_save dtheta start_bloc

dir_save = 'YourFileName';
load([ dir_save, 'UF.mat'])
load([ dir_save, 'P.mat'])
load([ dir_save, 'Receive.mat'])

start_bloc = 712;

na                  = UF.numAngles ;
dtheta              = (14*pi/180)/(na-1) ;

Aperture(1) = (128-UF.imaging_aperture)/2 + 1 ;
Aperture(2) = 128 - (128-UF.imaging_aperture)/2 ;


% clearvars -except P UF na dtheta filename UF dir_save
simulateMode = 2; % use mode=2 for post processing data collected by 'RunSetUp_ultrafast_AM_mpw_Doppler_acquire' and stored in 'RFdataHFR.mat'.


% Define system parameters.
filename = 'RunSetUp_uAM_Doppler_acquire'; % used to launch VSX automatically
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;
Resource.Parameters.Connector = 1;

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';


Trans.name = UF.Probe;
Trans.units = 'wavelengths';    % Explicit declaration avoids warning message when selected by default
Trans.frequency = UF.TwFreq ;   % emission frequency [MHz]
Trans = computeTransMG(Trans);  % L22-14v transducer is 'known'(added by user) transducer so we can use computeTrans.
Trans.maxHighVoltage = 26;      % set maximum high voltage limit for pulser supply.

% conversion in lambda
UF.Lambda = Resource.Parameters.speedOfSound/UF.RcvFreq*1e-3;   % [mm]
P.startDepth =  floor(UF.Depth(1)/UF.Lambda);	% Acquisition depth in wavelengths
P.endDepth =    ceil(UF.Depth(2)/UF.Lambda);    % This should preferrably be a multiple of 128 samples

% Specify PData structure array.
PData(1).PDelta =  [Trans.spacing/2, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = ceil((UF.imaging_aperture*Trans.spacing)/PData.PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(UF.imaging_aperture-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.
P.lateral_voxel_size = PData(1).PDelta(1)*100e-6 ;
P.axial_voxel_size = PData(1).PDelta(3)*100e-6 ;

% - specify Region structures.
PData(1).Region = struct('Shape',struct( ...
    'Name','Rectangle',...
    'Position',[0,0,P.startDepth],...
    'width',Trans.spacing*UF.imaging_aperture,...
    'height',P.endDepth-P.startDepth )) ;
PData(1).Region = computeRegions(PData(1));


% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = UF.rowsPerFrameRcvBuffer ;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 1;    % 1 superframe

% - Define a 2nd RcvBuffer to receive the process RF data from RcvBuffer 1.
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = UF.rowsPerFrameRcvBuffer ;
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 1;    % 1 super frame


% Interbuffer used to store IQ data
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;                      % one intermediate buffer needed.
Resource.InterBuffer(1).pagesPerFrame = 2*na*UF.numFrames;
Resource.InterBuffer(1).rowsPerFrame = PData(1).Size(1);	% this size allows for maximum range
Resource.InterBuffer(1).colsPerFrame = PData(1).Size(2);

Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 1;
Resource.ImageBuffer(1).rowsPerFrame = PData(1).Size(1);
Resource.ImageBuffer(1).colsPerFrame = PData(1).Size(2);


%% Specify Transmit waveform structure.
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
    UF.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    UF.startAngle = -fix(na/2)*dtheta;
end

% - Calculate TX.Delays so that there is no overlap in the delay times. Add 3 wls to separate waveforms
for n = 1:na   % na transmit events
    TX(n).Steer = [(UF.startAngle+(n-1)*dtheta),0.0];
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
        for k = Aperture(1) : Aperture(2)
            % Calculate pulse code for channel, based on delay.
            TW(2+i).PulseCode(rowstart:rowend,:,k) = PCode;
            TW(2+i).PulseCode(rowstart,1,k) = PCode(1,1) + round(TX(j).Delay(k)*wls2clks) - SumClks(k);
            % Add to SumClks the clks used in the added pulse code
            SumCols = sum(abs(TW(2+i).PulseCode(rowstart:rowend,1:4,k)));
            SumClks(k) = SumClks(k) + sum(SumCols);
        end
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [700 900 1000 1023 1023 1023 1023 1023];% 1000*ones(1,8);
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);



% Specify Receive structure arrays.

% For Doppler, use narrow bandwidth coefficients (50% BW) centered at
% 15.625; this is a copy of update function's default coef array for 1
% samples per wave
BPFDop = [ -0.00162 +0.00000 +0.00568 +0.00000 -0.01065 +0.00000 +0.01349 ...
    +0.00000 -0.00858 +0.00000 -0.00955 +0.00000 +0.04312 +0.00000 ...
    -0.08841 +0.00000 +0.13550 +0.00000 -0.17130 +0.00000 +0.18463 ];


% Specify Receive structure arrays.
% - We need na Receives for the acquisition frames and one frame for RcvBuffer 2, which will
%   be used for reconstruction.
maxAcqLength = ceil(sqrt(P(1).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
    'startDepth', P.startDepth, ...
    'endDepth', maxAcqLength,...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', UF.sampling_mode, ...
    'mode', 0, ...
    'callMediaFunc', 0), 1, 3*na*UF.numFrames + 2*na*UF.numFrames);


% - Set event specific Receive attributes for each frame

Apod_pulse_1_2 = -1 ;
if Apod_pulse_1_2 == 1
    display('!!!!!!!'),display('WARNING: Linear mode!!')
    display('Change value of [Apod_pulse_1_2] to -1 for AM mode'),display('!!!!!!!')
end

for i = 1:UF.numFrames
    
    Receive(na*(i-1)+1).callMediaFunc = 1;
    
    for p = 1:3
        
        
        for j = 1:na
            
            if p==1
                Receive(3*na*(i-1) + na*(p-1) + j).Apod(1:128) = Apod_pulse_1_2;   % Pulse 1 (Apod in reception = -1 )
                Receive(3*na*(i-1) + na*(p-1) + j).acqNum = 2*na*(i-1) + j; %
            end
            
            if p==2
                Receive(3*na*(i-1) + na*(p-1) + j).mode = 1;          % Pulse 2 (mode=1 because summation of the pulse done in the receive buffer directly)
                Receive(3*na*(i-1) + na*(p-1) + j).Apod(1:128) = Apod_pulse_1_2;  % (Apod in reception = -1 )
                Receive(3*na*(i-1) + na*(p-1) + j).acqNum = 2*na*(i-1) +j; %   Pulse 2 is writting on the same acquisition number than pulse 1
                
            end
            
            if p==3
                Receive(3*na*(i-1) + na*(p-1) + j).mode = 0;          % Pulse 3 (mode=0 because we remplace the previous data. We want to save separately pulse)
                Receive(3*na*(i-1) + na*(p-1) + j).acqNum = 2*na*(i-1) +na + j; %  Pulse 3 is writting on a different acquisition number than pulse 1 and pulse 2
            end
            
            
        end
        
    end
end

% Define 2*na Receives acq. for reconstruction from RcvBuffer 2: only one frame for Recon.

a = 3*na*UF.numFrames;
for f = 1:UF.numFrames
    for j = 1:2*na
        Receive(a+ 2*na*(f-1) + j).framenum = 1 ;
        Receive(a+ 2*na*(f-1) + j).bufnum = 2;
        Receive(a+ 2*na*(f-1) + j).acqNum = 2*na*(f-1) + j;
    end
end

%% _________________________________________________________________________%
%     Recon/ReconInfo for reconstruct and store the IQData of each acq      %

% Define ReconInfo structures.
% We need na ReconInfo structures for the acquisition of Pulse 1 and 2 (already summed in hardware)
% ... and na ReconInfo structures for the acquisition of Pulse 3
% ... Times:  UF.numFrames
% In total: UF.numFrames*2*na ReconInfo structures

ReconInfo = repmat(struct('mode', 'replaceIQ', ...  % default is to accumulate IQ data.
    'txnum', 1, ...
    'rcvnum', 1, ...
    'pagenum',[],...
    'regionnum', 1), 1, 2*na*UF.numFrames);


% Only one Recon because everything is stored
% in one super frame only.
Recon(1) = struct('senscutoff', 0.8, ...
    'pdatanum', 1, ...
    'rcvBufFrame',1, ...
    'IntBufDest', [1,1], ...
    'ImgBufDest', [0,0], ...
    'newFrameTimeout', UF.Time_loop*1e3,...  % needed for the software to wait for a new frame in the first rcv slot
    'RINums',1:UF.numFrames*2*na);


a = 3*na*UF.numFrames; % Reconstruction after Hadamard processing so in RcvBuffer #2
page = 1;

for frame = 1:UF.numFrames
    
    rcv_pulse = 1; % the index of the receive acquisition in buffer 2 (After hadamard processing)
    % Pulse 1 and 2 are already summed (in hardware) and sent in Rcv#1[2] in buffer2
    txNum = 1 ;
    for acqNum = 1:na %
        ReconInfo(2*na*(frame-1) + rcv_pulse) = struct('mode', 'replaceIQ', ...   % replace IQ data every time
            'txnum', txNum, ...
            'rcvnum', a + 2*na*(frame-1) + rcv_pulse, ...       % the index of the receive acquisition (here, we reconstruct the RF after hadamard processing)
            'pagenum', page,...
            'regionnum', 1);
        rcv_pulse = rcv_pulse + 1 ;
        txNum = txNum + 1 ;
        page  = page + 1;
    end
    
    
    txNum = 1 ;
    for acqNum = na+1:2*na
        ReconInfo(2*na*(frame-1) + rcv_pulse) = struct('mode', 'replaceIQ', ...   % replace IQ data every time
            'txnum', txNum, ...
            'rcvnum', a + 2*na*(frame-1) + rcv_pulse, ...       % the index of the receive acquisition
            'pagenum', page,...
            'regionnum', 1);
        rcv_pulse = rcv_pulse + 1;
        txNum = txNum + 1 ;
        page = page + 1;
    end
    
end

%_________________________________________________________________________%
%%                                                                        %

% Specify Process structure array.

% Save IQ data
Process(1).classname = 'External';     % process structure for 1st Doppler ensemble
Process(1).method = 'saveIQData_MPW';
Process(1).Parameters = {'srcbuffer','inter',... % name of buffer to process.
    'srcbufnum',1,...
    'srcframenum',1,...     % one inter buffer only
    'dstbuffer','none'};
EF(1).Function = text2cell('%#EF-1');


% Display block time
Process(2).classname = 'External';     % process structure for 1st Doppler ensemble
Process(2).method = 'dispBlockTime';
Process(2).Parameters = {'srcbuffer','none',...
    'dstbuffer','none'};
EF(2).Function = text2cell('%#EF-2');


% Hadamard summing
Process(3).classname = 'External';
Process(3).method = 'post_acq_Had_sum_AM';
Process(3).Parameters = {'srcbuffer','none',...  % name of buffer to process.
    'dstbuffer','receive',...
    'dstbufnum',2,...
    'dstframenum',1};
EF(3).Function = text2cell('%#EF-3');


% Specify SeqControl structure arrays.
maxRoundTrip = 3*sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2)*2; % Multiplane wave: longer transmission/reception time and back and forth: x2
maxTOF = round(2*maxRoundTrip/Trans.frequency);  %  x2 to be sure
time_compound = na*maxTOF;
time_ensemble = 1/(UF.dopFrameRate*1e-6)-(time_compound-maxTOF);
time_2_nextSEQ = UF.Time_loop*1e6- 1/(UF.dopFrameRate*1e-6)*UF.numFrames - time_ensemble;

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start, not used here
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = maxTOF;  % 200 usecs
SeqControl(3).command = 'triggerOut';
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'timeToNextAcq';  % time between frames
SeqControl(5).argument = 2*time_ensemble;  % 2 msec (<=> 500 Hz)
SeqControl(6).command = 'timeToNextAcq';  % time between blocks
SeqControl(6).argument = time_2_nextSEQ;  % 3 sec
SeqControl(7).command = 'sync';  % synchronisation soft hard, (no ping pong)
SeqControl(7).argument = UF.Time_loop*1e6;  % time out for sync  % needed for the software to wait for the hardware to end the last TTNA of the loop
SeqControl(8).command = 'loopCnt'; % - Set loop count. for looping on blocs
SeqControl(8).argument = UF.NbOfBlocs-1;  %
SeqControl(8).condition = 'counter1';
SeqControl(9).command = 'loopTst';  % loop test
SeqControl(9).argument = [];    % set apres
SeqControl(9).condition = 'counter1';

nsc = 10; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;

% set loop count
Event(n).info = 'start counter';
Event(n).tx = 0;   % use next TX structure.
Event(n).rcv = 0;
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = 8;
n = n+1;
SeqControl(9).argument = n;


Event(n).info = 'Acquire RF';
Event(n).tx = na + 1;   % use next TX structure.
Event(n).rcv = 1;
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = 2;
n = n+1;


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


Event(n).info = 'Call external function';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 3;
Event(n).seqControl = nsc;
SeqControl(nsc).command = 'markTransferProcessed';
SeqControl(nsc).argument = nsc-2; % reference previous transferToHost
nsc = nsc+1;
n = n+1;

Event(n).info = 'Reconstruct';
Event(n).tx = 0;        % no transmit
Event(n).rcv = 0;       % no rcv
Event(n).recon = 1;     % reconstruction
Event(n).process = 0;	% process
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Sequence Time Control';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 2;    % process
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'save IQ Data Process';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 1;    % process
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'loop everything';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % process
Event(n).seqControl = 9 ;
n = n+1;

% Rcv profile
RcvProfile.DCsubtract = 'on';               % substract DC signal if 'on'
RcvProfile.AntiAliasCutoff = UF.AntiAliasingFilter;    % antialiasing analogical filter cuttoff freq [MHz]
RcvProfile.PgaGain = 30;	% 24 ou 30	% analog gain in dBfor preamp  %#test
RcvProfile.LnaGain = 24;	% gain in dB of the fixed gain low noise amp
RcvProfile.LnaZinSel = 31;	% Force high-Z state for best Doppler sensitivity  %#test

%Check coherence of the VS structures
if max([Event(:).rcv])>length(Receive)
    warning('Probably a problem of receive indexing in the Event struct')
elseif max([Event(:).tx])>length(TX)
    warning('Probably a problem of TX indexing in the Event struct')
elseif max([Event(:).seqControl])>length(SeqControl)
    warning('Probably a problem of SeqControl indexing in the Event struct')
elseif max([ReconInfo(:).txnum])>length(TX)
    warning('Probably a problem of TX indexing in the ReconInfo struct')
    % elseif max([ReconInfo(:).pagenum])>UF.numFrames
    %     warning('Probably a problem of interbuffer pages management in the ReconInfo struct')
elseif max([ReconInfo(:).rcvnum])>length(Receive)
    warning('Probably a problem of receive indexing in the ReconInfo struct')
end

% UI controls

% Set TPCHighVoltage for profile one to 20V
UI(1).Statement = '[result,hv] = setTpcProfileHighVoltage(UF.ImgVoltage,1);';
UI(2).Statement = 'hv1Sldr = findobj(''Tag'',''hv1Sldr'');';
UI(3).Statement = 'set(hv1Sldr,''Value'',hv);';
UI(4).Statement = 'hv1Value = findobj(''Tag'',''hv1Value'');';
UI(5).Statement = 'set(hv1Value,''String'',num2str(hv,''%.1f''));';


% Save all the structures to a .mat file. & auto start
save(['MatFiles/',filename]);
disp([ mfilename ': NOTE -- Running VSX automatically!'])

VSX
return


%#EF-1
saveIQData_MPW(IData,QData)
IQData=complex(IData,QData);
persistent bloc_count;
global na UF dir_save start_bloc
if isempty(bloc_count)
    bloc_count = start_bloc;
end

Rcv_IQ = squeeze(IQData);
Rcv_IQ = reshape(Rcv_IQ, size(Rcv_IQ, 1), size(Rcv_IQ, 2), 2*na, []);

IQ_pulse1and2 = Rcv_IQ(:,:,2:na,:);
IQ_pulse3 = Rcv_IQ(:,:,2+na:2*na,:);


IQ_pulse1and2 = squeeze(mean(IQ_pulse1and2, 3));   % do the compounding (less disk space needed)
IQ_pulse3 = squeeze(mean(IQ_pulse3, 3));   % do the compounding (less disk space needed)

IQ_pulse1and2_re = real(IQ_pulse1and2);
IQ_pulse1and2_im = imag(IQ_pulse1and2);

IQ_pulse3_re = real(IQ_pulse3) ;
IQ_pulse3_im = imag(IQ_pulse3);

IQ_AM_neg = [IQ_pulse1and2_re IQ_pulse1and2_im];
IQ_AM_pos = [IQ_pulse3_re IQ_pulse3_im] ;

file_neg = sprintf('IQ_AMneg_%.3d.bin', bloc_count);
fid = fopen([dir_save file_neg],'w');   % fast save
fwrite(fid,IQ_AM_neg, 'double');
fclose(fid);

file_pos = sprintf('IQ_AMpos_%.3d.bin', bloc_count);
fid = fopen([dir_save file_pos],'w');   % fast save
fwrite(fid,IQ_AM_pos, 'double');
fclose(fid);
bloc_count = bloc_count+1;
return
%#EF-1



%#EF-2
dispBlockTime()
persistent time_bloc;
global dir_save;
if isempty(time_bloc)
    time_bloc = tic;
end
T = toc(time_bloc) % absolute time (ms)
D = datestr(now, 'HH-MM-SS.FFF'); % fUS image interval (s)
fidresult = fopen([dir_save '\timeLog.txt'], 'a');
fprintf(fidresult,'%s  %g\r\n', D, T);
fclose(fidresult);
%#EF-2


%#EF-3
RFOut = post_acq_Had_sum_AM()

persistent frame_num bloc_num HA
global dir_save na UF start_bloc
load([dir_save '/Receive.mat'])
if isempty(bloc_num)
    bloc_num = start_bloc;
end
if isempty(HA)
    HA = evalin('base','H');
end
HA = int16(HA);

fid = fopen([dir_save , sprintf('RData_%.3d.bin', bloc_num)], 'r');
RData = fread(fid, 'int16');
fclose(fid);

RData = reshape(RData,[],128);
RData = int16(RData);


clear RFOut
RFOut = zeros(size(RData,1),size(RData,2));  % make zeroed copy of RData
RFOut=int16(RFOut);


for f=1:UF.numFrames
    
    clear ss_1and2 es_1and2 ss_3 es_3
    % Determine startSample and endSample for the various acqNums.
    ss_1and2 = zeros(1,na);
    es_1and2 = zeros(1,na);
    
    ss_3 = zeros(1,na);
    es_3 = zeros(1,na);
    
    
    for j= 1:na
        
        % Pulses 1 and 2 : pulse 2 [ Receive(na)->Receive(2*na) ] is
        % summed in the hardware with pulse 1 [ Receive(1)->Receive(na)]
        % So we only take the Receive data after summing --> Receive(na)-> Receive(2*na)
        
        ss_1and2(j) = Receive(3*na*(f-1) + na+j).startSample;
        es_1and2(j) = Receive(3*na*(f-1) + na+j).endSample;
        
        % Pulse 3
        ss_3(j) = Receive(3*na*(f-1) + 2*na+j).startSample;
        es_3(j) = Receive(3*na*(f-1) + 2*na+j).endSample;
        
    end
    
    for i=1:na
        for j=1:na
            % Pulses 1 and 2
            RFOut(ss_1and2(i):es_1and2(i),:) = RFOut(ss_1and2(i):es_1and2(i),:) + HA(j,i)*RData(ss_1and2(j):es_1and2(j),:);
            % Pulse 3
            RFOut(ss_3(i):es_3(i),:) = RFOut(ss_3(i):es_3(i),:) + HA(j,i)*RData(ss_3(j):es_3(j),:);
        end
    end
end

bloc_num = bloc_num +1 ;

return
%#EF-3

