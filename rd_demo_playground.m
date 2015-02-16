% This script is intented as a simple introductory to compressed sensing
% via random demodulation using a multi-tone signal as a show case.

clear all;
close all;

%% Presentation options
%usePause = false; % Determines whether the script pauses between plotting of signals
usePause = false; % Determines whether the script pauses between plotting of signals
usePlot = false; 
plotSize = 400;


%% signal parameters
fmin = 220;                    % minimum frequency
fMax = 350;                    % maximum frequency of a single tone
mTones = 5;

% discrete representation parameters
fNyquist = 2*fMax;            % Nyquist freqyency
overSample = 6;              % overSampling factor (in order to represent the signal with higher definition)
fs = overSample*fNyquist;     % Sampling frequency
T = 1/fs;                     % Sample time
ncycle = 1;                   % number of periods of the signal with minumum frequency
t = (0:T:ncycle-T);           % Time vector

% subsampling settings
subSampling = 3;            % undersampling ratio
eSNRvec = [30, 40, 50, 60, 70];
inSNRvec = zeros(100, 5);
recSNR_cvxVec = zeros(100, 5);
recSNR_spgl1Vec = zeros(100, 5);
recSNR_nestaVec = zeros(100, 5);
recSNR_aihtVec = zeros(100, 5); 
cvx_time_vec = zeros(100, 5);
spgl1_time_vec = zeros(100, 5);
nesta_time_vec = zeros(100, 5);
aiht_time_vec = zeros(100, 5);


for jj=1:5
eSNR = eSNRvec(jj);                        % desired SNR of the signal

for iii=1:100

%% generating input signal
% Sum of  randomly generated sine waves
singleToneM = zeros(mTones,fs*ncycle);        % pre-allocate matrix of single-tone signals
fsig = zeros(mTones);
for ii=1:mTones
    fsig(ii) = randi([fmin, fMax],1);
    singleToneM(ii,:) =  4*abs(rand(1)).*sin(2*pi*fsig(ii)*t);  %generate sine-signal
end
signalIdeal = sum(singleToneM);% add singleTone signals to create multi-tone

% pollute signal with Additive White-Gaussian Noise (AWGN)
signal = awgn(signalIdeal, eSNR, 'measured'); % signal with AWGN


%% generate sparsity basis #Psi
%creating dense DFT matix    - Fourier basis
N = length(signal);               % length of the input signal
Psi = dftmtx(N);
% inverse of the #Psi
invPsi = conj(Psi)/N;

if usePlot
    % Plot - generated signal in time and sparse domain
    figure(1)                           % ## FIGURE 1 ##
    sSize =  400;
    subplot(2,1,1);                             % #1
    plot(fs*t(1:sSize),signal(1:sSize));
    title('input signal in time domain')


    sSize =  400;
    subplot(2,1,2);                             % #2
    fsignal = abs(signal*Psi);
    plot(fs*t(1:sSize),fsignal(1:sSize));
    title('input signal spectrum')
end

%PAUSE
if usePause
    pause
end

%% generate chipping sequence
chippingSequence = randsrc(1,N);  % ** |for the simplicity|
% ** In principle this should be generated according to desired
% ** frequency and then using [kron] adjusted to the size
% ** of the input signal - to perfrom modulation

if usePlot
    % Plot - pseudo-random sequence in time and sparse domain
    figure(2)                           % ## FIGURE 1 ##
    sSize =  400;
    %time
    subplot(2,1,1);                             % #1
    plot(fs*t(1:sSize),chippingSequence(1:sSize));
    title('pseudo-random sequence in time domain')
    % frquency
    sSize =  400;
    subplot(2,1,2);                             % #2
    fchippingSequence = abs(chippingSequence*Psi);
    plot(fs*t(1:sSize),fchippingSequence(1:sSize));
    title('pseudo-random  sequence in frequency')
end

%PAUSE
if usePause
    pause
end

%% modulate input signal with the chipping sequence
mixedSignal = signal.*chippingSequence;

if usePlot
    % Plot - generated signal in time and sparse domain
    figure(3)                           % ## FIGURE 1 ##
    sSize =  400;
    subplot(2,1,1);                             % #1
    plot(fs*t(1:sSize),mixedSignal(1:sSize));
    title('mudulated signal in time domain')

    %frequency
    sSize =  400;
    subplot(2,1,2);                             % #2
    fmixedSignal = abs(mixedSignal*Psi);
    plot(fs*t(1:sSize),fmixedSignal(1:sSize));
    title('modulated signal spectrum (spreaded spectrum)')
end

%PAUSE
if usePause
    pause
end

%% low-pass filter
%filtering settings
order = 7;             % filter order
fNorm = (0.5*(fMax/subSampling))/(fs/overSample);            % normalized cuf-off frequency

% matlab filter design tool: Butterworth filter design
[nom den] = butter(order, fNorm, 'low');    % order - filter order ...
... fNorm - normalized frequency
    
%impulse response of the filter
[ImpulseResponse ~] = impz(nom,den);
ResponseLength = length(ImpulseResponse);

% filtering
filterSignal = filter(nom, den, mixedSignal);

if usePlot
    % Plot - filtered signal in time and frequency domain
    figure(4)                           % ## FIGURE 4 ##
    sSize =  400;
    subplot(2,1,1);                             % #1
    plot(fs*t(1:sSize),filterSignal(1:sSize));
    title('filtered signal in time domain')

    %frequency
    sSize =  400;
    subplot(2,1,2);                             % #2
    ffilterSignal = abs(filterSignal*Psi);
    plot(fs*t(1:sSize),ffilterSignal(1:sSize));
    title('filtered signal spectrum')
end 

%PAUSE
if usePause
    pause
end

%% measure the signal (low-rate sampling)
Ts = overSample*subSampling;

% sampling
Y = downsample(filterSignal, Ts, Ts-1); % downsample filtered signal
M = length(Y);                          % number of measurments taken

if usePlot
    % Plot - low-rate sampling
    figure(4)
    sSize = 400; 
    subplot(2,1,1);                             % #1
    hold on
    stem(Ts-1:Ts:sSize,Y(1,1:floor(sSize/Ts)),'green')
end

%PAUSE
if usePause
    pause
end

if usePlot
    % Show the plot of measured signal on the original input signal
    figure(1)
    subplot(2,1,1);                             % #1
    hold on
    stem(Ts-1:Ts:sSize,signal(1,Ts-1:Ts:sSize),'red')
end

%PAUSE
if usePause
    pause
end

%% generate measurement matrix #Phi
ImpulseResponse =  ImpulseResponse(ResponseLength:-1:1);
Phi = zeros(M,N);             % pre-allocate measurement matrix

% for j = 1:M
%     if j*Ts-Ts+ResponseLength < N
%         Phi(j,j*Ts-Ts+1:j*Ts-Ts+ResponseLength) = chippingSequence(1,j*Ts-Ts+1:j*Ts-Ts+ResponseLength).*ImpulseResponse(1:ResponseLength);
%     else
%         Phi(j,j*Ts-Ts+1:end) = chippingSequence(1,j*Ts-Ts+1:end).*ImpulseResponse(1:length(chippingSequence(1,j*Ts-Ts+1:end)));
%     end
% end

% generating Phi matrix - so that the values that can't fit are truncated
% from the left upper corner *( the most recent acquisition starts from
% bottom right corner of the matrix )
chippingSequence = chippingSequence';

for j = 1:M
    if j*Ts < ResponseLength
        Phi(j,1:j*Ts) = chippingSequence(1:j*Ts,1).*ImpulseResponse(ResponseLength-j*Ts+1:ResponseLength,1);
    elseif j*Ts == ResponseLength
        Phi(j,1:ResponseLength) = chippingSequence(1:ResponseLength,1).*ImpulseResponse(1:ResponseLength,1);
    else
        Phi(j,j*Ts-ResponseLength+1:j*Ts) = chippingSequence(j*Ts - ResponseLength+1:j*Ts,1).*ImpulseResponse(1:ResponseLength,1);
    end
end

if usePlot
    figure(5)
    subplot(1,2,1)                              % #1
    imagesc(Phi);
    title('Phi matrix')

    subplot(1,2,2)                            % #2
    imshow(Y');
    title('observed/compressed vector Y')
end

%PAUSE
if usePause
    pause
end

A = Phi*invPsi;


%% Reconstruction%
tic
%% L1 minimization by CVX            %
% ==================================%

cvx_begin quiet
variable alphaa(N) complex;    % define sparse vector alpha
minimize(norm(alphaa,1));
subject to
Phi*invPsi*alphaa == Y';
cvx_end
%------------------------------------%
cvx_time = toc();

tic
%% L1 minimization by SPGL1            %
% ==================================%
options = spgSetParms('verbosity', 0);
alphaa_spgl1 = spg_bpdn( A, Y', 0.001, options);
%------------------------------------%
spgl1_time = toc();



% alphaa_aiht = AIHT(Y, A', length(Y), mTones);
tic;
alphaa_aiht = AIHT(transpose(Y), A, size(A,2), mTones);
aiht_time = toc();

%% L1 minimization by NESTA          %
% ==================================%
% tic
% alphaa_nesta = NESTA(A, A.', Y', 3e-3, 0.001);
% nesta_time = toc();
% 
% alphaa_aiht = AIHT(Y,A',length(Y),mTones);
% 
% start = make_start_x('CGIHT_Matrix',size(A, 1),size(A, 2),mTones,real(A),Y);
% alphaa_cgiht = CGIHT_Matrix(size(A, 1),size(A, 2),mTones,real(A),Y,start);


%% time domain reconstructed signal
recSignal = real(invPsi*alphaa)';
recSignal_spgl1 = real(invPsi*alphaa_spgl1)';
recSignal_aiht = real(invPsi*alphaa_aiht)';

%% reconstruction evaluation ( time domain )
% CVX
recError = signal - recSignal;
maxError = max(recError);
inSNR = db(var(signalIdeal)/var(signal-signalIdeal), 'power');
recSNR = db(var(signal)/var(recError), 'power');

%SPGL1
recError_spgl1 = signal - recSignal_spgl1;
maxError_spgl1 = max(recError_spgl1);
recSNR_spgl1 = db(var(signal)/var(recError_spgl1), 'power');

% %NESTA
% recError_nesta = signal - recSignal_nesta;
% maxError_nesta = max(recError_nesta);
% recSNR_nesta = db(var(signal)/var(recError_nesta), 'power');
% 

%AIHT
recSNR_aiht = signal - recSignal_aiht;
maxError_aiht = max(recSNR_aiht);
recSNR_aiht = db(var(signal)/var(recSNR_aiht), 'power');

inSNRvec(iii, jj) = inSNR;
recSNR_cvxVec(iii,jj) = recSNR;
recSNR_spgl1Vec(iii,jj) = recSNR_spgl1;
% recSNR_nestaVec(iii,jj) = recSNR_nesta;
recSNR_aihtVec(iii, jj) = recSNR_aiht;

aiht_time_vec(iii, jj) = aiht_time;
cvx_time_vec(iii,jj) = cvx_time;
spgl1_time_vec(iii,jj) = spgl1_time; 
% nesta_time_vec(iii,jj) = nesta_time;

if ~mod(ii, 5)
save('RD_cvx_spgl1_nesta_test_i100_j5_mbpr_test_05_02_2015.mat', 'inSNRvec', 'recSNR_cvxVec', 'recSNR_spgl1Vec', 'recSNR_aihtVec' ,'cvx_time_vec', 'spgl1_time_vec', 'aiht_time_vec' );
end

end 
end

% 
% %% print
% fprintf('\n===============RECONSTRUCTION====================\n')
% fprintf('\n Nyquist frequency = %d Hz \n', max(fsig(1:3,1))*2)
% fprintf('\n Sampling frequency = %d Hz \n', M)
% fprintf('\n Undersamping factor ~ %3.2f times below Nyquist \n' , (max(fsig(1:3,1))*2)/M)
% fprintf('\n Signal input SNR was:    %5.3f dB \n' , inSNR)
% fprintf(' Signal output  SNR with CVX is:   %5.3f dB \n' , recSNR)
% fprintf(' Quality drop with CVX is:         %5.3f dB \n' , inSNR - recSNR)
% fprintf(' Signal output  SNR with SPGL1 is:   %5.3f dB \n' , recSNR_spgl1)
% fprintf(' Quality drop with SPGL1 is:         %5.3f dB \n' , inSNR - recSNR_spgl1)
% 
% 
% if usePlot
%     %% plot
%     % comparisment of the input signal and its reconstruction
%     figure(6)
%     sSize = 400;
%     p = plot(fs*t(1:sSize),signal(1:sSize));
%     set(p,'Color','blue','LineWidth',1)
%     hold on
%     %PAUSE
%     if usePause
%         pause
%     end
%     p2 = plot(fs*t(1:sSize),recSignal(1:sSize));
%     set(p2,'Color','green','LineWidth',1)
%     hold on
%     p3 = plot(fs*t(1:sSize), recError(1:sSize));
%     set(p3,'Color','red','LineWidth',1)
%     title('input signal, reconstructed signal & reconstruction error')
% end
