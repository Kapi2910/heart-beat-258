% S0 - Load the prototype lowpass filter impulse response h0[n]:
% --------------------------------------------------------------
M= 32;
h = fir1(2*M-1, 1/M, "low");         % h[n] is the prototype lowpass filter of length 512
L = length(h);

% S1 - Create the 32 x 512 analysis filter-bank hha[k,n] by cosine modulation from protoype :
% -----------------------------------------------------------------------------------------
numbands = 32;                 % number of banks (channels)
n=0:L-1;

hha=zeros(numbands,L);  % bank of filters hha[k,n] = 32 x 512 array.

hha(1,:) = h;
hha(numbands,:) = (-1).^(0:2*numbands-1).*h;
for i = 1:numbands-2
    hha(i+1,:) = fir1(2*M-1, [i/M (i+1)/M], "bandpass");
end


% S2 - Create the 32-polyphase components hhap[k,m,n] , for each one of 32 analysis filters hha[k,n]:
% ---------------------------------------------------------------------------------------------------
numpoly = numbands;             % polyphase component number = decimation ratio = number of channels
hhap = zeros(numbands,numpoly, L/numpoly);  % hhap = 32 x 32 x 512/32 , 3D ANALYSIS filter bank array

M = numpoly;                    % polyphase system decimation ratio
for k=1:numbands
    for m = 1:numpoly
        hhap(k,m,:) = hha(k,m:M:end);       % create the m-th polyphase component of k-th channel filter
    end
end


% S3 - Design the 32 x 512  synthesis (complementary) filter bank :
% -----------------------------------------------------------------
numbands = 32;                  % number of banks
n=0:L-1;
hhs = zeros(numbands,L);        % bank of filters
hhs(1,:) = h;
hhs(numbands,:) = (-1).^(0:2*numbands-1).*h;
for i = 1:numbands-2
    hhs(i+1,:) = fir1(2*M-1, [i/M (i+1)/M], "bandpass");
end


% S4 - Obtain the 32-polyphase components hhsp[k,m,n] , for each one of 32 synthesis filters hhs[k,n]:
% ----------------------------------------------------------------------------------------------------
numpoly = numbands;             % polyphase component number = interpolation ratio = number of channels
hhsp = zeros(numbands,numpoly, L/numpoly);  % hhap = 32 x 32 x 512/32 , 3D ANALYSIS filter bank array
M = numpoly;                    % polyphase system decimation ratio
for k=1:numbands
    for m = 1:numpoly
        hhsp(k,m,:) = hhs(k,m:M:end);       % create the m-th polyphase component of k-th channel filter
    end
end


% S5 - Generate the test input signal
% -----------------------------------
N = 8*L;
wav_in = heart(1:N);        % pure sine tone

% S6 - Apply test signal to the filterbank: ANALYSIS STAGE :
% -----------------------------------------------------------
yyd = zeros( numbands, floor(N/numbands));   % decimated outputs..
M = numbands;
for k=1:1:numbands  
    b = squeeze(hhap(k,1,:));
    a = wav_in(1:M:end);
    temp = filtfilt(b,1,a);
    for m=2:M
        temp = temp + filtfilt(squeeze(hhap(k,m,:)),1,wav_in(M-m+2:M:end));   
    end
    yyd(k,:) = temp(L/(2*M) : L/(2*M)+N/numbands-1);
end

% S7 - Apply SYNTHESIS filterbank on the decimated signal :
% ---------------------------------------------------------
ys = zeros(1, N);
temp = zeros(1, N);
for k=1:numbands
    for m = 1:numpoly
        temp(m:numbands:end) = filtfilt(squeeze(hhsp(k,m,:))',1,yyd(k,:) );
    end

    ys = ys + temp;    
end
ys = numbands*ys;


% SX - DISPLAY RESULTS:
% ---------------------
L = length(h);
figure,subplot(2,1,1)
stem(0:L-1,h);title('The Prototype Lowpass Filter');
subplot(2,1,2)
plot(linspace(-1,1,N),20*log10(abs(fftshift(fft(h,N)))));
grid on;

figure
plot(linspace(-1,1,N),20*log10(abs(fftshift(fft(hha(1,:),N)))));
hold on
for k=2:numbands
    plot(linspace(-1,1,N),20*log10(abs(fftshift(fft(hha(k,:),N)))));
end
title('32 CHANNEL FILTERBANK');

figure,subplot(2,1,1)
plot(wav_in);title('input signal')
subplot(2,1,2)
plot(linspace(-1,1,4*N),20*log10(abs(fftshift(fft(wav_in,4*N)))));

figure,subplot(2,1,1)
plot(ys);title('Synthesized Back');
subplot(2,1,2)
plot(linspace(-1,1,4*N),20*log10(abs(fftshift(fft(ys,4*N)))));