numbands = 32;
numpoly = numbands;    % polyphase component number = decimation ratio = number of channels
M = numpoly;                    % polyphase system decimation ratio

h(1,:) = fir1(2*M-1, 1/M, "low");
h(M,:) = (-1).^(0:2*M-1).*h(1,:);
for i = 1:M-2
    h(i+1,:) = fir1(2*M-1, [i/M (i+1)/M], "bandpass");
end
L = length(h);
         
hhap = zeros(numbands,numpoly, L/numpoly);  % hhap = 32 x 32 x 512/32 , 3D ANALYSIS filter bank array

for k=1:numbands
    for m = 1:numpoly
        hhap(k,m,:) = h(k,m:M:end);       % create the m-th polyphase component of k-th channel filter
    end
end



hhsp = conj(permute(hhap, [2 1 3]));


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
    a = [wav_in(1:M:end),0];
    temp = conv(a, b);
    for m=2:M
        temp = temp + conv([0,wav_in(M-m+2:M:end)],squeeze(hhap(k,m,:)));   
    end
    yyd(k,:) = temp(L/(2*M)+1 : L/(2*M)+N/numbands);
end

% S7 - Apply SYNTHESIS filterbank on the decimated signal :
% ---------------------------------------------------------
ys = zeros(1, N);

for k=1:numbands
    temp = zeros(1, N+L-1);
    for m = 1:numpoly
        temp(m:numbands:end-31) = conv( yyd(k,:) , squeeze(hhsp(k,m,:)) );
    end

    ys = ys + temp(L/2+1:L/2+N);    
end
ys = numbands*ys;


% SX - DISPLAY RESULTS:
% ---------------------
L = length(h);
figure,subplot(2,1,1)
stem(0:L-1,h(1,:));title('The Prototype Lowpass Filter');
subplot(2,1,2)
plot(linspace(-1,1,N),20*log10(abs(fftshift(fft(h(1,:),N)))));
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
