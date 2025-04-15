load("matfiles/100.mat")
fs = double(fs);
addpath("PR_Nguyen/")
%%
% figure;
% plot(heart)
% grid on
%% Baseline Wandering

% This section compares the detrending to filtering for the removal of
% baseline wander.
x = detrend(heart, 1);
h2 = highpass(heart,0.67, fs);

t = (0:length(heart)-1)/fs;

% f_base=pi/8;
% f1=0.5*f_base;
% f2=2.5*f_base;
% f3=4.5*f_base;
% f4=6.5*f_base;
% 
% x=sin(2*pi*f1*t) +sin(2*pi*f2*t) +sin(2*pi*f3*t) +sin(2*pi*f4*t);
% 


figure;
tiledlayout('vertical')
nexttile
plot(t, heart, 'DisplayName',"Original Version")
hold on
grid on
plot(t, x, 'DisplayName',"Detrended Version")
grid on
hold off
xlim([4 8.333])
ylim([-1.5 1.5])
xlabel('Time (s)')
ylabel('mV')
legend()
figure
plot(t, heart, 'DisplayName',"Original Version")
hold on
grid on
plot(t, h2, 'DisplayName',"Highpass Version")
grid on
hold off
xlabel('Time (s)')
ylabel('mV')
xlim([4 8.333])
ylim([-1.5 1.5])
legend()

%% 

% block the signal into blocks of size for a 32-channel MLT is 64
M = 32;
L = 2*M;
N=ceil(length(x)/M);
%%
clear H
m = 2;
H(1,:) = fir1(m*M-1, 1/M, "low");
H(M,:) = (-1).^(0:m*M-1).*H(1,:);
for i = 1:M-2
    H(i+1,:) = fir1(m*M-1, [i/M (i+1)/M], "bandpass");
end

F = H;
%%
S2 = zeros(M, N);
y = zeros(M, length(x));
for i = 1:M
    temp1 = filtfilt(F(i,:), 1, x);
    temp2 = downsample(temp1, M);
    S2(i,:) = temp2/M;
    temp3 = upsample(temp2,M);
    temp4 = filtfilt(H(i,:),1,temp3);
    y(i,:) =temp4(1:length(x)); 
end
y = sum(y)*M;

p1 = sum(abs(S2(2:4,:)),1);
p2 = sum(abs(S2(2:5,:)),1);
p3 = sum(abs(S2(3:5,:)),1);
p4 = sum((S2(2:4,:).^2),1);
p5 = sum((S2(2:5,:).^2),1);
p6 = sum((S2(3:5,:).^2),1);
m1 = movsum(p1,2);
e1 = zeros(size(m1));
e1(m1>0.02) = m1(m1>0.02);
figure;
plot(t, x, 'DisplayName','Original Signal');hold on;
plot(t, y, 'DisplayName','Reconstructed Signal'); hold off;
xlim([0 3])
ylim([-1.2 1.5])
xlabel('Time (s)')
ylabel('mV')
grid on;
legend();
%%
figure
tiledlayout('vertical')
nexttile;
plot(t, heart);
xlim([0 3])
title('Original Signal')
ylabel('mV')
xticks([])
grid on;
nexttile;
stem(t(1:M:end), 32000*S2(2,:),'Marker','none','LineWidth',1.33);
ylabel('A.U.')
title('[5.625, 11.25]')
xlim([0 3])
xticks([])
grid on;
nexttile;
stem(t(1:M:end), 32000*S2(3,:),'Marker','none','LineWidth',1.33);
ylabel('A.U.')
xlim([0 3])
title('[11.25, 16.875]')
xticks([])
grid on;
nexttile;
stem(t(1:M:end), 1000*S2(4,:),'Marker','none','LineWidth',1.33);
xlabel('Time (s)')
ylabel('A.U.')
title('[16.875, 22.5]')
xlim([0 3])
grid on;

%%
figure;
tiledlayout('vertical')
nexttile;
plot(t, x, 'DisplayName','Signal');
ylim([-1.2 1.5])
xlim([0 3])
title('Original Signal')
ylabel('mV')
xticks([])
nexttile;
stem(t(1:M:end), 32000*p1,'Marker','none','LineWidth',1.33); 
ylim([0 750])
ylabel('A.U.')
title('Feature P_1')
xlim([0 3])
xticks([])
nexttile;
stem(t(1:M:end), 32000*m1,'Marker','none','LineWidth',1.33);
ylim([0 1000])
ylabel('A.U.')
title('MWI')
xlim([0 3])
xticks([])
nexttile;
stem(t(1:M:end), 32000*e1,'Marker','none','LineWidth',1.33);
ylim([0 1000])
ylabel('A.U.')
title('Beat Detection')
xlim([0 3])
xlabel('Time (s)')



%% DFT Filter Bank
W = exp(-1i*2*pi*(0:M-1)'*(0:M-1)/M);
H = W';
F = W;

S1 = zeros(M, N);
y = zeros(1, length(x));
yy = zeros(size(S1));
xx = buffer(x,32);
for i = 1:N
    temp1 = H*xx(:,i);
    % temp2 = downsample(temp1, M);
    % temp3 = amp_coeff(i)*temp2;
    % temp4 = upsample(temp3, M);
    % temp5 = filter(F(i,:), 1, temp4);
    S1(:,i) = temp1;
    yy(:,i) = F*S1(:,i)/M;
end
yy = reshape(real(yy),1,[]);
y = yy(1:length(x));

p1 = sum(abs(S1(2:4,:)),1);
m1 = movsum(p1,2);
e1 = m1> 1.3;
tiledlayout('vertical')
nexttile;
plot(t, x, 'DisplayName','Signal');
ylim([-1.2 1.2])
legend();
nexttile;
plot(t, y, 'DisplayName','Signal');
ylim([-1.2 1.2])
legend();
nexttile;
plot(t(1:32:end), p1,'DisplayName','Feature P_1'); 
ylim([0 3])
legend();
nexttile;
plot(t(1:32:end), m1,'DisplayName','MWI');
ylim([0 3])
legend();
nexttile;
plot(t(1:32:end), e1,'DisplayName','Event');
ylim([0 1.2])
legend();


% X=fftshift(fft(x, 512) ) ;
% Y=fftshift(fft(y, 512)) ;
% plot(x)
% % plot (linspace(-pi, pi, length(X) ) , abs(X) , 'linewidth', 1) ;
% hold on
% plot(y)
% % plot (linspace(-pi, pi, length(Y) ) , abs(Y) , 'r', 'linewidth', 2) ;
% title("Uniform DFT")
% legend (' input', 'output' )
% % xlim ( [0 pi])

%%
clear H
H(1,:) = fir1(M-1, 1/M, "low");
H(M,:) = (-1).^(0:M-1).*H(1,:);
for i = 1:M-2
    H(i+1,:) = fir1(M-1, [i/M (i+1)/M], "bandpass");
end

F = inv(H);
F = F/max(F,[],1);
%%
T = zeros(size(H));
for i=1:M
    T(i,:) = H(i,:).*F(i,:);
end
T = sum(T,1);

%%
S2 = zeros(M, length(x));
y = zeros(M, length(x));
for i = 1:M
    temp1 = filtfilt(H(i,:), 1, x);
    temp2 = downsample(temp1, M);
    % S2(i,:) = temp1;
    temp3 = upsample(temp2,M);
    temp4 = filtfilt(F(i,:),1,temp3)/M;
    y(i,:) =temp4(1:length(x)); 
end
y = sum(y)/M;

p1 = sum(abs(S2(2:4,:)),1);
p2 = sum(abs(S2(2:5,:)),1);
p3 = sum(abs(S2(3:5,:)),1);
p4 = sum((S2(2:4,:).^2),1);
p5 = sum((S2(2:5,:).^2),1);
p6 = sum((S2(3:5,:).^2),1);
m1 = movsum(p1,2);
e1 = m1>1.5;
tiledlayout('vertical')
nexttile;
plot(t, x, 'DisplayName','Signal');hold on;
plot(t, y, 'DisplayName','Recon'); hold off;
ylim([-1.2 1.2])
legend();
nexttile;
plot(t, p1,'DisplayName','Feature P_1'); 
ylim([0 3])
legend();
nexttile;
plot(t, m1,'DisplayName','MWI');
ylim([0 3])
legend();
nexttile;
plot(t, e1,'DisplayName','Event');
ylim([0 1.2])
legend();



%%
clear H
H(1,:) = fir1(2*M-1, 1/M, "low");
for i = 2:M
    H(i,:) = real((exp(-2i*pi*i/M).*(fft(H(1,:)))));
end


%%
S2 = zeros(M, length(x));
for i = 1:M
    temp1 = filtfilt(H(i,:), 1, heart);
    temp2 = downsample(temp1, M);
    S2(i,:) = temp1;
end


p1 = sum(abs(S2(2:4,:)),1);
m1 = movsum(p1,2);
e1 = m1>1.5;
tiledlayout('vertical')
nexttile;
plot(t, heart, 'DisplayName','Signal');
ylim([-1.2 1.2])
legend();
nexttile;
plot(t, p1,'DisplayName','Feature P_1'); 
ylim([0 3])
legend();
nexttile;
plot(t, m1,'DisplayName','MWI');
ylim([0 3])
legend();
nexttile;
plot(t, e1,'DisplayName','Event');
ylim([0 1.2])
legend();


%%
figure;
for i=1
[H1(i,:),W]=freqz(H(i,:),1);
hold on
plot(W/pi,20*log10(abs(H1(i,:))),'k','linewidth',1.2)
axis([0 1 -120 30])
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
end

