clc;
clear all;
close all;

%askhsh 1
fs=10000;
Wp=3000*2*pi;
Ws=4000*2*pi;
delta_p=3;
delta_s=30;

[n,Wc]= buttord(Wp,Ws,delta_p,delta_s,'s');

[z,p,k]=buttap(n);

disp("The order of the analog order of the signal is :" );
disp(n);
[num,den]=zp2tf(z,p,k);

[num1,den1]=lp2lp(num,den,Wc);

f= linspace(0,fs/2,2048);
H_analog=freqs(num1,den1,2*pi*f);
H_analog_dB=20*log(abs(H_analog));


[Num1,Den1]=bilinear(num1,den1,fs);
H_digital=freqz(Num1,Den1,f,fs);
H_digital_dB= 20*log(abs(H_digital));

figure(1)
plot(f,H_analog_dB,'--',f,H_digital_dB);
grid on;
legend('Analog filter', 'Digital filter');
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Magnitude (dB)','FontSize',14);
title('Butterworth Lowpass Filter','FontSize',15,'FontWeight','bold');


%gia delta_s=50dB
delta_s1=50;

[n,Wc]= buttord(Wp,Ws,delta_p,delta_s1,'s');

[z,p,k]=buttap(n);

disp("The order of the analog order of the signal is :" );
disp(n);
[num11,den11]=zp2tf(z,p,k);

[num12,den12]=lp2lp(num11,den11,Wc);

f= linspace(0,fs/2,2048);
H_analog1=freqs(num12,den12,2*pi*f);
H_analog_dB1=20*log(abs(H_analog1));


[Num12,Den12]=bilinear(num12,den12,fs);
H_digital1=freqz(Num12,Den12,f,fs);
H_digital_dB1= 20*log(abs(H_digital1));

figure(2)
plot(f,H_analog_dB1,'--',f,H_digital_dB1);
grid on;
legend('Analog filter', 'Digital filter');
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Magnitude (dB)','FontSize',14);
title('Butterworth Lowpass Filter','FontSize',15,'FontWeight','bold');

%ASKHSH 2
figure(3);

N1 = 2;
ripple = 3;
wc = 2;
ts = 0.2;
fs1 = 1/ts;
fc = wc/(2*pi);
fk = fc/(fs1/2);

[Num, Den] = cheby1(N1, ripple, fk, 'high');
H1 = freqz(Num, Den, 256);
H1_dB = 20*log(abs(H1));

N2 = 16;
[Num2, Den2] = cheby1(N2, ripple, fk, 'high');
H2 = freqz(Num2, Den2, 256);
H2_dB = 20*log(abs(H2));

f = linspace(0, 1, 256);

plot(f,H1_dB,f,H2_dB);
axis([-0.1 1.1 -500 20]);
grid on;

legend('N=2', 'N=16');
xlabel('W (rad/sec)','FontSize',14);
ylabel('Magnitude (dB)','FontSize',14);
title('Chebyshev Highpass Filter','FontSize',15,'FontWeight','bold');

%askisi 3
%first part
f1 = 500/pi;
f2 = 8000/pi;
f3 = 15000/pi;
fs = 10000;
Ts = 1/fs;
N=500;                  % samples   
n=0:N-1;
t_final = 500*Ts - Ts;
t = 0:Ts:t_final;
x1 = 1+cos(2*pi*f1*n*Ts) + cos(2*pi*f2*n*Ts)+cos(2*pi*f3*n*Ts);

figure(4)
stem(n,x1);
grid on;
xlabel('Frequency(Hz)','Fontsize',14);
ylabel('Values of function', 'Fontsize',14);
title('The first 500 samples of x(t)=1+cos(1000t)+cos(16000t)+cos(30000t)','Fontsize',12,'FontWeight','bold');

%fourier transform of the filter and the spectrum
faxis = [-fs/2:fs/N:(fs/2-1)];
X = fftshift(abs(fft(x1)));
figure(5);
subplot(3,1,1);
stem(faxis,X);
xlabel('Frequency (Hz)','fontweight','bold');
ylabel('Amplitude','fontweight','bold');
title('The Spectrum','FontSize',12)
grid on;

%filter 1
xfil = filter(Num1,Den1,x1);

X1 = fftshift(abs(fft(xfil)));

subplot(3,1,2);
stem(faxis,X1);
xlabel('Frequency (Hz)','fontweight','bold');
ylabel('Amplitude','fontweight','bold');
title('The Spectrum of the filtered signal, attenuation 30dB','FontSize',12)
grid on;


%filter 2
%figure(6)
xfil2 = filter(Num12,Den12,x1);

X2 = fftshift(abs(fft(xfil2)));
subplot(3,1,3);
stem(faxis,X2);
xlabel('Frequency (Hz)','fontweight','bold');
ylabel('Amplitude','fontweight','bold');
title('The Spectrum of the filtered signal, attenuation 50dB','FontSize',12)
grid on;


%3b
F1 = 0.75/pi;
F2  = 2.5/pi;
ts = 0.2;
fs1=1/ts;
t_final1 = 500*ts - ts;
t1 = 0:ts:t_final1;
x2 = 1+cos(2*pi*F1*n*ts) + cos(2*pi*F2*n*ts);

figure(6)

stem(n,x2);
grid on;
xlabel('Frequency(Hz)','Fontsize',12);
ylabel('Values of function', 'Fontsize',12);
title('The first 500 samples of x(t)=1+cos(1.5t)+cos(5t)','Fontsize',12,'FontWeight','bold');

%fourier transform of the filter and the spectrum
faxis1 = [-fs1/2:fs1/N:(fs1/2)-(fs1/N)];
X2 = fftshift(abs(fft(x2)));
figure(7)
subplot(2,1,1)
stem(faxis1,X2);
xlabel('Frequency (Hz)','fontweight','bold');
ylabel('Amplitude','fontweight','bold');
title('The Spectrum','FontSize',12)
grid on;

x2fil = filter(Num2,Den2,x2);

X2f = fftshift(abs(fft(x2fil)));
subplot(2,1,2)
stem(faxis1,X2f);
xlabel('Frequency (Hz)','fontweight','bold');
ylabel('Amplitude','fontweight','bold');
title('The Spectrum of the filtered signal','FontSize',12)
grid on;

