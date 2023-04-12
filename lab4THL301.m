% 4h ergasthriakh askhsh

clc;
clear all;

% Askhsh 1
wc = 0.4*pi;
fs = 100;
N = 21; %length
fc = wc/(2*pi);
wn = fc/(fs/2); %wn = 0.004

% Creating the two windows
RecWindow = rectwin(N);
HamWindow = hamming(N);

% Creating the two filters
filter1 = fir1(N-1, wn, RecWindow);
filter2 = fir1(N-1, wn, HamWindow);

[h1, w1] = freqz(filter1, 512);
[h2, w2] = freqz(filter2, 512);

figure(1);

plot(w1,abs(h1),w2,abs(h2));

legend('Hamming','Rectangular');
xlabel('Frequency','FontSize',10);
ylabel('Magnitude','FontSize',10);
title('Frequency response of filters using Hamming & Rectangular window','FontSize',12,'FontWeight','bold');

grid on;

%askisi 2
%2a
clear all;
clc;

Wc=0.5*pi;
Fs=100;
Fc=Wc/(2*pi);
N1=21;
N2=41;

%normalising the frequency
wn=Fc/(Fs/2);

%creating the filter
hammfil1=fir1(N1-1,wn,hamming(N1));

hammfil2=fir1(N2-1,wn,hamming(N2));

hannfil1=fir1(N1-1,wn,hanning(N1));

hannfil2=fir1(N2-1,wn,hanning(N2));


%getting the frequency response
[hHamm1,wHamm1]=freqz(hammfil1,N1);

[hHamm2,wHamm2]=freqz(hammfil2,N2);

[hHann1,wHann1]=freqz(hannfil1,N1);

[hHann2,wHann2]=freqz(hannfil2,N2);

figure (2);

subplot(1,2,1);
plot(wHamm1,abs(hHamm1));

xlabel('Frequency','FontSize',10);
ylabel('Magnitude', 'FontSize',10);
title('Hamming: N=21','FontSize',12,'FontWeight','bold');
grid on;

subplot(1,2,2);
plot(wHamm2,abs(hHamm2));
xlabel('Frequency','FontSize',10);
ylabel('Magnitude', 'FontSize',10);
title('Hamming: N=41','FontSize',12,'FontWeight','bold');
grid on;

figure(3);

subplot(1,2,1);
plot(wHann1,abs(hHann1));
xlabel('Frequency','FontSize',10);
ylabel('Magnitude', 'FontSize',10);
title('Hanning: N=21','FontSize',12,'FontWeight','bold');
grid on;

subplot(1,2,2);
plot(wHann2,abs(hHann2));
xlabel('Frequency','FontSize',10);
ylabel('Magnitude', 'FontSize',10);
title('Hanning: N=41','FontSize',12,'FontWeight','bold');
grid on;

%2b
Fs=100;
Ts=1/Fs;
w1=15;
w2=200;
f1=15/(2*pi);
f2=200/(2*pi);

N=512;
n=0:N-1;
x=sin(2*pi*n*f1*Ts)+0.25*sin(2*pi*n*f2*Ts);

f_axis=-Fs/2:Fs/N:Fs/2-Fs/N;
X = fftshift(abs(fft(x)));

fil1=filter(hammfil1,1,x);
X1=fftshift(fft(fil1));
fil2=filter(hammfil2,1,x);
X2=fftshift(fft(fil2));
fil3=filter(hannfil1,1,x);
X3=fftshift(fft(fil3));
fil4=filter(hannfil2,1,x);
X4=fftshift(fft(fil4));


figure(4);
subplot(3,1,1)
plot(f_axis,abs(X));
xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;

subplot(3,1,2);
plot(f_axis,abs(X1));

xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of filtered (Hamming) x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;

subplot(3,1,3);
plot(f_axis,abs(X2));

xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of filtered (Hamming) x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;

figure(5)

subplot(3,1,1)
plot(f_axis,abs(X));
xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;

subplot(3,1,2);
plot(f_axis,abs(X3));
xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of filtered (Hanning) x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;

subplot(3,1,3);
plot(f_axis,abs(X4));
xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of filtered (Hanning) x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;

%2c
%we follow the procedure form above but with fs = 50Hz
Fs=50;
Ts=1/Fs;
w1=15;
w2=200;
f1=15/(2*pi);
f2=200/(2*pi);

N=512;
n=0:N-1;
x=sin(2*pi*n*f1*Ts)+0.25*sin(2*pi*n*f2*Ts);

f_axis=-Fs/2:Fs/N:Fs/2-Fs/N;
X = fftshift(abs(fft(x)));

fil1=filter(hammfil1,1,x);
X1=fftshift(fft(fil1));
fil2=filter(hammfil2,1,x);
X2=fftshift(fft(fil2));
fil3=filter(hannfil1,1,x);
X3=fftshift(fft(fil3));
fil4=filter(hannfil2,1,x);
X4=fftshift(fft(fil4));

figure(6);
subplot(3,1,1)
plot(f_axis,abs(X));

xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;

subplot(3,1,2);
plot(f_axis,abs(X1));

xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of filtered (Hamming) x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;

subplot(3,1,3);
plot(f_axis,abs(X2));

xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of filtered (Hamming) x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;

figure(7)
subplot(3,1,1)
plot(f_axis,abs(X));

xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;

subplot(3,1,2);
plot(f_axis,abs(X3));

xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of filtered (Hanning) x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;


subplot(3,1,3);
plot(f_axis,abs(X4));

xlabel('Frequency(Hz)','Fontsize',10);
ylabel('Values of function', 'Fontsize',10);
title('Spectrum of filtered (Hanning) x=sin(15*t)+0.25*sin(200*t)','Fontsize',12,'FontWeight','bold');
grid on;