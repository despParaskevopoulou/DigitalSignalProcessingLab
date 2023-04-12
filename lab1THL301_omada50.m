%Askhsh 1
%part A

%Creating the first signal
clear all;
clear;
nx = -15:1:15;
x = 3*sinc(pi/4.*nx);
figure(1);
subplot(1,3,1);
stem(nx,x);
axis([-17 17 -2 3]);
grid on;
title('x signal','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

%Creating the second signal
ny = -15:1:15;
y = 2*cos(pi/6.*ny+pi/4);
subplot(1,3,2);
stem(ny,y);
axis([-17 17 -2 3]);
grid on;
title('y signal','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

%reversing the second signal
n = [nx(1)+ny(1):nx(end)+ny(end)];
y_rev = y(end:-1:1);
subplot(1,3,3);
stem(ny,y_rev);
axis([-17 17 -2 3]);
grid on;
title('signal y reversed','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

len = length(n);

%zero padding for the first signal
x0 = [zeros(1,length(y)-1) x zeros(1, length(y)-1)];

%convolution
for i=1:len
    y_rev0 = [zeros(1,(i-1)) y_rev zeros(1,(len-i))];
  h(i) = sum(x0.*y_rev0);
end;
figure(2);
subplot(1,2,1);
stem(n,h);
axis([-32 32 -8 8]);
grid on;
title('without conv function','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

%convolution with 'conv' function
h=conv(x,y);
subplot(1,2,2);
stem(n,h,'m');
axis([-32 32 -8 8]);
grid on;
title('with conv function','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

%part B
%using the previous signals x,y
X = fft(x,len);
Y = fft(y,len);
h_f = ifft(X.*Y);
figure(3);
subplot(1,2,1);
stem(n,h);
grid on;
axis([-32 32 -8 8]);
title('Convolution in Time','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

subplot(1,2,2);
stem(n,h_f,'m');
grid on;
axis([-32 32 -8 8]);
title('Multiplication in Frequency','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

%Askhsh 2
clear all;
clear;
dt = 0.0005;
t = 0:dt:0.5;
x = 5*cos(24*pi.*t) - 2*sin(1.5*pi.*t);
figure(4);
plot(t,x);
title('signal x','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');
grid on;

%sampling with fs = 48Hz
fs = 48;
ts = 1/fs;
t1 = 0:ts:0.5;
x1 = 5*cos(24*pi.*t1) - 2*sin(1.5*pi.*t1);
figure(5);
subplot(1,3,1);
plot(t,x);
hold on
plot(t1,x1,'r.');
hold off;
title('Sampling points(red) fs = 48','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

%sampling with fs = 24Hz
fs2 = 24;
ts2 = 1/fs2;
t2 = 0:ts2:0.5;
x2 = 5*cos(24*pi.*t2) - 2*sin(1.5*pi.*t2);
subplot(1,3,2);
plot(t,x);
hold on;
grid on;
plot(t2,x2,'r.');
hold off;
title('Sampling points(red) fs = 24','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

%sampling with fs = 12Hz
fs3 = 12;
ts3 = 1/fs3;
t3 = 0:ts3:0.5;
x3 = 5*cos(24*pi.*t3) - 2*sin(1.5*pi.*t3);
subplot(1,3,3);
plot(t,x);
hold on;
grid on;
plot(t3,x3,'r.');
hold off;
title('Sampling points(red) fs = 12','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

%sampling with fs = 50Hz
figure(6);
fsA = 50;
tsA = 1/fsA;
tA = 0:tsA:0.5;
xA = 5*cos(24*pi.*tA) - 2*sin(1.5*pi.*tA);
plot(t,x);
hold on;
grid on;
plot(tA,xA,'r.');
hold off;
title('Sampling points(red) with fs = 50','fontweight','bold');
xlabel('time','fontweight','bold');
ylabel('amplitude','fontweight','bold');

%Askhsh 3
%Part A
clc;
clear all;

figure(7);
dt = 0.0001;

%sampling frequency
fs = 100;
ts = 1/fs;

t_final = 128*ts - ts;
t = 0:dt:t_final;
x = (10*cos(2*pi*20.*t) - 4*sin((2*pi*40.*t)+5));
plot(t,x);
hold on;
grid on;

t2 = 0:ts:t_final;
x2 = (10*cos(2*pi*20.*t2) - 4*sin((2*pi*40.*t2)+5));
plot (t2, x2, 'r.');
title('Blue: signal x - - - - - Red: sampling (128 samples)','fontweight','bold');
xlabel('Time','fontweight','bold');
ylabel('Amplitude','fontweight','bold');
hold off;

F = [-fs/2:fs/128:fs/2-fs/128];
X = fftshift(abs(fft(x2)));
figure(8);
plot(F,X);
title('Frequency Spectrum of signal x','fontweight','bold');
xlabel('Frequency (Hz)','fontweight','bold');
ylabel('Amplitude','fontweight','bold');
grid on;

n = 0:0.01:0.5;
xn = sin(2*pi*(100/8)*n + 50);
figure(9)
stem(n,xn);
grid on;
xlabel('Time axis','Fontsize',12);
ylabel('Amplitude','Fontsize',12);
title('Signal xn','Fontsize',14);


for f0=100:125:475
    fs=8000; %sampling frequency
    Ts=1/fs;
    N=300;     
    n=0:N-1;
    x_3=sin(2*pi*f0*n*Ts+50);
  
    f_axis=-fs/2:fs/N:fs/2-1; %frequency axis
    Xf_3=fftshift(abs(fft(x_3)));  %Fourier transform
    figure; %10,12,14,16
    stem(f_axis,Xf_3);
    grid on;
    xlabel('Frequency(Hz)','Fontsize',12);
    ylabel('Values of function','Fontsize',12);
    title('The Fourier transform of x(t)=cos(2*pi*f0*t+50)','FontSize',14);
end



for f0=7525:125:7900
    fs=8000; %sampling frequency
    Ts=1/fs;
    N=300;     
    n=0:N-1;
    x_3=sin(2*pi*f0*n*Ts+50);
   
    f_axis=-fs/2:fs/N:fs/2-1; %frequency axis
    Xf_3=fftshift(abs(fft(x_3)));  %Fourier transform
    figure; %18,20,22,24
    stem(f_axis,Xf_3);
    grid on;
    xlabel('Frequency(Hz)','Fontsize',12);
    ylabel('Values of function','Fontsize',12);
    title('The Fourier transform of x(t)=cos(2*pi*f0*t+50)','FontSize',14);
end