%Lab 2 THL 301
%Askhsh 1
%b
clear all;
clc;
num=[0 0.2 0];
den=[1 -0.7 -0.18];
ts=0.1;
sys= tf(num,den,ts)
zeros = roots(num)
poles = roots(den)
zplane(zeros,poles);


%d
b = [0 0.2 0];
a = [1  -0.7 -0.18];
sys2 = tf(b,a,ts)
f=-pi:pi/128:pi;
figure(2);
freqz(b,a,f);
figure(3);
freqz(b,a);

%e
a1 = [1 -1.7 0.52 0.18];
sys1 = tf(b,a1,ts)
figure(4)
freqz(b,a1,f);
%title('lrrlim larlom');

%Askisi 2
%a
num1 = [4 -3.5 0];
den1 = [1 -2.5 1];
%sys1 = tf(num1,den1)

[r,p,k] = residuez(num1,den1)

syms z;
H1 = r(1)/(1-p(1)*z^(-1));
H2 = r(2)/(1-p(2)*z^(-1));
H=H1+H2;
pretty(H)

%b
hn=iztrans(H)
