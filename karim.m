clc
clear all
close all

load('greasy.mat')
s0=s';

fs = 16000;
x=1/sqrt(2);
t=0:1/fs:(length(s)-1)/fs;
f=0:fs/length(s):fs-(fs/length(s));
figure;
plot(t,s);
title('orginal signal vs time');

s=abs(fft(s));

figure;
plot(f,s(1:length(s/2)));
title('original signal representation in freq domain');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=0:0.01:0.5;
z=exp(2*pi*f*1i);

Ho=(1/sqrt(2))*(1+z.^-1);

figure;
plot(f,abs(Ho));
line([0.25 0.25],[1 0],'color','k');
line([0.25 0],[1 1],'color','k');
title('haar filter');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=0:0.01:0.5;
z=exp(2*pi*f*1i);

Go=(1/sqrt(2))*(1-z.^-1);

figure;
plot(f,abs(Go));
line([0.25 0.25],[1 0],'color','k');
line([0.25 0],[1 1],'color','k');
title('haar filter  complementary');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hlopass=[-1/sqrt(2) 1/sqrt(2)];

s1low=pconv(hlopass,s0);

g0 = -((-1).^(1:length(hlopass))).* hlopass;

%%%

hoHipass=[1/sqrt(2) 1/sqrt(2)];

s2high=pconv(hoHipass,s0);

g0 = -((-1).^(1:length(hoHipass))).* hoHipass;

figure;
plot(t,s1low);
title('low pass applied');

figure;
plot(t,s2high);
title('high pass applied');



