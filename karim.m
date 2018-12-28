clc
clear all
close all

 addpath /home/masters/Desktop/toolbox/toolbox/toolbox_signal
 addpath /home/masters/Desktop/toolbox/toolbox/toolbox_general

load('greasy.mat')
s0=s';

fs = 16000;
x=1/sqrt(2);
t=0:1/fs:(length(s)-1)/fs;
f=0:fs/length(s):fs-(fs/length(s));
% figure;
% plot(t,s);
% title('orginal signal vs time');

s=abs(fft(s));
% 
% figure;
% plot(f,s(1:length(s/2)));
% title('original signal representation in freq domain');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=0:0.01:0.5;
z=exp(2*pi*f*1i);

Ho=(1/sqrt(2))*(1+z.^-1);

% figure;
% plot(f,abs(Ho));
% line([0.25 0.25],[1 0],'color','k');
% line([0.25 0],[1 1],'color','k');
% title('haar filter');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=0:0.01:0.5;
z=exp(2*pi*f*1i);

Go=(1/sqrt(2))*(1-z.^-1);
% 
% figure;
% plot(f,abs(Go));
% line([0.25 0.25],[1 0],'color','k');
% line([0.25 0],[1 1],'color','k');
% title('haar filter  complementary');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hlopass=[-1/sqrt(2) 1/sqrt(2)];

s1low=pconv(hlopass,s0);

g0 = -((-1).^(1:length(hlopass))).* hlopass;

%%%

hoHipass=[1/sqrt(2) 1/sqrt(2)];

s2high=pconv(hoHipass,s0);

g0 = -((-1).^(1:length(hoHipass))).* hoHipass;

% figure;
% plot(t,s1low);
% title('low pass applied');
% 
% figure;
% plot(t,s2high);
% title('high pass applied');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract vector
n2=16;
sto =   s(1:n2:1024); %extracted vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   h0=[1 1]./sqrt(2);
   res=filterbank(s0,3,h0,0);
   s2bank=filterbank(res,3,h0,1);

% figure;
% plot(t,res);
% title('decomposition signal')
% 
% figure;
% plot(t,s2bank);
% title('bank filter reconstructed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%2nd part
name = 'lena';
n = 256;
f0 = load_image(name);
f0 = rescale(crop(f0,n));

% figure;
% subplot(3,1,1)
% 
% imageplot(f0);
% title('original image');

%%%%%%%%%%%%%%%%%%%%%%%5
imageNEw=fwt_or_2d(0,f0,3,h0);

% 
% subplot(3,1,2)
% imageplot(imageNEw);
% title('decomposed image');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imageRecon=fwt_or_2d(1,imageNEw,3,h0);

% 
% subplot(3,1,3)
% imageplot(imageRecon);
% title('Reconstructed image');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transfer function for DDaubechies’
HD=(326/675)+((1095/1309)*z.^-1)+((648/2891)*z.^-2)+((675/5216)*z.^-3);
HardThresh = @(x,t)x.*(abs(x)>t);

vn=0.1;


f0noise=f0+vn*randn(size(f0)); %noised the image
ns=256;
jmin=3;
optimi=1;
% figure;
% imshow(f0noise,[]);

w=perform_wavelet_transf(f0noise,jmin,1,optimi);
% figure;
% imageplot(w,[]);
wt=perform_thresholding(w,vn*2,'soft');
resf=perform_wavelet_transf(wt,jmin,-1,optimi);
% figure;
% imageplot(resf);

%%
%3.2
rho = .3;
Lambda = rand(n,n)>rho;
Phi = @(f)f.*Lambda;
y = Phi(f0);
figure;
imageplot(y); %damaged image
%%%%%%%%%%%%%%%%%%

%T=linspace(5,0.1,30);

SoftThresh = @(x,T)x.*max(0, 1-T./max(abs(x),1e-10));
temp=y;
Jmax = log2(n)-1;
Jmin = Jmax-3;
options.ti = 0; % use orthogonality.
Psi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);
SoftThreshPsi = @(f,T)Psi(SoftThresh(PsiS(f),T));


for i=1:5 %50

    R=SoftThreshPsi(temp,.1);
    temp = y+ not(Lambda).*R;



end



% 
% figure;
% imageplot([ y clamp(temp)]);
% %imageplot(clamp(SoftThreshPsi(f0,.1)) );
% figure;
 T = linspace(5,0.5,30);
% plot(T, SoftThresh(T,0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


K = @(f)grad(f);
KS = @(u)-div(u);
Amplitude = @(u)sqrt(sum(u.^2,3));
F = @(u)sum(sum(Amplitude(u)));
ProxF = @(u,lambda)max(0,1-lambda./repmat(Amplitude(u), [1 1 2])).*u;
ProxFS = @(y,sigma)y-sigma*ProxF(y/sigma,1/sigma);
ProxG = @(f,tau)f + Phi(y - Phi(f));
L = 8;
sigma = 10;
tau = .9/(L*sigma);
theta = 1;
f = y;
g = K(y)*0;
f1 = f;

figure;

for i=1:300
        fold = f;
        g = ProxFS( g+sigma*K(f1), sigma);
        f = ProxG(f-tau*KS(g), tau);
        f1 = f + theta * (f-fold);
        drawnow
        %plot(F(K(f*i)));
        imageplot(f1);
        
end


























% Res=fwt_or_2d(0,f0noise,3,h0);
% 
% lph = Res([1:ns/(2.^3), 1:ns/(2.^3)]);
% 
% 
% wfirstScale=[ lph(ns/2:ns, 1:ns) , lph(1:ns/2),((ns/2+1):ns) ];
% 
% MAD=median(abs(wfirstScale(:)));
% sigma=MAD/0.6745;
% tnew=3*sigma;



 %imageDaubechies=fwt_or_2d(0,f0noise,3,lph);
 
 
% 
% figure;
% subplot(3,1,2)
% imshow(imageDaubechies);
% title('decomposed image with Daubechies’ transfer fn ');






