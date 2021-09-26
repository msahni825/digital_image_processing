%% perform spatial convolution using FFT
close all
clear all
clc
a=[1 1 0 1 ];b=[2 0 1 0];  % time domain
c=[1 1 0 1 0 0 0]; d = [2 0 1 0 0 0 0];  % zero padded
fftres = real(ifft(fft(c).*fft(d)))  % convolution via FFT - correct
fftrescenter = real(ifft(ifftshift((fftshift(fft(c)).*fftshift(fft(d))))))  % convolution via FFT - correct

fftresincorrect = real(ifft(fft(a).*fft(b))) %convoltion via FFT - incorrect

convres = conv(a,b) % spatial conv, should be equal to convolution via FFT - correct

subplot(4,1,1)
stem(a)
axis([-1 8 0 3])
title('a')
subplot(4,1,2)
stem(b)
axis([-1 8 0 3])
title('b')
subplot(4,1,3)
stem(convres)
axis([-1 8 0 3])
title('spatial convolution')
subplot(4,1,4)
stem(fftresincorrect)
axis([-1 8 0 3])
title('convolution obtained via FFT without padding')

figure

subplot(5,1,1)
stem(c),axis([-1 8 0 3])
title('c')
subplot(5,1,2)
stem(d),axis([-1 8 0 3])
title('d')
subplot(5,1,3)
stem(convres),axis([-1 8 0 3])
title('spatial convolution')
subplot(5,1,4)
stem(fftres),axis([-1 8 0 3])
title('convolution obtained via FFT with padding')
subplot(5,1,5)
stem(fftrescenter),axis([-1 8 0 3])
title('convolution obtained via centered FFT with padding')

%% 2d conv via FFT
 clc    
 a=[0 -2; 2 3]; 
 b = [1 2; 3 2];
 
 c = [0 -2 0 ;2 3 0 ;0 0 0 ];
 d = [1 2 0;3 2 0; 0 0 0];
 
 fftres = ifft2(fft2(c).*fft2(d))% convolution via FFT - correct
 
 fftresincorrect = ifft2(fft2(a).*fft2(b))%convoltion via FFT - incorrect
 
 convres = conv2(a,b) % spatial conv, should be equal to convolution via FFT - correct
 

%% Note the time
a = randi(512, 512);
b = randi(512, 512);
c=padarray(a,[511 511],0,'post');
d=padarray(b,[511 511],0,'post');
tic
fftres = (ifft2(fft2(c).*fft2(d)));
fftresimage = real(fftres(1:512,1:512));
toc
% fftresincorrect = ifft2(fft2(a).*fft2(b));
tic
convres = conv2(a,b);
toc

diff = mean(mean(abs(fftres) - abs(convres)))



%% another example

a=[1 2; 3 4]; 
b = [1 2 3; 4 5 6;7 8 9];
c = [5 6 0 0 4;8 9 0 0 7;0 0 0 0 0;2 3 0 0 1];

d = [1 2 0 0 0;3 4 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
fftres = ifft2(fft2(c).*fft2(d))% convolution via FFT - correct

fftresincorrect = ifft2(fft2(a,3,3).*fft2(b))%convoltion via FFT - incorrect

 
