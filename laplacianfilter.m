clear all
close all
clc
f=imread('Fig0338(a)(blurry_moon).tif');
% f=double(f)./max(max(double(f)));

%% zeropad image
M=size(f,1); N=size(f,2); % nr of rows/columns of image f
C=3; D=3; % nr of rows/columns of kernel h
P=M+C-1; Q=N+D-1; % nr of rows/columns after padding
fp=zeros(P,Q); % zero padding: start with zeroes
fp(1:M,1:N)=f; % insert f into image fp

%% define filter
hp=zeros(P,Q); % Construct filter matrix hp, same size as fp.
hp(1,1)=-8; hp(2,1)=1; hp(1,2)=1; % Center is at (1,1)
hp(P,1)=1; hp(1,Q)=1; % Indices modulo P or Q
hp(P,2) = 1; hp(2,Q) = 1;hp(2,2)=1;hp(P,Q) = 1;

%% perform DFT on both image and filter
Fp=fft2(double(fp), P, Q); % FFT of image fp
Hp=fft2(double(hp./8), P, Q); % FFT of kernel hp

%% derived from spatial filter with horizontal and vertical. diagonals are zeros
for k = 0:P-1
    for l = 0:Q-1
        Hp(k+1,l+1) = 2*cos(2*pi*k/P) + 2*cos(2*pi*l/Q) -4;
    end
end
%Hp = 1- Hp./(min(min(Hp)));

%% visualize filter magnitude response
H = fftshift(Hp); % centered FFT of kernel
F1 = abs(H); % Get the magnitude
F1 = log(F1+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
F1 = mat2gray(F1,[0 1]);
imshow(F1);%visualize kernel in Fourier domain
title('centered magnitude spectrum')
figure

%% Hadamard product
Gp=Fp .* Hp;

%% Inverse DFT, real and top left crop size of image
gp=ifft2(Gp); % Inverse FFT, division by 8 is to ensure that values are in [-1,1]
gp=real(gp); % Take real part
g=gp(1:M, 1:N);
imshow(g,[])  %%% Laplacian image
title('Laplacian image')
 
%% sharpening f - Laplaian image
fim =double(f);
gsharp =mat2gray(fim  - g)*255;
figure,imshow(uint8(gsharp))   %%%Sharp image
title('sharpened image')


%% show original image
figure
imshow(f,[])   %%  Original image
title('original image')
