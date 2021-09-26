% close all
figure
clear all
close all

%% read the image
 footBall=imread('x5.bmp');
% footBall=imcrop(rgb2gray(imread('my.jpg')),[1 1 511 511]);
% footBall=imcrop(rgb2gray(imread('newsprint.jpg')),[1 1 147 147]);
% footBall = footBall1(1:64,1:64);
% footBall=imcrop(rgb2gray(imread('Fig0451(a)(satellite_original).jpg')),[1 1 511 511]);

%% create an image
[M,N]=size(footBall);
% M=512;N=512;
% footBall(1:M,1:N/2)=zeros(M,N/2);
% footBall(1:M,N/2 + 1:N)=255*ones(M,N/2);
figure
imshow(footBall)
title('original image')
fim = double(footBall);
% fim = (fim - min(min(fim)))./(max(max(fim))- min(min(fim)));
% footBall = fim; % for Laplacian
% % Convert to grayscale
% F=fftshift(fft2(double(footBall)));
% S2 = log(1 + abs(F));
% figure,imshow(S2,[])

%% horizontal line image
v=zeros(511,511);
v(256,:)=1;
%v(250,:)=1;
%v(240,:)=1;
figure, imshow(v,[])
Fd=(fft2(double(v)));
Sd = log(1 + abs((Fd)));
figure,imshow(Sd,[])


%% display FT of some simple signals
v=ones(1,511);
c=diag(v);
figure, imshow(c,[])
Fd=(fft2(double(c)));
Sd = log(1 + abs(Fd));
figure,imshow(Sd,[])

%% translation invariance for magnitude response
v=zeros(511,511);
v(500,:)=1;
figure, imshow(v,[])
Fd=fftshift(fft2(double(v)));
Sd = log(1 + abs(Fd));
figure,imshow(Sd,[])



close all
%Determine good padding for Fourier transform
% PQ = paddedsize(size(footBall));

%% generate the filter in frequency domain, first define the meshgrid
PQ = [3*M 3*N];
u = 0:(2*M-1);
v = 0:(2*N-1);
[U,V] = meshgrid(v, u);

%% Create a Lowpass filter 5% the width of the Fourier transform
D0 = 0.05*PQ(1);
[H, Htime ] = lpfilter('gaussian', PQ(1), PQ(2), D0); %%%% non cetered filter

%% Calculate the discrete Fourier transform of the image
F=fftshift(fft2(double(footBall),size(H,1),size(H,2))); %%%FT of the image for same size as that of filter
imshow(mat2gray(log(1+abs(F))))
imtool(mat2gray(log(1+abs(F))))
title('magnitude spectrum of original image')
%% using manual padding and centering via multiplication with (-1).^(U+V)
% padfootBall = padarray(footBall, [512 512], 0, 'post');
% F=(fft2(double(padfootBall).*((-1).^(U+V)),size(H,1),size(H,2)));
% imshow(mat2gray(log(1+abs(F))))
% title('magnitude spectrum of original image using manual padding')


%% Apply the highpass filter to the Fourier spectrum of the image ie Hadamard product
LPFS_football = (H).*F; %%%%%%%%apply filter in Fourier domain 


%% convert the result to the spatial domain. 

LPF_football=real(ifft2(ifftshift((LPFS_football)))); %%%%%filtered image
% LPF_football=real(ifft2(((LPFS_football)))).*((-1).^(U+V)); 
%  LPF_football=real(ifft2(((LPFS_football)))); 


%% Crop the image to undo padding
LPF_football=LPF_football(1:size(footBall,1), 1:size(footBall,2));  %%%%%%%%%get only relevant portion of the image
%Display the blurred image
figure, imshow(LPF_football, [])  % display filtered image
title('filtered image')
% foot = padarray(footBall,[64 64],0,'post');
% % realH = padarray(realHtime,[63 63],0,'post');
% T = ifft2(fft2(double(foot).*fft2(realHtime)));
% rT = real(T);
% foot = padarray(foot,[127 127],0,'post');
% realH = padarray(realHtime,[127 127],0,'post');
% Z = (conv2(double(footBall),(realHtime)));figure,imshow(Z,[])
% Display the Fourier Spectrum 
% Move the origin of the transform to the center of the frequency rectangle.


%% display FT of original and filtered image

Fcf=(LPFS_football);
% use abs to compute the magnitude and use log to brighten display
S=log(1+abs(Fcf)); %magnitude FFT of filtered image
figure, imshow(S,[])
title('magnitude spectrum of filtered image')

figure
imshow(footBall)
title('original image')

