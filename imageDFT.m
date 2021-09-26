%  f=imread('motionanvgFFT.jpg');
f=imread('x5.bmp');
[M,N] = size(f);
fzeropad = zeros(M+100,N+100);
fzeropad(30:M+30-1,10:N+10-1)=f;
f=fzeropad;
close all

imshow(f,[])
title('orig')
F = fft2(f);  % FFT
S = abs(F); % magnitude of FFT

phi = atan2(imag(F), real(F)); % phase of FFT taninv(real/imag)
F1 = S.*exp(1i*phi); %% reconstruction of image using it Fourier transform

origf = real(ifft2(F1));
figure,imshow(origf,[])
title('reconstructed')

%% Magnitude spectrum
f=imread('x5.bmp');

close all

imshow(f,[])

%%% find fft and display magnitude spectrum
F = fft2(f);
S = log(1 + abs(F));
S = mat2gray(S); 
figure,imshow(S,[])
title('magnitude spectrum')


phi = atan2(imag(F), real(F));
Phase = mat2gray(phi); 
figure,imshow(Phase,[])
title('Phase spectrum')

%%% find fft centered and display magnitude spectrum
Fc = fftshift(F);
Sc = log(1 + abs(Fc));
Sc = mat2gray(Sc); 
figure,imshow(Sc,[])
title('centered magnitude spectrum')

phic = atan2(imag(Fc), real(Fc));
Phasec = mat2gray(phic); 
figure,imshow(Phasec,[])
title('Centered Phase spectrum')


%% reconstruct using magnitude only, phase only and DC 
phi = atan2(imag(F), real(F));

F1 = abs(F);            % use magnitude only
origfMag = real(ifft2(F1));
figure,imshow(log(1+origfMag),[])

F2 = exp(1i*phi);       % use phase only
origfPh = real(ifft2(F2));
figure,imshow(log(1+origfPh),[])

FDc = zeros(size(F1,1),size(F1,2));     % used DC only
FDc(1,1) = F1(1,1);
F3 = FDc.*exp(1i*phi);
origfDc = real(ifft2(F3));
figure,imshow(uint8(origfDc))


%%  reconstruct using the centered Fourier transform
F = ifftshift(Fc);
phi = atan2(imag(F), real(F));
F1 = abs(F).*exp(1i*phi);

origf = real(ifft2(F1));
figure,imshow(origf,[])
title('reconstruct using centered FFT')
% reconstruct using magnitude only
F = fft2(f);
phi = atan2(imag(F), real(F));
F1 = abs(F);

origf = real(ifft2(F1));
origf = 255*(origf - min(min(origf)))./(max(max(origf)) - min(min(origf)));
% origf = mat2gray(origf);
figure,imshow((origf),[])
title('reconstruct using magnitude spectrum')
% reconstruct using phase only
F = fft2(f);
phi = atan2(imag(F), real(F));
F1 = exp(1i*phi);

origf = real(ifft2(F1));
figure,imshow(origf,[])
title('reconstruct using phase spectrum')