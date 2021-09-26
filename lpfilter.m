function [H, Htime] = lpfilter(type, M, N, D0, n)
%LPFILTER Computes frequency domain lowpass filters
%   H = LPFILTER(TYPE, M, N, D0, n) creates the transfer function of
%   a lowpass filter, H, of the specified TYPE and size (M-by-N).  To
%   view the filter as an image or mesh plot, it should be centered
%   using H = fftshift(H).
%
%   Valid values for TYPE, D0, and n are:
%
%   'ideal'    Ideal lowpass filter with cutoff frequency D0.  n need
%              not be supplied.  D0 must be positive
%
%   'btw'      Butterworth lowpass filter of order n, and cutoff D0.
%              The default value for n is 1.0.  D0 must be positive.
%
%   'gaussian' Gaussian lowpass filter with cutoff (standard deviation)
%              D0.  n need not be supplied.  D0 must be positive.

% Use function dftuv to set up the meshgrid arrays needed for 
% computing the required distances.
[U, V] = dftuv(M, N);

% Compute the distances D(U, V).
D = sqrt((U - size(U,1)/2 ).^2 + (V- size(V,2)/2 ).^2);

%% only for nothc filter
% y = 50;x=50;x1 = 50; y1 = 100;
% D = sqrt((U - size(U,1)/2 - y).^2 + (V- size(V,2)/2 - x).^2);
% D1 = sqrt((U - size(U,1)/2 + y).^2 + (V- size(V,2)/2 + x).^2);
%  
% D2 = sqrt((U - size(U,1)/2 - y).^2 + (V- size(V,2)/2 + x).^2);
% D3 = sqrt((U - size(U,1)/2 + y).^2 + (V- size(V,2)/2 - x).^2);
% 
% D4 = sqrt((U - size(U,1)/2 - y1).^2 + (V- size(V,2)/2 - x1).^2);
% D5 = sqrt((U - size(U,1)/2 + y1).^2 + (V- size(V,2)/2 + x1).^2);
% 
% D6 = sqrt((U - size(U,1)/2 - y1).^2 + (V- size(V,2)/2 + x1).^2);
% D7 = sqrt((U - size(U,1)/2 + y1).^2 + (V- size(V,2)/2 - x1).^2);
 
% Begin fiter computations.
switch type
case 'ideal'
  % H = 1 - (double(D <=D0)&double(D >=D0 -3));
  H = double(D <=D0);
%   H1 = double(D<=D0+10);
%   H2 = double(D<=D0-15);
%   H = H1-H+H2;
%   H2 = double(D2<=D0);
%   H3 = double(D3<=D0);
%   H = 1 - (H+H1+H2+H3);   % notch filter
%    H = ones(1024,1024);
%    H(1:512-25,508:516) = 0;H(512+25:1024,508:516) = 0;
case 'btw'
   if nargin == 4
      n = 1;
   end
   n = 1;
   H = 1./(1 + (D./D0).^(2*n));
   %H1 = 1./(1 + (D1./D0).^(2*n));
%    H =1-H.*H1;   % notch filter
case 'gaussian'
   H = exp(-(D.^2)./(2*(D0^2)));
   
   %% For Bandpass
%    H1 = exp(-(D.^2)./(2*((D0+15)^2)));
%    H = H1-H;

%% For notch filter
%   H1 = exp(-(D1.^2)./(2*(D0^2)));
%   H2 = exp(-(D2.^2)./(2*(D0^2)));
%   H3 = exp(-(D3.^2)./(2*(D0^2)));
%   
%   H4 = exp(-(D4.^2)./(2*(D0^2)));
%   H5 = exp(-(D5.^2)./(2*(D0^2)));
%   H6 = exp(-(D6.^2)./(2*(D0^2)));
%   H7 = exp(-(D7.^2)./(2*(D0^2)));
%  %   H = 1 - (H+H1+H2+H3);   % notch filter
%   H = 1-(H+H1+H2+H3+H4+H5+H6+H7);   % notch filter
case 'laplacian'
   H = -4*(pi^2)*(D.^2);
otherwise
   error('Unknown filter type.')
end

%%%%view the magnitude spectrum
S=log(1+abs((H)));
figure, imshow(S,[])
figure, mesh(abs((H)))
colormap('default') 

%%%%time domain
Htime = real(ifft2(ifftshift(H))); % time domain filter in [0:M-1,0:N-1]
realHtime = ((Htime - min(min(Htime)))./(max(max(Htime))) - min(min(Htime)));% time domain filter with center at origin
figure,imshow(fftshift(realHtime),[])
figure, mesh(fftshift(realHtime))
colormap('default')
figure

% this is for visualization only.
 