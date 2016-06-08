function [MM, NN, immax]=fastreg_Liheng(standimage,compimage)

% Liheng, May 7th 2016
% add immax in outputs

% function [m n]=fastreg(standimage,compimage):
% Inputs
% standimage: the first image
% compimage:  the second image. It should be the same size as the first image

% Outputs
% m: the shift in X
% n: the shift in Y


[M N]=size(standimage);
standimage=double(standimage);
compimage=double(compimage);
R0=2;
R1=R0+5;
% pixel level registration 
[M0, N0, temp, immax]=regsurf(standimage,compimage);

M0=floor(M0-M/2-1);
N0=floor(N0-N/2-1);

% % MM=M0+1;
% % NN=N0+1;
MM=M0;
NN=N0;

function [m, n, im, immax]=regsurf(standimage,compimage)
s=fft2(standimage);
c=ifft2(compimage);
sc=s.*c;
im=abs(fftshift(ifft2(sc)));
immax = max(im(:));
[M0, N0]=find(im==immax);
m=round(mean(M0));
n=round(mean(N0));

