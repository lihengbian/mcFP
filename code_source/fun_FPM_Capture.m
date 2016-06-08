function [sample, f_sample, im_capture, fprob_real, fprob_save, dkx, dky, kx, ky, Masks, overlapratio, motionxy, NAfil] = fun_FPM_Capture(samplename, noise, n, M_factor, fprob_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, Jan. 5th, 2016
% Contact: lihengbian@gmail.com
% This function implements the simulation of FPM measurements in the case with sample motion.
% Thanks to Xiaoze Ou for offering code samples.

% Inputs:
% samplename.amplitude: amplitude image name
% samplename.phase: phase image name or 0
% noise.type: 'Guassian', 'poisson', 'speckle' or 'fluctuation'
% noise.variance is the variance of the noise
% n: pixels of the HR sample
% M_factor: magnification factor between HR and LR
% fprob_flag: if using real aberrated pupil function

% Outputs:
% sample, f_sample: groundtruth of HR sample and its spatial spectrum
% fprob_real: real aberrated pupil function
% fprob_save: ideal pupil function (all ones in the NA circle)
% dkx, dky, kx, ky: wave vector
% Masks: the upper left pixel location of each sub-spectrum
% overlapratio: overlap ratio between two adjacent sub-spectra
% motionxy: groundtruth of motion shift
% NAfil: pixel number of the radius in NA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeronum = 0;
modulus = imresize(im2double(imread(samplename.amplitude)), [n-2*zeronum, n-2*zeronum]);
if samplename.phase == 0
    phase = zeros([n, n]);
else
    phase = imresize(im2double(imread(samplename.phase)), [n-2*zeronum, n-2*zeronum]);
    phase = phase - min(min(phase));
    phase = phase/max(max(phase));
    phase = ((phase)*pi-pi/2)/2;
end
sample = modulus .* exp(1j*phase);
sample = padarray(sample,[zeronum zeronum]);
f_sample=fftshift(fft2(sample)); % low frequency in center

%% system parameters
if strcmp(noise.type, 'LEDp')
    LEDp=noise.variance;
else
    LEDp=4;%LED pitch distance, mm
end
H=84.8; %height from LED to sample, mm
wlength=0.625e-6; %wavelength of red light, m
NA=0.08; %numerical aperture of the objective
psize=0.2e-6; %pixel size in the reconstructed image, m

overlapratio = overlap(NA, H,LEDp)

snum = 7; % half LED numbers
arraysize = snum*2+1;%LED arraysize
L = arraysize^2;

%% calculate kx, ky for LED (center (0,0))
centerx=0;centery=0;
xlocation=zeros(1,(2*snum+1)^2);
ylocation=zeros(1,(2*snum+1)^2);
lightupseq=gseq(2*snum+1)-0.01;% problem for floor (6/3)=2, I want 1

for tt=1:(2*snum+1)^2
    xi=centerx-snum+round(mod(lightupseq(1,tt),2*snum+1))-1;
    yi=centery-snum+floor(lightupseq(1,tt)/(2*snum+1)); % problem for floor (6/3)=2, I want 1
    xlocation(1,tt)=xi-centerx;
    ylocation(1,tt)=yi-centery;
end;
clear tt xi yi lightupseq
dkx=2*pi/(psize*n);
dky=2*pi/(psize*n);
kx=2*pi./wlength*(xlocation*LEDp./sqrt(xlocation.^2.*LEDp.^2 + ylocation.^2.*LEDp.^2 + H.^2)); % kx = 2pi/lambda*sin(theta)   pixel shift = psize*n/lambda*sin(theta)
ky=2*pi./wlength*(ylocation*LEDp./sqrt(xlocation.^2.*LEDp.^2 + ylocation.^2.*LEDp.^2 + H.^2));

tempkx = kx;
tempky = ky;
tempkx(abs(kx)>3.152e6 | abs(ky)>3.152e6) = [];
tempky(abs(kx)>3.152e6 | abs(ky)>3.152e6) = [];
kx = tempkx;
ky = tempky;

L = length(kx)

%% capture image
NAfil=round(NA*(1/wlength)*n*psize); % pixel number of the radius in NA
mask=zeros(n,n);
[km, kn]=meshgrid(1:n,1:n);
mask((((km-n/2)/NAfil).^2+((kn-n/2)/NAfil).^2)<1)=1;
fprob_save=mask;

def=0.1;ax=-0.21;ay=0.13;cx=0.0;cy=0.0;

zn=def*gzn(n,NAfil*2,0,2)+ax*gzn(n,NAfil*2,2,2)...
    +ay*gzn(n,NAfil*2,-2,2)+cx*gzn(n,NAfil*2,1,3)...
    +cy*gzn(n,NAfil*2,-1,3);%gzn(n,NAfil*2,0,0)+

if fprob_flag == 1
    fprob_real=mask.*exp(pi*1j.*zn);
else
    fprob_real = mask;
end

hwidth = n/M_factor/2;

Masks = zeros([L,2]);
for i = 1:size(kx,2)
    Masks(i,1) = n/2+round(kx(1,i)/dkx)-hwidth+1;
    Masks(i,2) = n/2+round(ky(1,i)/dky)-hwidth+1;
end
fmaskpro = fprob_real(n/2-hwidth+1:n/2+hwidth,n/2-hwidth+1:n/2+hwidth);
if strcmp(noise.type,'samplemotion') == 1
    [im_capture, motionxy] = A_Function_Real_samplemotion(f_sample, Masks, n/M_factor, n/M_factor, fmaskpro, sample, noise.variance);
end
im_capture = abs(im_capture).^2;

end