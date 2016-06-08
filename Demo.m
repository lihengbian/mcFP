%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, May 18th, 2016. Contact me: lihengbian@gmail.com.
% This demo does the simulation of FPM in the case with sample motion, and use mcFP to reconstruct the HR complex image and unknown motion shift.
% Ref: Liheng Bian, Guoan Zheng, et al., "Motion-corrected Fourier ptychography".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
clc;
close all;
addpath(genpath(pwd));

%% FPM simulation parameters
samplename.amplitude = 'Lena_512.png';
samplename.phase = 'Map_512.tiff';

n = 256; % pixels of high resolution image

noise.type = 'samplemotion';
noise.variance = (round(n*0.060)).^2

M_factor = 8; % magnification factor (HR/LR)

%% FPM simulation to generate captured LR images
[sample, f_sample, im_capture, fprob_real, fprob_save, dkx, dky, kx, ky, Masks, overlapratio, motionxy, NAfil] = fun_FPM_Capture(samplename, noise, n, M_factor, 0);

newfolder = ['reconstruction/Noise_' num2str(noise.type) '_' num2str(noise.variance) '_Amp_' samplename.amplitude '_Phase_' num2str(samplename.phase) '_n_' num2str(n)];
mkdir(newfolder);

save([newfolder '/sample.mat'],'sample');
save([newfolder '/f_sample.mat'],'f_sample');
save([newfolder '/im_capture.mat'],'im_capture');
save([newfolder '/fprob_real.mat'],'fprob_real');
save([newfolder '/fprob_save.mat'],'fprob_save');
save([newfolder '/dkx.mat'],'dkx');
save([newfolder '/dky.mat'],'dky');
save([newfolder '/kx.mat'],'kx');
save([newfolder '/ky.mat'],'ky');
save([newfolder '/Masks.mat'],'Masks');
save([newfolder '/overlapratio.mat'],'overlapratio');
save([newfolder '/motionxy.mat'],'motionxy');

amp_groundtruth = abs(sample);
phase_groundtruth = angle(sample);
phase_groundtruth = phase_groundtruth + pi/2;
phase_groundtruth = phase_groundtruth/pi;
imwrite(uint8(255*(amp_groundtruth)),[newfolder '/sample_amp.jpg'],'jpg');
imwrite(uint8(255*(phase_groundtruth)),[newfolder '/sample_phase.jpg'],'jpg');

%% mcFP reconstruction
shift_all_x = 3*(-sqrt(noise.variance):sqrt(noise.variance)); % users can use denser sampling for more accurate motion recovery
shift_all_y = shift_all_x;

% without motion guess (namely conventional AP)
shiftguess_flag = 0;
[im_reconst_0] = fun_mcFP(im_capture, n, M_factor, fprob_save, dkx, dky, kx, ky, NAfil, newfolder, shiftguess_flag, shift_all_x, shift_all_y, sample, f_sample, motionxy);


% with motion guess (mcFP)
shiftguess_flag = 1;
[im_reconst_1, motionxy_guess] = fun_mcFP(im_capture, n, M_factor, fprob_save, dkx, dky, kx, ky, NAfil, newfolder, shiftguess_flag, shift_all_x, shift_all_y, sample, f_sample, motionxy);

