function [xx_c, motionxy] = A_Function_Real_samplemotion(xx, Masks, n1_LR, n2_LR, fmaskpro, sample, variance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, Feb 4th, 2016.
% This function operates the linear transform on the signal (xx_c = Axx, Y = |xx_c|^2) with sample motion
% xx: n1 * n2, original signal (HR spatial spectrum)
% Masks: Masks: L * 2 (each point indicates the left-upper point of the LR image in the SR spectrum)
% xx_c: n1_LR * n2_LR * L, sampling output (without abs)
% fmaskpro: pupil function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = size(Masks,1);
[n1, n2] = size(xx);

xx_c = zeros(n1_LR,n2_LR,L);

motionxy = zeros(L,2);
for k = 1:L
    if k == 1
        motionxy(k,:) = [0, 0];
    else
        motionxy(k,:) = round(randn([1,2])*sqrt(variance));
% %         motionxy(k,:) = round(sqrt(variance));
    end

    sample_now = circshift(sample,motionxy(k,:));
    xx=fftshift(fft2(sample_now)); % low frequency in center
    index_x = Masks(k,1);
    index_y = Masks(k,2);
    xx_c(:,:,k) = xx(index_x:index_x+n1_LR-1 ,index_y:index_y+n2_LR-1); % low frequency in center
    xx_c(:,:,k) = xx_c(:,:,k) .* fmaskpro;
    xx_c(:,:,k) = xx_c(:,:,k) /n1/n2*n1_LR*n2_LR; % solve scaling problem by Mfactor
end

xx_c = ifftshift(xx_c,1);
xx_c = ifftshift(xx_c,2); % high frequency in center
xx_c = ifft2( xx_c );

end