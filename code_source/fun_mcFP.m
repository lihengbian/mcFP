function [im_reconst_compare, motionxy_guess] = fun_mcFP(im_capture, n, M_factor, fprob_save, dkx, dky, kx, ky, NAfil, newfolder, shiftguess_flag, shift_all_x, shift_all_y, sample, f_sample, motionxy, motionxyflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, Apr 26th, 2016.
% Ref: Liheng Bian, Guoan Zheng, et al., "Motion-corrected Fourier ptychography".
% Contact: lihengbian@gmail.com
% This function implements the motion-corrected Fourier ptychographic reconstruction (mcFP).

% Inputs:
% im_capture: captured images
% n: pixels of the HR sample
% M_factor: magnification factor between HR and LR
% fprob_save: pupil function
% dkx, dky, kx, ky: wave vector
% NAfil: pixel number of the radius in NA
% newfolder: folder to save figs and results
% shiftguess_flag: if guess the motion shift
% shift_all_x, shift_all_y: the search space of motion shift

% Unnecessary Inputs:
% sample, f_sample: groundtruth of HR sample and its spatial spectrum (for calculating reconstruction error)
% motionxy: groundtruth of motion shift
% motionxyflag: if using the groundtruth motion shift to reconstruct

% Outputs:
% im_reconst: recontructed HR image
% motionxy_guess: estimated motion shift for each LR image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization

alpha = 1; % AP updating parameter
iterationnum = 2; % iteration number
save_sep_AP = 5; % iterations to save

hwidth = n/M_factor/2; % half pixels of LR

f_sample_one = ones(n, n);
sample_one = ifft2(ifftshift(f_sample_one));

if size(fprob_save,1) == size(im_capture,1)
    fprob_rec = fprob_save;
else
    fprob_rec = fprob_save(n/2-hwidth+1:n/2+hwidth,n/2-hwidth+1:n/2+hwidth);
end

% creat sub-folder to save results
newfolder = [newfolder '/mcFP'];
newfolder = [newfolder '_shiftguess_' num2str(shiftguess_flag)];
if shiftguess_flag == 1 && exist('motionxyflag','var') == 1 && motionxyflag == 1
    newfolder = [newfolder '_Real'];
end
mkdir(newfolder);

% ini £¨using the LR normal image£©
temp = sum(sum(im_capture,1),2);
lr_loc = temp == max(temp);
f_reconst = fftshift(fft2(imresize(sqrt(im_capture(:,:,lr_loc)),[n, n]))); % low frequency in center
im_reconst_compare = ifft2(ifftshift(f_reconst));
amp_reconst = abs(im_reconst_compare);
phase_reconst = angle(im_reconst_compare);
phase_reconst = phase_reconst + pi/2;
phase_reconst = phase_reconst/pi;

% calculate reconstruction error
if exist('sample','var')
    AP_Relerrs_s(1) = norm(sample - exp(-1i*angle(trace(sample'*im_reconst_compare))) * im_reconst_compare, 'fro')/norm(sample,'fro'); % Initial rel. error in spatial domain
    gcf_AP_error_s = figure;
else
    AP_Relerrs_s = 0;
end


motionxy_guess = zeros(size(kx,2),2);
error_pre = 10^10*ones(size(kx,2),1);

kdis_last = 0;
i_last = 1;
j = 1;
if shiftguess_flag == 1
    newfolder_AP = [newfolder '/iteration_' num2str(j)];
    mkdir(newfolder_AP);
else
    newfolder_AP = newfolder;
end


imwrite(uint8(255*(log(abs(f_reconst))/15)),[newfolder_AP '/recover_AP_famp_0.jpg'],'jpg');
imwrite(uint8(255*(amp_reconst)),[newfolder_AP '/recover_AP_amp_0.jpg'],'jpg');
imwrite(uint8(255*(amp_reconst/max(max(amp_reconst)))),[newfolder_AP '/recover_AP_amp_norm_0.jpg'],'jpg');
imwrite(uint8(255*(phase_reconst)),[newfolder_AP '/recover_AP_phase_0.jpg'],'jpg');

%% begin updating
while 1
    %% if one loop ends, set i_last=1, and guess from the first one again
    if i_last == size(kx,2)
        % save results
        gcf_AP = figure;
        subplot(1,3,1),imshow(amp_reconst,[]); title(['AP amplitude j=' num2str(j)]); colorbar;
        subplot(1,3,2),imshow(phase_reconst,[0,1]); title('AP phase'); colorbar;
        subplot(1,3,3),imshow(log(abs(f_reconst)),[]); title('AP spatial spectrum'); colorbar;
        suptitle(['mcFP ' num2str(shiftguess_flag)]);
        saveas(gcf_AP, [newfolder_AP '/AP_R_j_' num2str(j) '.fig']);
        save([newfolder_AP '/motionxy_guess_' num2str(j) '.mat'],'motionxy_guess');
        save([newfolder '/im_reconst_' num2str(j) '_compare.mat'],'im_reconst_compare');
        
        % set to begining
        j = j + 1;
        if j>iterationnum
            break
        end
        kdis_last = 0;
        i_last = 1;
        if shiftguess_flag == 1
            error_pre = 10^10*ones(size(kx,2),1);

            newfolder_AP = [newfolder '/iteration_' num2str(j)];
            mkdir(newfolder_AP);
        end
    end
    i_last_last = i_last;
    
    %% determine next update sequence (overlap)
    if shiftguess_flag == 1
        for i = i_last + 1:size(kx,2)
            kx1 = n/2+round(kx(1,i)/dkx)-hwidth+1;
            kx2 = kx1+hwidth*2-1;
            ky1 = n/2+round(ky(1,i)/dky)-hwidth+1;
            ky2 = ky1+hwidth*2-1;

            if sqrt(round(kx(1,i)/dkx)^2+round(ky(1,i)/dky)^2) - kdis_last < NAfil
                i_last = i;
            else
                break;
            end
        end
        kdis_last = sqrt(round(kx(1,i_last)/dkx)^2+round(ky(1,i_last)/dky)^2);
    else
        i_last = size(kx,2);
    end

    %% choose the optimum motion shift guess
    if shiftguess_flag == 1
        if exist('motionxyflag','var') == 1 && motionxyflag == 1
            motionxy_guess = motionxy;
        else
            for i = i_last_last+1:i_last
                kx1 = n/2+round(kx(1,i)/dkx)-hwidth+1;
                kx2 = kx1+hwidth*2-1;
                ky1 = n/2+round(ky(1,i)/dky)-hwidth+1;
                ky2 = ky1+hwidth*2-1;

                shift_all_x_temp = shift_all_x;
                shift_all_y_temp = shift_all_y;
                for shiftnum_1 = 1:length(shift_all_x_temp)
                    for shiftnum_2 = 1:length(shift_all_y_temp)
                        x0_now = shift_all_x_temp(shiftnum_1);
                        y0_now = shift_all_y_temp(shiftnum_2);
                        addphase = fftshift(fft2(circshift(sample_one,[x0_now,y0_now])));
                        f_reconst_temp = f_reconst.*addphase;
                        phi_j = f_reconst_temp(kx1:kx2,ky1:ky2).*fprob_rec;
                        f_phi_j = ifft2(ifftshift(phi_j/M_factor^2));
                        error_now(shiftnum_1,shiftnum_2) = sum(sum(abs(abs(f_phi_j).^2 - im_capture(:,:,i))));
                    end
                end
                [shiftnum_1, shiftnum_2] = find(error_now == min(min(error_now)));
                if error_now(shiftnum_1(1),shiftnum_2(1)) < error_pre(i)
                    x0 = shift_all_x_temp(shiftnum_1(1));
                    y0 = shift_all_y_temp(shiftnum_2(1));
                    motionxy_guess(i,1) = x0;
                    motionxy_guess(i,2) = y0;
                    error_pre(i) = error_now(shiftnum_1(1),shiftnum_2(1));
                end
            end
        end
    end                
    
    %% using determined sequence (including previous ones) to reconstruct
    for k = 1:10 % update 10 times
        for i = 1:i_last

            kx1 = n/2+round(kx(1,i)/dkx)-hwidth+1;
            kx2 = kx1+hwidth*2-1;
            ky1 = n/2+round(ky(1,i)/dky)-hwidth+1;
            ky2 = ky1+hwidth*2-1;

            x0 = motionxy_guess(i,1);
            y0 = motionxy_guess(i,2);

            addphase = fftshift(fft2(circshift(sample_one,[x0,y0])));

            % compensate the motion offset
            f_reconst_temp = f_reconst.*addphase;
            phi_j = f_reconst_temp(kx1:kx2,ky1:ky2).*fprob_rec;
            f_phi_j = ifft2(ifftshift(phi_j));
            
            % correct intensity
            if k>3
                tt(i) = (mean(mean(abs(f_phi_j)))/mean(mean(M_factor^2*sqrt(im_capture(:,:,i)))));
                im_capture(:,:,i) = im_capture(:,:,i)*(tt(i)^2);
            end
            
            Phi_j = sqrt(im_capture(:,:,i)).*exp(1j*angle(f_phi_j));
            phi_j_p = fftshift(fft2(Phi_j))*M_factor^2; % scaling factor

            % update the HR spectrum
            f_dif = (phi_j_p-phi_j);
            temp = conj(fprob_rec)./(max(max((abs(fprob_rec)).^2))).*f_dif;
            f_reconst(kx1:kx2,ky1:ky2) = f_reconst(kx1:kx2,ky1:ky2)+alpha*temp .* conj(addphase(kx1:kx2,ky1:ky2));
        end
        
        % in case of global shift
        im_reconst_compare = ifft2(ifftshift(f_reconst));
        [MM, NN]=fastreg_Liheng(abs(sample), abs(im_reconst_compare));
        im_reconst_compare = circshift(im_reconst_compare,[MM, NN]);
 
        % save results
        amp_reconst = abs(im_reconst_compare);
        phase_reconst = angle(im_reconst_compare);
        phase_reconst = phase_reconst + pi/2;
        phase_reconst = phase_reconst/pi;
        
        if exist('sample','var')
            AP_Relerrs_s(1+length(AP_Relerrs_s)) = norm(sample - exp(-1i*angle(trace(sample'*im_reconst_compare))) * im_reconst_compare, 'fro')/norm(sample,'fro');

            set(groot,'CurrentFigure',gcf_AP_error_s);
            plot(AP_Relerrs_s); title(['mcFP ' num2str(shiftguess_flag) ' relative error spatial']);
            pause(0.01);
        else
            AP_Relerrs_s(1+length(AP_Relerrs_s)) = 0;
        end

        if mod(length(AP_Relerrs_s),save_sep_AP) == 0
            imwrite(uint8(255*(log(abs(f_reconst))/15)),[newfolder_AP '/recover_AP_famp_' num2str(length(AP_Relerrs_s)) '.jpg'],'jpg');
            imwrite(uint8(255*(amp_reconst)),[newfolder_AP '/recover_AP_amp_' num2str(length(AP_Relerrs_s)) '.jpg'],'jpg');
            imwrite(uint8(255*(amp_reconst/max(max(amp_reconst)))),[newfolder_AP '/recover_AP_amp_norm_' num2str(length(AP_Relerrs_s)) '.jpg'],'jpg');
            imwrite(uint8(255*(phase_reconst)),[newfolder_AP '/recover_AP_phase_' num2str(length(AP_Relerrs_s)) '.jpg'],'jpg');
            fprintf(['mcFP_' num2str(shiftguess_flag) ' j ' num2str(j) ' AP ' num2str(length(AP_Relerrs_s)) '\n']);
        end
    end

end

if exist('sample','var')
    save([newfolder '/AP_Relerrs_s.mat'],'AP_Relerrs_s');
    save([newfolder '/motionxy_guess.mat'],'motionxy_guess');
    saveas(gcf_AP_error_s, [newfolder '/AP_Relerrs_s.fig']);
end


end