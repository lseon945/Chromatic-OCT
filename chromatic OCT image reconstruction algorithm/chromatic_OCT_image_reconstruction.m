clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chromatic OCT Version
% 
% v20241010 drafted by Seung Eon Lee
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters
n_alines = 2000; % number of alines
galvoshift = 1786; 
axial_pixel_resol = 0.7251; % sample in water

%%% Frame view
frame_num = 32;

%%% OPTIONS
ChG = true; % true: ChG algorithm, false: FT
contrast = [53 93]; 

hann_win = false;   window = hann(4096);

additional_discom_2nd_coeff = 0;  
additional_discom_3rd_coeff = 0;


%%%% INI Parameters
%%% Filename
filelist = dir;
for i = 1 : length(filelist)
    if (length(filelist(i).name) > 5)
        if strcmp(filelist(i).name(end-3:end),'.ini')
            dfilenm = filelist(i).name(1:end-4);
        end
    end
end
clear filelist;

%%% Please check ini file in the folder
fid = fopen(strcat(dfilenm,'.ini'));
config = textscan(fid,'%s');
fclose(fid); 

for i = 1 : length(config{:})
    if (strfind(config{1}{i},'acqHeight'))
        eq_pos = strfind(config{1}{i},'=');
        n_alines0 = str2double(config{1}{i}(eq_pos+1:end));
    end
    if (strfind(config{1}{i},'acqWidth'))
        eq_pos = strfind(config{1}{i},'=');
        n_scans = str2double(config{1}{i}(eq_pos+1:end));
    end
end

% Size parameters
n_len = 2^ceil(log2(n_scans));
n2_len = n_len/2;

clear fid config eq_pos x y theta rho term i;

%%%% Background data from --.background
bgname = [dfilenm,'.background'];
fid = fopen(bgname,'r','l');
oct_bg = fread(fid,'uint16');
fclose(fid);

bgmm = reshape(oct_bg(1:n_scans*n_alines0),n_scans,n_alines0);
bgm = mean(bgmm(:,1:20), 2);
clear bgname fid oct_bg oct_bg;

%%%% K-linearization coefficeint data from spectrometer company
pixel = 1:1:4096;
C0 = 640.757;
C1 = 0.0830411;
C2 = -1.00828e-06;
C3 = -7.54088e-11;
lambda_p = C0 + C1*pixel + C2*pixel.^2 + C3*pixel.^3;
k_p = 2*pi./lambda_p;
k_lin = linspace(2*pi/641, 2*pi/959, 4096);
clear C0 C1 C2 C3 lambda_p 

%%%% Data for dispersion compensation from d1 mirror
name = 'd1.bin';
fid = fopen(name,'r','l');
d1_data = fread(fid,'uint16');
fclose(fid);

d1=reshape(d1_data,n_scans,n_alines0);

I_s = imfilter(d1,fspecial('average', [100 1]),'replicate');
d1_nobg = d1 - I_s;
d1_nobg = fft(d1_nobg, n_len);
d1_nobg(n2_len+1:end,:) = 0; % mirror image   
d1_nobg = ifft(d1_nobg);    
d1_klin = interp1(k_p, d1_nobg, k_lin, 'pchip'); % k-linearization of d1 data

phase = unwrap(angle(d1_klin));
phase = mean(phase,2);
index = transpose(1:1:4096);
p=polyfit(index,phase,1);
discom_ph = polyval(p, index);

dispersion_d1 = phase - discom_ph; 

clear fid d1 d1_data d1_nobg d1_klin I_s phase discom_ph p

%%%% Dispersion compensation
discomidex=transpose(1:1:n_scans);
discom2=exp(1i*additional_discom_2nd_coeff*((discomidex-n_scans/2)/n_scans).^2);
discom3=exp(1i*additional_discom_3rd_coeff*((discomidex-n_scans/2)/n_scans).^3);

discom = exp(-1i*dispersion_d1);
dsp = discom .* discom2 .* discom3;

clear discomidex discom discom2 discom3 dispersion_d1 dispersion_from_already_known_p2 dispersion_from_d1_mirror only_additional_discom dsp_re dsp_im bgmm bg0 bg1 bg2 bg1_im bg1_re f0 f1 i;

%%%% Data Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname11 = sprintf([dfilenm,'.data']);
fid_dx = fopen(fname11,'r','l');
fseek(fid_dx,0,'eof');
n_frame = floor(ftell(fid_dx)/2/n_scans/n_alines);


%%% load data
fid_dx = fopen(fname11,'r','l');
fseek(fid_dx, (frame_num-1)*2*n_scans*n_alines, 'bof');
ft = fread(fid_dx,n_scans*n_alines,'uint16=>single');  
clear fid_dx
        
%%% Reshape data
f_oct_raw = reshape(ft,n_scans,n_alines);
fringe_raw = zeros(4096, size(f_oct_raw,2));
fringe_raw(1:4:4096,:)=f_oct_raw(1:2:2048,:);
fringe_raw(2:4:4096,:)=f_oct_raw(2:2:2048,:);
fringe_raw(3:4:4096,:)=f_oct_raw(2049:2:end,:);
fringe_raw(4:4:4096,:)=f_oct_raw(2050:2:end,:);
clear ft

%%% BG removal 
bg = repmat(bgm, 1, n_alines);
data_signal = fringe_raw - bg;
clear fringe_raw

%%% I_s_removal
I_s = imfilter(data_signal,fspecial('average', [51 1]),'replicate');
data_signal = data_signal - I_s;
clear I_s

%%% hanning window   
if hann_win == true
    data_signal = [data_signal .* repmat(window,[1 n_alines]); zeros(n_len-n_scans,n_alines)];
else
    data_signal = [data_signal ; zeros(n_len-n_scans,n_alines)];
end

%%% Complex fringe
dft_signal = fft(data_signal, n_len);       % FFT real data
dft_signal(n2_len+1:end,:) = 0;         % mirror image removal
scan_signal_zoom1 = ifft(dft_signal);       % IFFT to get complex fringes
clear dft_signal data_signal

%%% K-linearization
dft_complex1 = interp1(k_p, scan_signal_zoom1, k_lin, 'pchip');

%%% Dispersion compensation
discom_final = repmat(dsp,[1 n_alines]);
fringe_final = discom_final.*dft_complex1;
clear discom_final dft_complex1

%% Fourier transform

if ChG == false
    %%% FFT & log to make final image
    image_linear = fft(fringe_final, n_len);
    image_linear = image_linear(1:n2_len,:);
    
    image_dB = 10*log10(abs(image_linear).^2);    
    
    clear scan_signal_zoom1 dft_complex* a_line_linear;   
    
    %%% Scaling 8 bit image
    xyim1 = 255 * mat2gray(image_dB,contrast);
    xyim1 = circshift(xyim1',-galvoshift)';
    xyim2 = medfilt2(xyim1,[3 3],'symmetric');
    xyim2 = uint8(xyim2);
    clear xyim1;
    
    %%% Show image
    figure(); set(gcf,'Position',[420 180 516 778]); 
    imshow(xyim2);
    title("ORIGINAL, contrast ["+num2str(contrast(1)) +", "+num2str(contrast(2)) + "]");
    
    fclose('all');    
    clear frame;
end


%% ChG algorithm

if ChG == true
    z_lambda_c = 500; % z position of the focus of center wavelength relative to reference mirror (unit: um)
    z_delta_lambda = 646.7017; % chromatic focal shift (unit: um)
    gaussian_radius_ChG3 = 840.69; % width of gaussian window

    corr_freq = ((0:n2_len-1) * 0.5 / n2_len)';
    corr_freq = repmat(corr_freq,1,n_len);

    corr_win = repmat(1:n_len,n2_len,1);
    corr_win = exp(1i*2*pi*corr_freq.*corr_win);
    
    z_lambda1 = z_lambda_c - z_delta_lambda/2;
    zi_lambda1 = (z_lambda1 / axial_pixel_resol);
    zi_delta_lambda = (z_delta_lambda / axial_pixel_resol);

    zi = (1:n2_len)';
    pi_window_center = 4096/(zi_delta_lambda) * (zi - zi_lambda1);

    edge_size = gaussian_radius_ChG3; 
    pi_window_center(pi_window_center > 4096-edge_size) = 4096-edge_size;
    pi_window_center(pi_window_center < edge_size) = edge_size;

    [Pi, ~] = meshgrid(1:n_len, 1:n2_len);

    gaussian_win_ChG = exp((-(Pi-pi_window_center).^2)./gaussian_radius_ChG3.^2);

    corr_win_ChG3 = corr_win.*gaussian_win_ChG;

    image_linear_ChG3 = zeros(n2_len,n_alines);
    fringe_final_single = single(fringe_final);
    corr_win_ChG3_r = single(real(corr_win_ChG3));
    corr_win_ChG3_r_m = -corr_win_ChG3_r;
    corr_win_ChG3_i = single(imag(corr_win_ChG3));
    fringe_final_r = single(real(fringe_final));
    fringe_final_i = single(imag(fringe_final));

    for i_aline = 1:n_alines
        Aline_ChG3_r = corr_win_ChG3_r * fringe_final_r(:,i_aline) + corr_win_ChG3_i * fringe_final_i(:,i_aline);
        Aline_ChG3_i = corr_win_ChG3_r_m * fringe_final_i(:,i_aline) + corr_win_ChG3_i * fringe_final_r(:,i_aline);
        Aline_ChG3 = Aline_ChG3_r + 1i*Aline_ChG3_i;
        image_linear_ChG3(:,i_aline) = Aline_ChG3;
    end

    image_dB_ChG = 10*log10(abs(image_linear_ChG3).^2);

    xyim1 = 255 * mat2gray(image_dB_ChG,contrast);
    xyim1 = circshift(xyim1',-galvoshift)';
    xyim2 = medfilt2(xyim1,[3 3],'symmetric');
    clear xyim1;

    %%% Show image
    xyim2=uint8(xyim2);
    figure(); set(gcf,'Position',[420 180 516 778]); 
    imshow(xyim2);
    title("ChG algorithm, contrast ["+num2str(contrast(1)) +", "+num2str(contrast(2)) + "]");

    clear frame;                
end
