%% This program was used to produce system PSF for the line confocal and SIM imaging.

clear all;
% set the path to save the simulated PSF
input_path = 'H:\TripleViewConfocal\20210415_ImageFormation\'

pixelSize = 97.5;  %nm
excWavelength = 530; % nm
detWavelength = 590; % nm
RI = 1.333; lambdaExc = excWavelength/pixelSize;lambdaDet = detWavelength/pixelSize;
NA = 1.2;  % 0.8 for top views; 1.2 for bottom view
slit_width = 5; % pixels 5 @ 97.5 nm

nx = 257; ny = 257; nz = 257;  % define the pixel size for the final PSF.
PSF_ExcLine = zeros(ny, nx, nz);

[X,Y] = meshgrid(linspace(-1,1,nx)/2,linspace(-1,1,ny)/2);
mask_NA = X.^2 + Y.^2 <= (NA/lambdaExc)^2;
w0=0.001;
mask_Gaussian = exp(-X.^2/w0^2);
mask_line = mask_NA.*mask_Gaussian;
mask_point = mask_NA;
g = 2*pi*RI/lambdaExc.*sqrt(1-(lambdaExc/RI)^2*(X.^2+Y.^2));
z = -floor(nz/2):floor(nz/2); % here z step size 1 is equivalent to pixel size
for k=1:length(z)
   PSF_ExcLine(:,:,k) = (abs(fftshift(ifft2(mask_line.*exp(1j*z(k)*g))))).^2;
end

%Producing Excitation Line
Excitation_Line_Coherent = PSF_ExcLine/max(PSF_ExcLine(:))*65535;
WriteTifStack(Excitation_Line_Coherent,[input_path, 'Test_Line_Coherent.tif'],'32'); 

%creating Detection PSF
PSF_DetPoint = zeros(ny, nx, nz);
mask_NA = X.^2 + Y.^2 <= (NA/lambdaDet)^2;
g = 2*pi*RI/lambdaDet.*sqrt(1-(lambdaDet/RI)^2*(X.^2+Y.^2));
z = -floor(nz/2):floor(nz/2); % here z step size 1 is equivalent to pixel size
for k=1:length(z)
   PSF_DetPoint(:,:,k) = (abs(fftshift(ifft2(mask_NA.*exp(1j*z(k)*g))))).^2;
end
PSF_DetPoint = PSF_DetPoint/max(PSF_DetPoint(:))*65535;

Slit = zeros(ny,nx,nz);
Slit(round(ny/2)-floor(slit_width/2):round(ny/2)+floor(slit_width/2),:,round(nz/2)) = 1000;
PSF_Confocal = abs(fftshift(ifftn(fftn(Excitation_Line_Coherent).*fftn(Slit)))).*PSF_DetPoint ;
WriteTifStack(PSF_Confocal,[input_path, 'LineConfocal_System_PSF.tif'],'32');

%loop a bead to calculate the phase Images
Beads = zeros(ny, nx);  
Beads(round(ny/2),round(nx/2)) = 1000;
phase_number = 5;
SIM_data = zeros(ny, nx, nz, phase_number);
S_2D = Slit(:,:,round(nz/2));
for k=1:nz
    L = Excitation_Line_Coherent(:,:,k);
    P = PSF_DetPoint(:,:,k);
    OTF = fft2(circshift(P,-floor([ny, nx]/2)));   
    for m =1:5
        for i=-round(ny/2)+m:phase_number:round(ny/2)
            E = imtranslate(L,[0, i]);
            S = imtranslate(S_2D, [0, i]);
            B = abs(ConvFFT3_S(E.*Beads, OTF)).*S; % blur 3D image
            SIM_data(:,:,k, m) =  SIM_data(:,:,k, m) + B;
        end
    end
end

for m = 1:5
    WriteTifStack(SIM_data(:,:,:,m),[input_path,'I_Phase', num2str(m),'.tif'],'32');
end


%calcuate the effective PSFs for 1D SIM mode
I_shrink = zeros(ny*2,nx*2,nz);
I_average = zeros(ny,nx,nz);
initial_phases = [8 7 6 5 4] + 4;
for z=1:nz
   slices = squeeze(SIM_data(:,:,z,:));
   I_shrink(:,:,z) = Shrink(slices, phase_number, initial_phases);
   I_average(:,:,z) = mean(slices,3);
end

I_shrink = I_shrink/max(I_shrink(:));
WriteTifStack(I_shrink,[input_path,'PSF_Shrink_Beads.tif'],'32');
I_average = I_average/max(I_average(:));
WriteTifStack(I_average,[input_path,'PSF_Average_Beads.tif'],'32');


function data_Blur = Blur(data, PSF)

[ny, nx, nz] = size(data);
[nyP,nxP,nzP] = size(PSF);
PSF = align_size(PSF, ny, nx, nz);

PSF = PSF/sum(PSF(:));
g = gpuDevice(1); reset(g); wait(g);

gpu_OTF = fftn(circshift(single(gpuArray(PSF)),-floor([ny, nx, nz]/2)));
gpu_data = gpuArray(data);
gpu_Blur = abs(ConvFFT3_S(gpu_data, gpu_OTF));

data_Blur = gather(gpu_Blur);          
  
end
