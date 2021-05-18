 function [I_Shrink, I_Average] = ConfocalSIM_Rescale(data, zoom_flag)
 if nargin <=1
     zoom_flag =0 ;
 end
 [ny, nx, nz] = size(data);
 phase_number = 5;
 
for k=1:phase_number 
    I_Phase{k} = max(data(:,:,k:phase_number:end),[],3);   
    Intensity{k} = sum(I_Phase{k},2);     
    Phase(k)  = find_phase(Intensity{k}, phase_number) + phase_number;
end
       
 for k=1:phase_number
    if Phase(k) > Phase(1)
        Phase(k) = Phase(k) - phase_number; 
    end
    Phase(k) = Phase(k) + k - 1;
 end

 PhaseM = mean(Phase);
 P = [PhaseM, PhaseM-1, PhaseM-2, PhaseM-3, PhaseM-4];
 
 I_Shrink = zeros(2*ny, 2*nx, nz/phase_number);
 if zoom_flag == 1
    I_Average = zeros(2*ny, 2*nx, nz/phase_number); 
 else
     I_Average = zeros(ny, nx, nz/phase_number);
 end

for z=1:nz/phase_number
    slices = data(:,:,(z-1)*phase_number+1:z*phase_number); 
    I_Shrink(:,:,z) = Shrink(slices, phase_number, P);
    if zoom_flag == 1
        I_Average(:,:,z) = imresize(squeeze(sum(slices,3)),[2*ny, 2*nx],'bilinear');
    else
        I_Average(:,:, z) = sum(slices,3);
    end
end
% 
