function [stack, cropBox] = Confocal_PreProcessing_StageMode(stack,view, Step, SIM_Mode, crop_Mode, cropBox, file_path)

w = warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('off', 'MATLAB:imagesci:tifftagsread:expectedTC1_TranslatedagDataFormat');

if nargin <= 6
    SaveFlag = 0;
else 
    SaveFlag = 1;
end

if SIM_Mode == 1
    ResA = 0.1625/2; ResB = 0.1625/2; ResC = 0.0975/2;  
else
    ResA = 0.1625; ResB = 0.1625; ResC = 0.0975;    
end

 MzA = ResA/ResC;
 MzB = ResB/ResC;
 shear = Step/ResC/1.414;
 [ny, nx, nz] = size(stack); 
 
  MxA = ResA/ResC;
  MxB = ResB/ResC; 
  Mz =  shear;
 
% interpolate and rotation view A;
if view == 'ViewA'     
  tic;   
  theta = 45;
  stagescanning_matix = [1 0 0 0; 0 1 0 0; 0 -shear 1 0; 0 0 0 1];
  tform = affine3d(stagescanning_matix);
  [ny, nx, nz] = size(stack);
  stack = imresize3(stack, [ny*MxA nx*MxA nz], 'linear');   
  stack = imwarp(stack,tform);
  
  stack = flip(flip(permute(stack, [2, 1, 3]),2),3);  %1x;2y; 3z
  [ny, nx, nz] = size(stack);
  stack = imresize3(stack, [ny, nx, round(nz*Mz)], 'linear');
  stack = imrotate3(stack,theta,[0,1,0],'linear','loose'); % cannot run on gpu, out of memory issue even for 200 MB
        if crop_Mode == 1  % auto crop          
             MP = squeeze(max(stack, [], 1));
%             WriteTifStack(MP,[path_in, 'A_MP_XZ.tif'],'32');
            [MP_out, cropBox] = boundary_crop(MP); 
      %      WriteTifStack(MP_out,'F:\LineConfocal\20201002_Brain_HE\brain2\A_MP_XZ_crop.tif','32');
            
            stack = stack(:,cropBox(1):cropBox(2),cropBox(3):cropBox(4));  
        elseif crop_Mode == -1  % crop based on the box
                stack = stack(:,cropBox(1):cropBox(2),cropBox(3):cropBox(4));     
        end 
       
  disp(['ViewA pre-processing takes ', num2str(toc), ' s']);
end

% interpolate and rotation view B;
if view == 'ViewB'
    tic;
    tic;   
  theta = -45;
  stagescanning_matix = [1 0 0 0; 0 1 0 0; 0 shear 1 0; 0 0 0 1];
  tform = affine3d(stagescanning_matix);
  [ny, nx, nz] = size(stack);
  stack = imresize3(stack, [ny*MxB nx*MxB nz], 'linear');   
  stack = imwarp(stack,tform);
    
    stack = flip(flip(permute(stack, [2, 1, 3]),2),3); 
    [ny, nx, nz] = size(stack);
    stack = imresize3(stack, [ny nx round(nz*Mz)], 'linear');
          
    [ny, nx, nz] = size(stack);
     stack = imrotate3(stack,theta,[0,1,0],'linear','loose'); % cannot run on gpu, out of memory issue even for 200 MB
        if crop_Mode == 1     
            MP = squeeze(max(stack, [], 1));
            [MP_out, cropBox] = boundary_crop(MP);
             stack = stack(:,cropBox(1):cropBox(2),cropBox(3):cropBox(4));  
        elseif crop_Mode == -1
            stack = stack(:,cropBox(1):cropBox(2),cropBox(3):cropBox(4));   
        end           
    end
    
        
    disp(['ViewB pre-processing takes ', num2str(toc), ' s']);


end


