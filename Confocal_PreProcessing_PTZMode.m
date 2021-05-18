function [stack_out, cropBox] = Confocal_PreProcessing_PTZMode(stack,view, Step, SIM_Mode, crop_Mode, cropBox)

w = warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('off', 'MATLAB:imagesci:tifftagsread:expectedTC1_TranslatedagDataFormat');

if nargin <= 6
    SaveFlag = 0;
else 
    SaveFlag = 1;
end

if SIM_Mode == 1
    ResA = 0.1625/2; ResB = 0.1625/2; ResC = 0.0975/2;  
%     StepA = 0.25; StepB = 0.25; StepC = 0.25; %z step size in um 
else
    ResA = 0.1625; ResB = 0.1625; ResC = 0.0975;  
%     StepA = 1; StepB = 1; StepC = 0.5; %z step size in um
end

ds = 4;  % downsampling factor

  MxA = ResA/ResC;
  MxB = ResB/ResC; 
  Mz = Step/ResC;
 
if view == 'ViewC'
    tic
   [ny, nx, nz] = size(stack);
    stack_out = imresize3(stack, [ny nx round(Mz*nz)], 'linear'); 
    if crop_Mode == 1  % auto crop
            MP = squeeze(max(stack_out, [], 3));
            [stack_down_out, cropBox] = boundary_crop(MP);  
            stack_out = stack_out(cropBox(1):cropBox(2),cropBox(3):cropBox(4),:); 
        elseif crop_Mode == -1  % crop based on the box
           %stack_out = stack_out(cropBox(1):cropBox(2),cropBox(3):cropBox(4),cropBox(5):cropBox(6)); 
           if ~isempty(cropBox)
              stack_out = stack_out(cropBox(1):cropBox(2),cropBox(3):cropBox(4),:);    
           end
        end
    
    disp(['ViewC pre-processing takes ', num2str(toc), ' s']);

 end

% interpolate and rotation view A;
if view == 'ViewA'     
    stackA = flip(flip(permute(stack, [2, 1, 3]),2),3);  %1x;2y; 3z
    [ny, nx, nz] = size(stackA);
    stack_out = imresize3(stackA, [ny*MxA nx*MxA round(nz*Mz)], 'linear');
         stack_out = imrotate3(stack_out,45,[0,1,0],'linear','loose'); % cannot run on gpu, out of memory issue even for 200 MB
        if crop_Mode == 1  % auto crop
            MP = squeeze(max(stack_out, [], 1));
            [stack_down_out, cropBox] = boundary_crop(MP);  
            stack_out = stack_out(:,cropBox(1):cropBox(2),cropBox(3):cropBox(4));  
        elseif crop_Mode == -1  % crop based on the box
                stack_out = stack_out(cropBox(1):cropBox(2),cropBox(3):cropBox(4),cropBox(5):cropBox(6));     
        end
            
  
  disp(['ViewA pre-processing takes ', num2str(toc), ' s']);
end

% interpolate and rotation view B;
if view == 'ViewB'
    tic;
    stackB = flip(flip(permute(stack, [2, 1, 3]),2),3); 
    [ny, nx, nz] = size(stackB);
    stack_out = imresize3(stackB, [ny*MxB nx*MxB round(nz*Mz)], 'linear');
          
    [ny, nx, nz] = size(stackB);
 
        stack_out = imrotate3(stack_out,-45,[0,1,0],'linear','loose'); % cannot run on gpu, out of memory issue even for 200 MB
        if crop_Mode == 1
            MP = squeeze(max(stack_out, [], 1));
            [stack_down_out, cropBox] = boundary_crop(MP);
            %ImageJ_formatted_TIFF.WriteTifStack(MP,'stack_down_out.tif','16');           
            stack_out = stack_out(:,cropBox(1):cropBox(2),cropBox(3):cropBox(4));  
        elseif crop_Mode == -1
            stack_out = stack_out(cropBox(1):cropBox(2),cropBox(3):cropBox(4),cropBox(5):cropBox(6)); 
        end          
   
    
    disp(['ViewB pre-processing takes ', num2str(toc), ' s']);


end


