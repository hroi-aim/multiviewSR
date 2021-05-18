clear all;

% define the input and output folder of the data
input_path = 'Z:\Yicong\ConfoalLineScanning\20200802_Actin_SIM\Deep Learning\Processing\';
output_path = 'Z:\Yicong\ConfoalLineScanning\20200802_Actin_SIM\Deep Learning\Processing\6degree_decon\';

% define the PSF files for A, B, C view
File_PSF_A = 'H:\TripleViewConfocal\NewPSF20200804\PSF_A_SIM_48nm.tif';
File_PSF_B = 'H:\TripleViewConfocal\NewPSF20200804\PSF_B_SIM_48nm.tif';
File_PSF_C = 'H:\TripleViewConfocal\NewPSF20200804\PSF_C_SIM_48nm.tif';

angle_number = 6;
angles = linspace(-90,90,angle_number+1); 
angles = angles(1:end-1);

% define the file name for each view as:
% 'Reg_DL_', num2str(angles(angle)),'_A/B/C_', num2str(i),

g = gpuDevice(1); reset(g); wait(g);
g.FreeMemory

tic
for i = 1:100  % time  
data = ReadTifStack([input_path, 'Reg_DL_0_C_', num2str(i),'.tif']);
[ny, nx, nz] = size(data);
nyP = 350; nxP = 350; nzP=nz;
CropY = 30; CropX = 30; CropZ =30; OverlapY = 32; OverlapX = 32; OverlapZ=32;
TileSizeY = nyP - (2*CropY + OverlapY);
TileSizeX = nxP - (2*CropX + OverlapX);
TileY = ceil(ny/TileSizeY);
TileX = ceil(nx/TileSizeX);
nyT = (TileY-1)*TileSizeY + nyP;
nxT = (TileX-1)*TileSizeX + nxP;

disp(['there are ', num2str(TileX*TileY),' tiles!']);

iteration_number = 15; 

decon_mode = 1; % 1 for joint decon; 2 for additive decon
WB_Decon = 0;   % 1 for traditional; 0 for WB.

PSFA_1DSIM = single(ReadTifStack(File_PSF_A));
PSFB_1DSIM = single(ReadTifStack(File_PSF_B));
PSFC_1DSIM = single(ReadTifStack(File_PSF_C));
PSFA_1DSIM = align_size(PSFA_1DSIM, nyP, nxP, nzP);
PSFB_1DSIM = align_size(PSFB_1DSIM, nyP, nxP, nzP);
PSFC_1DSIM = align_size(PSFC_1DSIM, nyP, nxP, nzP);
PSFA_1DSIM  = PSFA_1DSIM /sum(PSFA_1DSIM (:));
PSFB_1DSIM  = PSFB_1DSIM /sum(PSFB_1DSIM (:));
PSFC_1DSIM  = PSFC_1DSIM /sum(PSFC_1DSIM (:));



for angle=1:length(angles)
     PSFA{angle} = imrotate3(PSFA_1DSIM, angles(angle), [-1,0, 1], 'crop');  
     PSFB{angle} = imrotate3(PSFB_1DSIM, angles(angle), [1, 0, 1], 'crop');
     PSFC{angle} = imrotate3(PSFC_1DSIM, angles(angle), [0, 0, 1], 'crop');   
end

   tic
   for angle=1:length(angles)
      disp(['loading files - angles: ', num2str(angles(angle))]);
     ImA{angle}  = align_size(max(single(ImageJ_formatted_TIFF.ReadTifStack([input_path, 'Reg_DL_', num2str(angles(angle)),'_A_', num2str(i),'.tif'])),0.01), nyT, nxT, nzP);%, 2,0.01);
     ImB{angle}  = align_size(max(single(ImageJ_formatted_TIFF.ReadTifStack([input_path, 'Reg_DL_', num2str(angles(angle)),'_B_', num2str(i),'.tif'])),0.01), nyT, nxT, nzP);%, 2,0.01);
     ImC{angle}  = align_size(max(single(ImageJ_formatted_TIFF.ReadTifStack([input_path, 'Reg_DL_', num2str(angles(angle)),'_C_', num2str(i),'.tif'])),0.01), nyT, nxT, nzP);%, 2,0.01);    
   end
   toc
   
  for TileN = 1:TileX*TileY
      tic
      disp(['processing Tile ', num2str(TileN)]);
      [x, y] = ind2sub([TileX TileY], TileN); 
      sy = (y-1)*TileSizeY+1; %start position y of Tile N  Base View
      sx = (x-1)*TileSizeX+1; %strat position x of Tile N
     
      ey = min(sy+nyP-1, nyT); %end positio y of Tile N
      ex = min(sx+nxP-1, nxT); %end positio x of Tile N 
      syT(x,y) = sy; eyT(x,y) = ey; sxT(x,y) = sx; exT(x,y) = ex;
      CropSize{x,y} = [ey-sy+1, ex-sx+1];  
             
      Estimate = zeros(nyP,nxP, nzP);
      for angle=1:length(angles)
         dataA{angle} = gpuArray(ImA{angle}(sy:ey,sx:ex,:))+0.01;
         dataB{angle} = gpuArray(ImB{angle}(sy:ey,sx:ex,:))+0.01;
         dataC{angle} = gpuArray(ImC{angle}(sy:ey,sx:ex,:))+0.01;
         Estimate = Estimate + dataA{angle} + dataB{angle} + dataC{angle};
         %Estimate = Estimate + dataC{angle};
      end
       
        Estimate =  Estimate/length(angles)/3 + 0.01;
   
    if decon_mode == 1 % joint decon
        mode = 'Joint_Decon_';
       % it=0;
        for iteration = 1:iteration_number  
            %disp(['running iteration:', num2str(iteration)]);  
             for angle=1:length(angles)
                for view = 1:3       
                    if view==1
                       data = dataA{angle};
                       PSF = gpuArray(PSFA{angle});
                     %  PSF_bp = PSFA_bp{angle};
                    elseif view == 2
                        data = dataB{angle};
                        PSF = gpuArray(PSFB{angle});
                     %   PSF_bp = PSFB_bp{angle};
                    elseif view == 3
                        data = dataC{angle};
                        PSF = gpuArray(PSFC{angle});
                       % PSF_bp = PSFC_bp{angle};
                    end
                                      
%                     OTF = fftn(circshift(gpuArray(PSF),-floor([nyP, nxP, nzP]/2)));
%                     OTF_conj = fftn(circshift(gpuArray(PSF_bp),-floor([nyP, nxP, nzP]/2)));    
                    OTF = fftn(circshift(PSF,-floor([nyP, nxP, nzP]/2)));
                    OTF_conj = fftn(circshift(flip(flip(flip(PSF,1),2),3),-floor([nyP, nxP, nzP]/2)));
                    %OTF_conj = fftn(circshift(gpuArray(flipPSF(PSF)),-floor([nyP, nxP, nzP]/2)));               
                    Temp = abs(ConvFFT3_S(Estimate, OTF));
                    Temp = data./Temp;
                    Temp = abs(ConvFFT3_S(Temp, OTF_conj));
                    Estimate = max(Estimate.*Temp,0.01);                  
                  %  g.FreeMemory
                end          
              end
        
%             data_Decon(:,:,iteration) = gather(gpu_Estimate); 
%              if mod(iteration,30) == 0
%                  it = it+1;
%                      
%                  %results = Estimate;
%                 % results = results/max(results(:))*60000;
%                 %% WriteTifStack(results, [output_path, mode,'_iteration', num2str(iteration), 'time_', num2str(i),'.tif'], '16');
%              end
        end
        Decon{x,y} = gather(Estimate); 
    end
    toc
  end
                      
    for TileN = 1: TileX*TileY
     [x, y] = ind2sub([TileX TileY], TileN);   
    if x<=TileX-1
        overlap_right{x,y} = OverlapX;
        crop_right{x,y} = CropX;
    else
        overlap_right{x,y} = 0;
        crop_right{x,y} = 0;
    end
    
    if y<=TileY-1
        overlap_bottom{x,y} = OverlapY;
        crop_bottom{x,y} = CropY;
    else
        overlap_bottom{x,y} = 0;
        crop_bottom{x,y} = 0;
    end
    
    if x>=2
        overlap_left{x,y} = OverlapX;
        crop_left{x,y} = CropX;
    else
        overlap_left{x,y} = 0;
        crop_left{x,y} = 0;
    end  
    
    if y>=2
        overlap_top{x,y} = OverlapY;
        crop_top{x,y} = CropY;
    else
        overlap_top{x,y} = 0;
        crop_top{x,y} = 0;
    end
        
    W{x,y} = Weigth(CropSize{x,y}, [crop_top{x,y}, overlap_top{x,y}], [crop_bottom{x,y}, overlap_bottom{x,y}], [crop_left{x,y}, overlap_left{x,y}], [crop_right{x,y}, overlap_right{x,y}]);  
    W_3D{x,y} =  repmat(W{x,y},[1,1, nzP]);
    end
    
  
%     for k=1:size(Decon,3)
    TileMerge = zeros(nyT, nxT, nzP);
    for TileN = 1: TileX*TileY
          [x, y] = ind2sub([TileX TileY], TileN);  
          TileMerge(syT(x,y):eyT(x,y),sxT(x,y):exT(x,y),:) = TileMerge(syT(x,y):eyT(x,y),sxT(x,y):exT(x,y),:) + Decon{x,y}.*W_3D{x,y}; % top lef is positive     
    end   
    TileMerge = align_size(TileMerge, ny, nx, nz);
    TileMerge = TileMerge/max(TileMerge(:))*60000;
      
    WriteTifStack(TileMerge,[output_path, mode,'_iteration', num2str(iteration_number), 'time_', num2str(i),'.tif'], '32');
  
  end
 toc

    
          







