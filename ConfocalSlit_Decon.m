clear all;

w = warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('off', 'MATLAB:imagesci:tifftagsread:expectedTC1_TranslatedagDataFormat');

PreProcessing = 1; StepA = 1; StepB = 1; StepC=0.5; % step size for each view
Scanning_Mode = 0; % 1 stage; 0 OBJ
RegFlag = 1; % 
CoarseRegFlag = 1;
RegAllTimePoints = 1; % 1 for all time points; 0 for one matrix to apply all time points. 
TripleViewDeconFlag = 1; % 1 for triple-view decon; 0 for not to run decon -> 
DualViewDeconFlag = 0; % 1 for dual-view decon
SingleViewDeconFlag = 1;% 1 for viewcC - decon; 0 for not
WriteCropFlag = 0; % saving cropped results. 
SIM_Mode = 0; % 1 for SIM mode 49nm per pixel; 0 for diffraction limit mode; 98 nm per pixel
iteration = 30; % iteration number for WB
decon_method = 1; % 1 for addtive decon; 2 joint/alterating alteranating; 3 for additive decocn with WB; 4 for alternating decon with WB. 
 
file_path = 'I:\20210511_FlyWing_HtlGAL4_UASNLSGFP_Phalloidin647\Wing5\';
mkdir([file_path, '\Processing'])
mkdir([file_path, '\Processing\AdditiveDeconViewABC']); mkdir([file_path, 'Processing\Crop']);
mkdir([file_path, '\Processing\DeconViewC']);
mkdir([file_path, '\Processing\ViewC']);mkdir([file_path, 'Processing\ViewA']);mkdir([file_path, 'Processing\ViewB']);

libPath = 'E:\Yicong_Program\CUDA_DLL\spimfusion_DLL\';
libName = 'libapi_20200130';
% configurations
BgReg = 0;  % 8000 is normalization 
AffineRegDeconFlag = 0; % 1 for based on decon; 0 for background subtraction.
regChoice = 2; % *** registration choice: regChoice                
affMethod = 4; % affine registration method: only if regChoice == 2, 3, 4
                % 1: translation only; 
                % 2: rigid body; 
                % 3: 7 degrees of freedom (translation, rotation, scaling equally in 3 dimemsions)
                % 4: 9 degrees of freedom(translation, rotation, scaling); 
                % 5: 12 degrees of freedom; 
                % 6: rigid body first, then do 12 degrees of freedom
                % 7 ??? default one 
%int reg3d(float *h_reg, float *iTmx, float *h_img1, float *h_img2, unsigned int *imSize1, unsigned int *imSize2, int regChoice, int affMethod,
%	bool inputTmx, float FTOL, int itLimit, int deviceNum, int gpuMemMode, bool verbose, float *records);
 
if TripleViewDeconFlag ==1 | SingleViewDeconFlag==1 | DualViewDeconFlag==1 | AffineRegDeconFlag ==1
     tic;
if SIM_Mode == 1
    %SIM PSF 48 nm per pixel
    PSFA = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_A_SIM_48nm.tif'));
    PSFB = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_B_SIM_48nm.tif'));
    PSFC = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_C_SIM_48nm.tif'));
    PSFA_bp = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_A_SIM_48nm_bp.tif'));
    PSFB_bp = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_B_SIM_48nm_bp.tif'));
    PSFC_bp = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_C_SIM_48nm_bp.tif'));
else 
    PSFA = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_A_DL_97nm.tif'));
    PSFB = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_B_DL_97nm.tif'));
    PSFC = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_C_DL_97nm.tif'));
    PSFA_bp = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_A_DL_97nm_bp.tif'));
    PSFB_bp = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_B_DL_97nm_bp.tif'));
    PSFC_bp = single(ReadTifStack('H:\TripleViewConfocal\NewPSF20200804\PSF_C_DL_97nm_bp.tif'));
end
disp(['loading PSF files take ', num2str(toc), ' s']);
end

color = [488,561];
time = 1:1; %,3,7,8,9];
for t=1:length(time)
for c=1:length(color)
    
disp(['Loading raw data files.....#', num2str(time(t))]);
tic;
filenameA = ['CoarseReg_Raw_Wing5_A_', num2str(time(t)),'_',num2str(color(c)), 'nm.tif']; 
filenameB = ['CoarseReg_Raw_Wing5_B_', num2str(time(t)),'_',num2str(color(c)), 'nm.tif']; 
filenameC = ['CoarseReg_Reg_Wing5_C_', num2str(time(t)),'_',num2str(color(c)), 'nm.tif']; 
viewA = max(single(ReadTifStack([file_path, filenameA])),0);
viewB = max(single(ReadTifStack([file_path, filenameB])),0);
viewC = max(single(ReadTifStack([file_path, filenameC])),0); 

disp(['loading files take ', num2str(toc), ' s']);

if PreProcessing == 1
    tic;
%     BgA = single(ReadTifStack('F:\LineConfocal\20200722_integrin\AVG_bgA_1.tif'));
%     BgB = single(ReadTifStack('F:\LineConfocal\20200722_integrin\AVG_bgB_1.tif'));
%     BgC = single(ReadTifStack('F:\LineConfocal\20200722_integrin\AVG_bgC_1.tif'));
 %   BgA = 0; BgB =0; BgC = 0;
    BgA = 95; BgB = 95; BgC = 95;
    viewA = max(viewA - BgA,0);
    viewB = max(viewB - BgB,0);
    viewC = max(viewC - BgC,0);
    
    if t>=1 & c==1
     cropModeA = 1; cropBoxA = []; 
     cropModeB = 1; cropBoxB = []; % 1 autocrop; 0 no crop; -1 based on input. 
     cropModeC = 0; cropBoxC = []; % 1 autocrop; 0 no crop; -1 based on input. 
    else
      cropModeA = -1;  
      cropModeB = -1; 
      cropModeC = -1;
    end
%      [viewC, cropBoxC] =  boundary_crop(viewC);
  %    [viewA, cropBoxA] =  boundary_crop(viewA);
  %    [viewB, cropBoxB] =  boundary_crop(viewB);    
    % initial crop settix = 101; cropCy = 19; sizeCx = 760; sizeCy= 790;
%     cropAx = 180; cropAy = 163; sizeAx = 620; sizeAy= 420;
%     cropBx = 212; cropBy = 106; sizeBx = 620; sizeBy= 440;
%     cropCx = 1; cropCy = 1; sizeCx = 960; sizeCy= 1000;
%     viewA = viewA(cropBy:cropAy+sizeAy-1,cropAx:cropAx+sizeAx-1,:);
%     viewB = viewB(cropBy:cropBy+sizeBy-1,cropBx:cropBx+sizeBx-1,:);
%     viewC = viewC(cropCy:cropCy+sizeCy-1,cropCx:cropCx+sizeCx-1,:);
   [viewC,cropBoxC] = Confocal_PreProcessing_PZTMode(viewC, 'ViewC', StepC, SIM_Mode, cropModeC, cropBoxC);
   ImageJ_formatted_TIFF.WriteTifStack(viewC,[file_path, '\Processing\Reg_', filenameC],'16');
   
   if Scanning_Mode == 1  % moving Stage PZT
        [viewA, cropBoxA] = Confocal_PreProcessing_StageMode(viewA, 'ViewA', StepA, SIM_Mode, cropModeA, cropBoxA);
        ImageJ_formatted_TIFF.WriteTifStack(viewA,[file_path, '\Processing\ViewA\Raw_', filenameA(1:end-4), '.tif'],'16');
    	[viewB, cropBoxB] = Confocal_PreProcessing_StageMode(viewB, 'ViewB', StepB, SIM_Mode, cropModeB, cropBoxB);
        ImageJ_formatted_TIFF.WriteTifStack(viewB,[file_path, '\Processing\ViewB\Raw_', filenameB(1:end-4), '.tif'],'16');
   elseif Scanning_Mode == 0  % moving OBJ PZT 
        viewA = Confocal_PreProcessing_PZTMode(viewA, 'ViewA', StepA, SIM_Mode, cropModeA);
      %  ImageJ_formatted_TIFF.WriteTifStack(viewA,[file_path, '\Processing\ViewA\Raw_', filenameA(1:end-4), '.tif'],'16');
        viewB = Confocal_PreProcessing_PZTMode(viewB, 'ViewB', StepB, SIM_Mode, cropModeB); 
      %  ImageJ_formatted_TIFF.WriteTifStack(viewB,[file_path, '\Processing\ViewB\Raw_', filenameB(1:end-4), '.tif'],'16');
   end
       
   disp(['Data pre-processing takes ', num2str(toc), ' s']);
   
end

%TxID = fopen([file_path, '\processing\', 'record.txt'], 'a');
if CoarseRegFlag == 1     
    if t==1 & c==1
        tic;
         ds = 4;  % downsampling factor, for coarse registration - phasor
       if AffineRegDeconFlag == 1
         viewA_Down = max(SingleViewDecon(viewA(1:ds:end,1:ds:end,1:ds:end), PSFA(1:ds:end,1:ds:end,1:ds:end)),1.1);
         viewB_Down = max(SingleViewDecon(viewB(1:ds:end,1:ds:end,1:ds:end), PSFB(1:ds:end,1:ds:end,1:ds:end)),1.1);
         viewC_Down = max(SingleViewDecon(viewC(1:ds:end,1:ds:end,1:ds:end), PSFC(1:ds:end,1:ds:end,1:ds:end)),1.1);
        [RegC_Down, RegA_Down, ceA, MtxA, PhasorA_ini, PhasorA_end] = CoarseRegistration_Phasor(viewC_Down, viewA_Down);% M       
        [RegC_Down, RegB_Down, ceB, MtxB, PhasorB_ini, PhasorB_end] = CoarseRegistration_Phasor(viewC_Down, viewB_Down); 
       else
         viewA_Down = max(viewA(1:ds:end,1:ds:end,1:ds:end) - BgReg,0);%         
         viewB_Down = max(viewB(1:ds:end,1:ds:end,1:ds:end) - BgReg,0); 
         viewC_Down = max(viewC(1:ds:end,1:ds:end,1:ds:end) - BgReg,0);
         viewA_Down = max(viewA_Down - min(viewA_Down(:)),1.1);
         viewB_Down = max(viewB_Down - min(viewB_Down(:)),1.1);
         viewC_Down = max(viewC_Down - min(viewC_Down(:)),1.1);
  
        [RegC_Down, RegA_Down, ceA, MtxA, PhasorA_ini, PhasorA_end] = CoarseRegistration_Phasor(viewC_Down, viewA_Down);% M    
        disp(['coefA: ', num2str(ceA)]);
        [RegC_Down, RegB_Down, ceB, MtxB, PhasorB_ini, PhasorB_end] = CoarseRegistration_Phasor(viewC_Down, viewB_Down); 
        disp(['coefB: ', num2str(ceB)]);
 %      ImageJ_formatted_TIFF.WriteTifStack(RegA_Down ,[file_path, 'DownA_', num2str(color(c)), 'nm.tif'],'16');
 %      ImageJ_formatted_TIFF.WriteTifStack(RegB_Down ,[file_path, 'DownB_', num2str(color(c)), 'nm.tif'],'16');   
%        ImageJ_formatted_TIFF.WriteTifStack(RegC_Down ,[file_path, 'DownC_', num2str(color(c)), 'nm.tif'],'16');        
       end

        MtxA = ds*MtxA;
        MtxB = ds*MtxB;
       
%         ImageJ_formatted_TIFF.WriteTifStack(PhasorA_ini ,[file_path, 'PhasorA_ini_Down_', num2str(color(c)), 'nm.tif'],'32');
%         ImageJ_formatted_TIFF.WriteTifStack(PhasorA_end ,[file_path, 'PhasorA_end_Down_', num2str(color(c)), 'nm.tif'],'32');
%         ImageJ_formatted_TIFF.WriteTifStack(PhasorB_ini ,[file_path, 'PhasorB_ini_Down_', num2str(color(c)), 'nm.tif'],'32');
%         ImageJ_formatted_TIFF.WriteTifStack(PhasorB_end ,[file_path, 'PhasorB_end_Down_', num2str(color(c)), 'nm.tif'],'32');
          
        disp(['Coarse registration for downsampled data takes ', num2str(toc), ' s']);
    end
    
    tic;
     viewA = align_size1(viewA, size(viewC), ceil(size(viewC)/2+[MtxA(2), MtxA(1), MtxA(3)]));
     viewB = align_size1(viewB, size(viewC), ceil(size(viewC)/2+[MtxB(2), MtxB(1), MtxB(3)]));
%   
%     
      ImageJ_formatted_TIFF.WriteTifStack(viewA ,[file_path, 'CoarseReg_', filenameA(1:end-4), '.tif'],'16');%
      ImageJ_formatted_TIFF.WriteTifStack(viewB ,[file_path, 'CoarseReg_', filenameB(1:end-4), '.tif'],'16');   
      ImageJ_formatted_TIFF.WriteTifStack(viewC ,[file_path, 'CoarseReg_', filenameC(1:end-4), '.tif'],'16');
 
    disp(['Coarse registration for Whole data takes ', num2str(toc), ' s']);
end

if RegFlag==1 | TripleViewDeconFlag==1 | DualViewDeconFlag == 1 | SingleViewDeconFlag == 1;


[nyA, nxA, nzA] = size(viewA);
[nyB, nxB, nzB] = size(viewB);
[nyC, nxC, nzC] = size(viewC);

% viewA = viewA/mean(viewA, 'all')*1.0;
% viewB = viewB/mean(viewB, 'all')*1.0;
% viewC = viewC/mean(viewC, 'all')*1.0;
%  ma  = max(viewA(:)); mb  = max(viewB(:)); mc  = max(viewC(:)); mm = max(max(ma,mb),mc);
%   viewA = viewA/ma * 10000 + 1;
%   viewB = viewB/mb * 10000 + 1;
%   viewC = viewC/mc * 10000 + 1;

% definding tiles of number


VolumeTH = 1.3; % 1.3 GB for registration only; 0.5 GB for triple-view decon; 1GB for dual-view decon
FinalVolume = VolumeTH; 
fileSize = nyC*nxC*nzC*2/10^9; %(in GB, 16bit)
TileN = ceil(fileSize/VolumeTH); 

while FinalVolume >=VolumeTH %   
% TileX = 1; TileY=5; TileZ=1;
if TileN == 1
    TileX = 1; TileY = 1; TileZ = 1;
else
    %TileZ = max(floor(nzC/350),1);   % maximal z is set to 512]
    TileZ = 1;
    TileXY = ceil(TileN/TileZ);
    ratio = nxC/nxC;
    if ratio>=1.5
        TileX = max(ceil(ratio*sqrt(TileXY), TileXY));
        TileY = ceil(TileXY/TileX);
    elseif ratio<=0.5
       TileY = max(ceil(1/ratio*sqrt(TileXY), TileXY));
       TileX = ceil(TileXY/TileY);
    else
        TileX = ceil(sqrt(TileXY));
        TileY = ceil(TileXY/TileX);
    end  
end

TileN = TileX * TileY * TileZ;
   
if TileZ == 1
    CropZ = 0; OverlapZ = 0;
    CropSizeZ = nzC;
    TileSizeZ = nzC;
else
    CropZ = 40; OverlapZ = 40;
    TileSizeZ = ceil(nzC/TileZ);
    CropSizeZ = TileSizeZ + 2*CropZ + OverlapZ;
end

if TileX == 1
    CropX = 0; OverlapX = 0;
    CropSizeX = nxC
    TileSizeX = nxC;
else
    CropX = 40; OverlapX = 40;
    TileSizeX = ceil(nxC/TileX);
    CropSizeX = TileSizeX + 2*CropX + OverlapX;
end

if TileY == 1
    CropY = 0; OverlapY = 0;
    CropSizeY = nyC
    TileSizeY = nyC;
else
    CropY = 40; OverlapY = 40;
    TileSizeY = ceil(nyC/TileY);
    CropSizeY = TileSizeY + 2*CropY + OverlapY;
end

FinalVolume = CropSizeX*CropSizeY*CropSizeZ*2/(10^9);

if FinalVolume >=VolumeTH
    TileN = TileN + 1;
end

end

disp(['Total tile number is ', num2str(TileN) '; splitting into (xyz): ', num2str(TileX), 'x', num2str(TileY), 'x', num2str(TileZ)]);
disp(['Each volume ', num2str(CropSizeX), 'x', num2str(CropSizeY), 'x', num2str(CropSizeZ)]);
disp(['Final sub-volume size is  ', num2str(FinalVolume), '  GB']);
 
% TxtID = fopen([file_path, '\processing\', 'record.txt'], 'a');
% fprintf(TxtID , ['Crop into Tile Number: X', num2str(TileX), '_Y', num2str(TileY), '_Z', num2str(TileZ), '\n']);
% fclose(TxtID);
   
 for TileN = 1:TileX*TileY*TileZ
    tic
    disp(['processing Tile ', num2str(TileN)]);
    [x, y, z] = ind2sub([TileX TileY TileZ], TileN); 
    
    syC = (y-1)*TileSizeY+1; %start position y of Tile N  Base View
    sxC = (x-1)*TileSizeX+1; %strat position x of Tile N
    szC = (z-1)*TileSizeZ+1; 
    eyC = min((y-1)*TileSizeY+CropSizeY, nyC); %end positio y of Tile N
    exC = min((x-1)*TileSizeX+CropSizeX, nxC); %end position of Tile N 
    ezC = min((z-1)*TileSizeZ+CropSizeZ, nzC);
    cyC = round((eyC+syC-1)/2); % center position y of Tile N
    cxC = round((exC+sxC-1)/2); % center position x of Tile N
    czC = round((ezC+szC-1)/2);
    CropSize{x,y,z} = [eyC-syC+1, exC-sxC+1,ezC-szC+1];
    sy(x,y,z) = syC; ey(x,y,z) = eyC; sx(x,y,z) = sxC; ex(x,y,z) = exC;
    sz(x,y,z) = szC; ez(x,y,z) = ezC;
    
    cA = [cyC, cxC, czC, 1]; % *oTmxA
    syA = max(cA(1) - round(CropSizeY/2) + 1, 1);    
    sxA = max(cA(2) - round(CropSizeX/2) + 1,  1); 
    szA = max(cA(3) - round(CropSizeZ/2) + 1, 1); 
    eyA = min(cA(1) + round(CropSizeY/2), nyA);
    exA = min(cA(2) + round(CropSizeX/2), nxA);
    ezA = min(cA(3) + round(CropSizeZ/2), nzA); 
    
    cB = [cyC, cxC, czC, 1]; % *oTmxA
    syB = max(cB(1) - round(CropSizeY/2) + 1, 1);
    sxB = max(cB(2) - round(CropSizeX/2) + 1, 1); 
    szB = max(cB(3) - round(CropSizeZ/2) + 1, 1); 
    eyB = min(cB(1) + round(CropSizeY/2), nyB);
    exB = min(cB(2) + round(CropSizeX/2), nxB);
    ezB = min(cB(3) + round(CropSizeZ/2), nzB); 
    subVolumeC = viewC(syC:eyC, sxC:exC,szC:ezC); 
    subVolumeA = viewA(syC:eyC, sxC:exC,szC:ezC);     
    subVolumeB = viewB(syC:eyC, sxC:exC,szC:ezC); 
%   subVolumeA = viewA(syA:eyA, sxA:exA,szA:ezA); 
%   subVolumeB = viewB(syB:eyB, sxB:exB,szB:ezB); 
% subVolumeA = viewA;
% subVolumeB = viewB;

     if RegFlag == 1 
         if (t==1 & c==1) | (RegAllTimePoints == 1)
           % regChoice = 2;  affMethod = 7;          
            % based on decon to register
            if AffineRegDeconFlag == 1
             Decon_subVolumeA = SingleViewDecon_WB(subVolumeA, PSFA, PSFA_bp); 
             Decon_subVolumeB = SingleViewDecon_WB(subVolumeB, PSFB, PSFB_bp); 
             Decon_subVolumeC = SingleViewDecon_WB(subVolumeC, PSFC, PSFC_bp); 
             [RegA, oTmxA{x,y,z}, nccA, mmrA] = reg3d_CUDA(Decon_subVolumeC, Decon_subVolumeA, libPath, libName, regChoice, affMethod);
             [RegB, oTmxB{x,y,z}, nccB, mmrB] = reg3d_CUDA(Decon_subVolumeC, Decon_subVolumeB, libPath, libName, regChoice, affMethod); 
            else
             [RegA, oTmxA{x,y,z}, nccA, mmrA] = reg3d_CUDA(subVolumeC-BgReg, subVolumeA-BgReg, libPath, libName, regChoice, affMethod);
             [RegB, oTmxB{x,y,z}, nccB, mmrB] = reg3d_CUDA(subVolumeC-BgReg, subVolumeB-BgReg, libPath, libName, regChoice, affMethod);
%              TxtID = fopen([file_path, '\processing\', 'record.txt'], 'a');
%             fprintf(TxtID , ['Fine Regisration A TileX', num2str(x), ' TileY', num2str(y), '- NCC: ', num2str(nccA), '  Memory:', num2str(mmrA),'\n']);
%             fprintf(TxtID , ['Fine Regisration B TileX', num2str(x), ' TileY', num2str(y), '- NCC: ', num2str(nccB), '  Memory:', num2str(mmrB),'\n']);
%             fclose(TxtID); 
            end
         end
            [subVolumeRegA{x,y,z}, oTmxA{x,y,z}, nccA, mmrA] = reg3d_CUDA(subVolumeC, subVolumeA, libPath, libName, 0, 0, 1, oTmxA{x,y,z});
            [subVolumeRegB{x,y,z}, oTmxB{x,y,z}, nccB, mmrB] = reg3d_CUDA(subVolumeC, subVolumeB, libPath, libName, 0, 0, 1, oTmxB{x,y,z});
             subVolumeRegC{x,y,z} = subVolumeC;
  
     elseif  RegFlag == 0
             subVolumeRegC{x,y,z} = subVolumeC; 
             subVolumeRegA{x,y,z} = subVolumeA; 
             subVolumeRegB{x,y,z} = subVolumeB; 
     end
     
    if WriteCropFlag == 1
           ImageJ_formatted_TIFF.WriteTifStack(subVolumeRegA{x,y,z},[file_path, '\Processing\Crop\RawA_X',num2str(x),'_Y', num2str(y),'_Z', num2str(z), '.tif'],'16');
           ImageJ_formatted_TIFF.WriteTifStack(subVolumeRegB{x,y,z},[file_path, '\Processing\Crop\RawB_X',num2str(x),'_Y', num2str(y),'_Z', num2str(z), '.tif'],'16');
           ImageJ_formatted_TIFF.WriteTifStack(subVolumeRegC{x,y,z},[file_path, '\Processing\Crop\RegC_X',num2str(x),'_Y', num2str(y),'_Z', num2str(z), '.tif'],'16');
    end
     
    if TripleViewDeconFlag == 1 
          switch decon_method
                 case 1
                     stackABC_Decon{x,y,z} = TripleViewAdditiveDecon(subVolumeRegA{x,y,z},subVolumeRegB{x,y,z}, subVolumeRegC{x,y,z},PSFA, PSFB, PSFC,1, iteration);
                 case 2
                     stackABC_Decon{x,y,z} = TripleViewJointDecon(subVolumeRegA{x,y,z},subVolumeRegB{x,y,z}, subVolumeRegC{x,y,z},PSFA, PSFB, PSFC,1, iteration);
                 case 3
                     stackABC_Decon{x,y,z} = TripleViewAdditiveDecon_WB(subVolumeRegA{x,y,z},subVolumeRegB{x,y,z}, subVolumeRegC{x,y,z},PSFA, PSFB, PSFC, PSFA_bp, PSFB_bp, PSFC_bp, 1, iteration);
                 case 4
                     stackABC_Decon{x,y,z} = TripleViewJointDecon_WB(subVolumeRegA{x,y,z},subVolumeRegB{x,y,z}, subVolumeRegC{x,y,z},PSFA, PSFB, PSFC, PSFA_bp, PSFB_bp, PSFC_bp, 1, iteration);
          end
     end
               
    if SingleViewDeconFlag == 1
         stackC_Decon{x,y,z} = SingleViewDecon(subVolumeRegC{x,y,z}, PSFC); 
%          stackA_Decon{x,y,z} = SingleViewDecon(subVolumeRegA{x,y,z}, PSFA); 
%          stackB_Decon{x,y,z} = SingleViewDecon(subVolumeRegB{x,y,z}, PSFB); 
    end
      
    if DualViewDeconFlag == 1
         stackAB_Decon{x,y,z} = DualViewJointDecon(subVolumeRegA{x,y,z},subVolumeRegB{x,y,z}, PSFA, PSFB, 1, 10);
    end
        
    if WriteCropFlag == 1 & TripleViewDeconFlag == 1
         ImageJ_formatted_TIFF.WriteTifStack(stackABC_AdditiveDecon{x,y,z},[file_path, '\Processing\TripleViewDecon_X',num2str(x),'_Y', num2str(y),'_Z', num2str(z),'.tif'],'16'); 
    end
    
  
    if WriteCropFlag == 1 & SingleViewDeconFlag == 1
        ImageJ_formatted_TIFF.WriteTifStack(stackC_Decon{x,y,z},[file_path, '\Processing\SingleViewDecon_X',num2str(x),'_Y', num2str(y),'_Z', num2str(z), '.tif'],'16'); 
    end
    
end

disp(['Registration/deconvolution/writing time for this tile takes ', num2str(toc), ' s']);

tic
  %W = Weigth(nsize, top, bottom, left, right)  [crop, overlap];
clear viewA; clear viewB; clear viewC;  % clean variable, free memroy

for TileN = 1: TileX*TileY*TileZ
[x, y, z] = ind2sub([TileX TileY TileZ], TileN);
clear subVolumeRegA{x,y,z};
clear subVolumeRegB{x,y,z};
clear subVolumeRegC{x,y,z};
end

for TileN = 1: TileX*TileY*TileZ
    TileN
     [x, y, z] = ind2sub([TileX TileY TileZ], TileN); 
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
    
    if z == 1
        z_overlap_top{z} = 0;
        z_crop_top{z} = 0;
        z_overlap_bottom{z} = OverlapZ;
        z_crop_bottom{z} = CropZ;
    elseif z == TileZ
        z_overlap_top{z} = OverlapZ;
        z_crop_top{z} = CropZ;
        z_overlap_bottom{z} = 0;
        z_crop_bottom{z} = 0;
    else
        z_overlap_top{z} = OverlapZ;
        z_crop_top{z} = CropZ;
        z_overlap_bottom{z} = OverlapZ;
        z_crop_bottom{z} = CropZ;
    end     
   
    
    W{x,y} = Weigth(CropSize{x,y,z}, [crop_top{x,y}, overlap_top{x,y}], [crop_bottom{x,y}, overlap_bottom{x,y}], [crop_left{x,y}, overlap_left{x,y}], [crop_right{x,y}, overlap_right{x,y}]);  
    W_3D{x,y,z} =  repmat(W{x,y},[1,1, CropSize{x,y,z}(3)]);
    W_Z{x,y,z} = WeightZ(CropSize{x,y,z},[z_crop_top{z}, z_overlap_top{z}],[z_crop_bottom{z},z_overlap_bottom{z}]);
    W_3D{x,y,z} = W_3D{x,y,z}.*W_Z{x,y,z};
    %ImageJ_formatted_TIFF.WriteTifStack(W_3D{x,y,z},[file_path, '\Processing\W_x', num2str(x),'_y', num2str(y), '_z', num2str(z), '.tif'],'32');
end

disp(['Creating weight image for stictching takes ', num2str(toc), ' s']);

if RegFlag == 1
    tic
       %stitching ViewA
        TileMerge = zeros(nyC,nxC, nzC);
             for TileN = 1: TileX*TileY*TileZ
                [x, y, z] = ind2sub([TileX TileY TileZ], TileN);      
                 TileMerge(sy(x,y,z):ey(x,y,z),sx(x,y,z):ex(x,y,z),sz(x,y,z):ez(x,y,z)) = TileMerge(sy(x,y,z):ey(x,y,z),sx(x,y,z): ex(x,y,z),sz(x,y,z):ez(x,y,z)) + subVolumeRegA{x,y,z}.*W_3D{x,y,z}; % top lef is positive     
             end   
        WriteTifStack(TileMerge,[file_path, '\Processing\Reg_', filenameA],'16'); 
         
       %stitching ViewB
        TileMerge = zeros(nyC,nxC, nzC);
             for TileN = 1: TileX*TileY*TileZ
                [x, y, z] = ind2sub([TileX TileY TileZ], TileN);      
                 TileMerge(sy(x,y,z):ey(x,y,z),sx(x,y,z):ex(x,y,z),sz(x,y,z):ez(x,y,z))= TileMerge(sy(x,y,z):ey(x,y,z),sx(x,y,z): ex(x,y,z),sz(x,y,z):ez(x,y,z)) + subVolumeRegB{x,y,z}.*W_3D{x,y,z}; % top lef is positive     
             end   
        WriteTifStack(TileMerge,[file_path, '\Processing\Reg_', filenameB],'16'); 
       
   disp(['Writing stitched, registered A and B views take ', num2str(toc), ' s']);
end

 tic
if TripleViewDeconFlag == 1
       %%stitching triple-view decon
        TileMerge = zeros(nyC,nxC, nzC);
             for TileN = 1: TileX*TileY*TileZ
                [x, y, z] = ind2sub([TileX TileY TileZ], TileN);      
                 TileMerge(sy(x,y,z):ey(x,y,z),sx(x,y,z):ex(x,y,z),sz(x,y,z):ez(x,y,z)) = TileMerge(sy(x,y,z):ey(x,y,z),sx(x,y,z): ex(x,y,z),sz(x,y,z):ez(x,y,z)) + stackABC_Decon{x,y,z}.*W_3D{x,y,z}; % top lef is positive     
             end
             
        WriteTifStack(TileMerge,[file_path, '\Processing\AdditiveDeconViewABC\TripleViewDecon_C', num2str(color(c)),'_t', num2str(time(t)),'.tif'],'32');
end

if DualViewDeconFlag == 1
       TileMerge = zeros(nyC,nxC, nzC);
             for TileN = 1: TileX*TileY*TileZ
                 [x, y, z] = ind2sub([TileX TileY TileZ], TileN);      
                 TileMerge(sy(x,y,z):ey(x,y,z),sx(x,y,z):ex(x,y,z),sz(x,y,z):ez(x,y,z)) = TileMerge(sy(x,y,z):ey(x,y,z),sx(x,y,z): ex(x,y,z),sz(x,y,z):ez(x,y,z)) + stackAB_Decon{x,y,z}.*W_3D{x,y,z}; % top lef is positive     
             end   
        WriteTifStack(TileMerge,[file_path, '\Processing\DualViewDecon_', filenameC(1:end-4), '.tif'],'16'); 
end

         %stitching Decon ViewC
if SingleViewDeconFlag == 1   
        TileMerge = zeros(nyC,nxC, nzC);
             for TileN = 1: TileX*TileY*TileZ
                 [x, y, z] = ind2sub([TileX TileY TileZ], TileN);      
                 TileMerge(sy(x,y,z):ey(x,y,z),sx(x,y,z):ex(x,y,z),sz(x,y,z):ez(x,y,z)) = TileMerge(sy(x,y,z):ey(x,y,z),sx(x,y,z): ex(x,y,z),sz(x,y,z):ez(x,y,z)) + stackC_Decon{x,y,z}.*W_3D{x,y,z}; % top lef is positive     
             end   
        WriteTifStack(TileMerge,[file_path, '\Processing\DeconViewC\SingleViewDecon_C', num2str(color(c)),'_t', num2str(time(t)),'.tif'],'32'); 
end
    disp(['Writing stitched, registered, view C decon and triple-view decon take ', num2str(toc), ' s']);
end

end
end
