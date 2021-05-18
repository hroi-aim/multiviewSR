# Copyright 2020 DRVision Technologies LLC.
# SPDX-License-Identifier: CC-BY-NC-4.0


from rcan.utils import apply, get_model_path, normalize, load_model

import argparse
import keras
import numpy as np
import tifffile
import os
from scipy.ndimage import zoom
from scipy.ndimage import rotate

#parser = argparse.ArgumentParser()
#parser.add_argument('-m', '--model_dir', type=str, required=True)
#parser.add_argument('-i', '--input', type=str, required=True)
#parser.add_argument('-o', '--output', type=str, required=True)
#parser.add_argument('-g', '--ground_truth', type=str)
#parser.add_argument('-b', '--bpp', type=int, choices=[8, 16, 32], default=32)
#args = parser.parse_args()


model_dir = 'G:\\LineConfocal\\20201211_WormNeuron\\Training\\ViewA\\RCAN_Model'
predition_dir = 'Y:\\Yicong\LineConfocal\\20201202_DCR8528_whole_worm\\worm1\\Reconstruction\\Crop_DiffractionLimit_2x'
test_folder = 'ViewA'
output_folder = test_folder + '_RCAN'
try:
    DL_path = predition_dir + '\\' + output_folder
    if not os.path.exists(DL_path):
        os.makedirs(DL_path)
except OSError:
    print ("Creation of the directory %s failed" % DL_path)
else:
    print ("Successfully created the directory %s " % DL_path)

input_labels=os.listdir(predition_dir + '\\' + test_folder)
maxlen = len(input_labels)
angle = np.linspace(-90,90, 7)

model_path = get_model_path(model_dir)
print('Loading model from', model_path)
#model = keras.models.load_model(str(model_path), compile=False)
model = load_model(str(model_path), input_shape=(32, 128, 128))

for i in range(0,maxlen):
    Predition_File = predition_dir + '\\' + test_folder + '\\'+ input_labels[i]
    print('Loading raw image from', Predition_File)
    input_data = normalize(tifffile.imread(Predition_File))
    print('Applying model')

    ndim = input_data.ndim
    if ndim == 3:
        (nz, ny, nx) = input_data.shape
        nx1 = int(np.sqrt(nx*nx + ny*ny))
        input_data_pad = np.zeros((nz,nx1,nx1))
        print(input_data_pad.shape)
        input_data_pad[:,round(nx1/2)-round(ny/2):round(nx1/2)-round(ny/2)+ny,round(nx1/2)-round(nx/2):round(nx1/2)-round(nx/2)+nx] = input_data
    else:
        (ny, nx) = input_data.shape
        nx1 = int(np.sqrt(nx*nx + ny*ny))
        input_data_pad = np.zeros((nx1,nx1))
        print(input_data_pad.shape)
        input_data_pad[round(nx1/2)-round(ny/2):round(nx1/2)-round(ny/2)+ny,round(nx1/2)-round(nx/2):round(nx1/2)-round(nx/2)+nx] = input_data       
    
    for k in range(len(angle)-1):
        if ndim == 3:
            input_data1 = rotate(input_data_pad, -1*angle[k], axes=(1,2), reshape=False)
            result = apply(model, input_data1, verbose=True)
            denoised_result = rotate(result, angle[k], axes=(1,2), reshape=False)
            final_result = denoised_result[:,round(nx1/2)-round(ny/2):round(nx1/2)-round(ny/2)+ny,round(nx1/2)-round(nx/2):round(nx1/2)-round(nx/2)+nx]
        else:
            input_data1 = rotate(input_data_pad, -1*angle[k], reshape=False)
            result = apply(model, input_data1, verbose=True)
            denoised_result = rotate(result, angle[k], reshape=False)
            #denoised_result = result
            result1 = denoised_result[round(nx1/2)-round(ny/2):round(nx1/2)-round(ny/2)+ny,round(nx1/2)-round(nx/2):round(nx1/2)-round(nx/2)+nx]

        inmin, inmax = result1.min(),result1.max()
        final_result =(result1-inmin)/(inmax-inmin)*65535
        
        #final_result = np.clip(65535 *  final_result, 0, 65535)
        tifffile.imsave(predition_dir + '\\' + output_folder + '\\DL_' + str(int(angle[k])) + '_' + input_labels[i], final_result,imagej=True)
        
        print('Saving output image', input_labels[i])
        #    result.append(normalize(tifffile.imread(args.ground_truth)))

        #result = np.stack(result)
        #if result.ndim == 4:
        #    result = np.transpose(result, (1, 0, 2, 3))


