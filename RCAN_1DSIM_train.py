# Copyright 2020 DRVision Technologies LLC.
# SPDX-License-Identifier: CC-BY-NC-4.0

from rcan.data_generator import DataGenerator
from rcan.losses import mae
from rcan.metrics import psnr
from rcan.model import build_rcan
from rcan.utils import normalize, staircase_exponential_decay

import argparse
import functools
import itertools
import json
import jsonschema
import keras
import numpy as np
import pathlib
import tifffile
import tqdm
import tqdm.keras
import sys
import os


def load_data_dir(input_dir,gt_dir):
    data = []
    input_labels = os.listdir(input_dir)
    gt_labels = os.listdir(gt_dir)
    
    maxlen_input = len(input_labels)
    maxlen_gt = len(gt_labels)

    if maxlen_input != maxlen_gt:
        raise ValueError(
            'the number of raw and GT images are not the same')
    for i in range(0,maxlen_input):
            if input_labels[i] [-4:]!= gt_labels[i][-4:]:
                raise ValueError(
                    'the name of raw and GT images are not the same')
        
            raw = tifffile.imread(input_dir + '\\' + input_labels[i]);
            gt =  tifffile.imread(gt_dir + '\\' + gt_labels[i]);
            if raw.shape != gt.shape:
                raise ValueError(
                        'Raw and GT images must be the same size: ')          
            data.append([normalize(m) for m in [raw, gt]])

    return data

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
output_dir = 'D:\YicongWu\LineConfocalData\Ivans_cell\EB3-TdTomatog\RCAN_Model'
input_dir = 'D:\YicongWu\LineConfocalData\Ivans_cell\EB3-TdTomato\Input'
gt_dir = 'D:\YicongWu\LineConfocalData\Ivans_cell\EB3-TdTomato\GT'

epochs = 300;
steps_per_epoch = 400;

num_channels = 32
num_residual_blocks = 5
num_residual_groups = 5
channel_reduction = 8

print('Loading training data')

training_data = load_data_dir(input_dir,gt_dir)
validation_data = None
ndim = training_data[0][0].ndim
input_shape = (16, 256, 256) if ndim == 3 else (256, 256)

for p in itertools.chain(training_data, validation_data or []):
    if p[0].ndim != ndim:
        raise ValueError('All images must have the same number of dimensions')

for p in itertools.chain(training_data, validation_data or []):
    input_shape = np.minimum(input_shape, p[0].shape)

print('Building RCAN model')
print('  - input_shape =', input_shape)
print('  - num_channels =', num_channels)
print('  - num_residual_blocks =', num_residual_blocks)
print('  - num_residual_groups =', num_residual_groups)
print('  - channel_reduction =', channel_reduction)

model = build_rcan(
    (*input_shape, 1),
    num_channels=num_channels,
    num_residual_blocks=num_residual_blocks,
    num_residual_groups=num_residual_groups,
    channel_reduction=channel_reduction)

model.compile(
    optimizer=keras.optimizers.Adam(lr=1e-4),
    loss=mae,
    metrics=[psnr])

data_gen = DataGenerator(
    input_shape,
    1,
    intensity_threshold=0, #0.25
    area_ratio_threshold=0, #0.05
    transform_function=None)

training_data = data_gen.flow(*list(zip(*training_data)))

if validation_data is not None:
    validation_data = data_gen.flow(*list(zip(*validation_data)))
    checkpoint_filepath = 'weights_{epoch:03d}_{val_loss:.8f}.hdf5'
else:
    checkpoint_filepath = 'weights_{epoch:03d}_{loss:.8f}.hdf5'

#output_dir = pathlib.Path(args.output_dir)
#output_dir.mkdir(parents=True, exist_ok=True)

print('Training RCAN model')
model.fit_generator(
    training_data,
    epochs= epochs, #config['epochs'],
    steps_per_epoch= steps_per_epoch, #config['steps_per_epoch'],
    validation_data=validation_data,
    validation_steps=steps_per_epoch, #config['steps_per_epoch'],
    verbose=0,
    callbacks=[
        keras.callbacks.LearningRateScheduler(
            staircase_exponential_decay(epochs // 4)),
        keras.callbacks.ModelCheckpoint(
           # str(output_dir / checkpoint_filepath),
            output_dir + '\\' + checkpoint_filepath,
            monitor='loss' if validation_data is None else 'val_loss',
            save_best_only=True),
        keras.callbacks.TensorBoard(
            log_dir=str(output_dir),
            write_graph=False),
        tqdm.keras.TqdmCallback(
            tqdm_class=functools.partial(tqdm.tqdm, dynamic_ncols=True))
    ])
