# -*- coding: utf-8 -*-
"""
Generate crops around soma and contrast-enhanced images for
1) assistance for human to assign cell type labels
2) ML classifier to make prediction
"""

import numpy as np
from skimage import io, filters, measure
from skimage.feature import hog
import tifffile
import imagecodecs
from skimage.io import imsave
import imageio
import cv2
from pathlib import Path
import os
import shutil
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PIL import Image
from PIL import ImageDraw
import seaborn as sns
import sys
from tqdm import tqdm
import pandas as pd
from collections import Counter

def viz_img(tif_np):
    """visualize image"""
    plt.imshow(tif_np, cmap='gray')
    plt.show()
    
def crop_around_point(image, center, crop_size):
    """ generate crops around soma"""
    # top and left corner
    top_left = (center[0] - crop_size // 2, center[1] - crop_size // 2)
    # bottom and right corner
    bottom_right = (center[0] + crop_size // 2, center[1] + crop_size // 2)
    # crop the image
    cropped = image[top_left[1]:bottom_right[1], top_left[0]:bottom_right[0],:]
    
    return cropped

def preprocess_img(img, soma_x, soma_y, r):
    """ perform cropping, histogram normalization, min-max scaling of image"""
    # crop an area around the soma
    ch1_tif_cropped = crop_around_point(img, (soma_y, soma_x), r)
   
    # histogram normalization
    histogram, _ = np.histogram(ch1_tif_cropped.flatten(), bins=256, range=(0, 256))
    cdf = histogram.cumsum()
    cdf_normalized = (cdf - cdf.min()) * 255 / (cdf.max() - cdf.min())
    ch1_cropped_normalized = cdf_normalized[ch1_tif_cropped]
    ch1_cropped_normalized = ch1_cropped_normalized.astype(np.uint8)
    
    ch1_cropped_normalized_r = np.copy(ch1_cropped_normalized)
    ch1_cropped_normalized_r[ch1_cropped_normalized_r < 200] = 0
    
    return ch1_cropped_normalized


####################################################################################################
################## generate enhanced crops around soma (for each brain at a time) ##################
####################################################################################################
# specify the brain
brain = '...' 
section_list = [rf'0{s}' for s in [i  for i in range(1, 10)]] # define slices

# define output directory
dest_dir_base = Path(r'...')
dest_dir = dest_dir_base/rf'{brain}'
dest_dir.mkdir(exist_ok=True)

# define the directory of morphometric data
morpho_dir = Path(r'...') 
morpho_df = pd.read_csv(morpho_dir/'....csv')
morpho_neuron_list = list(morpho_df['file_path'])
print(len(morpho_neuron_list))

failed = []
neuron_names = []
for section in section_list:
    
    # define the directory that contains ch1 images
    dir = Path(rf'...{section}...')

    r_around_soma = 500

    for folder in tqdm([n for n in os.listdir(dir) if not n.startswith('x_')]):
        # print(folder)
        
        # restrict to neurons that are in the morphometric dataset
        neuron_name  = '_'.join(folder.split('_')[1:])
        if neuron_name in morpho_neuron_list:
        
            folder_path = dir/folder/'ch1'
            try:
                ch1_tif_file = list(folder_path.rglob('*.tif'))[0]
                
                ch1_tif_file_name = ch1_tif_file.name.replace('.tif', '')
                # get the soma location in local space
                x, y, z, radii = 1024, 1024, 150, 25
                soma_coord = (1024, 1024, 150)
        
                # only load particular z ranges (+- n slices around the soma)
                with tifffile.TiffFile(ch1_tif_file) as tif:
                    ch1_tif = tif.asarray(key = slice(z-3, z+3+1)) # depth, height, width
                    # swap axies
                    ch1_tif = np.swapaxes(ch1_tif, 0, 2) # width (x), height (y), depth (z)
                    
                # process image (without min-max scaling to keep between 0~255)
                ch1_tif_cropped = preprocess_img(ch1_tif, soma_coord[0], soma_coord[1], r_around_soma)
                ch1_tif_cropped_to_save = ch1_tif_cropped[:,:,3].squeeze()
                # draw a circle around the soma points
                pil_img = Image.fromarray(ch1_tif_cropped_to_save)
                pil_img = pil_img.convert('RGB')
                
                draw = ImageDraw.Draw(pil_img)
                draw.ellipse([(r_around_soma/2 - radii, r_around_soma/2 - radii), (r_around_soma/2 + radii, r_around_soma/2 + radii)], outline='yellow', width=1)
                pil_img.save(dest_dir/rf"{folder}_processed.png", 'PNG')
                
                neuron_names.append(neuron_name)
                
            except Exception as e:
                print(e)
                print(folder)
                failed.append(folder)

# save a .csv file with list of neurons for human to assign gold-standard cell type labels
dest_df = pd.DataFrame(neuron_names, columns=['file_path'])   
dest_df.to_csv(dest_dir/rf'neuron_list_{brain}.csv', index=None)
