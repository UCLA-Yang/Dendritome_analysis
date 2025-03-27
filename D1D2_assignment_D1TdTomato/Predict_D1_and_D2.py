# -*- coding: utf-8 -*-
"""
load the trained models, and use a voting scheme to decide the label (D1/D2)
for the output of each model, use probability instead of hard labels
for the voted results, use the weighted sum of probabilities for final prediction
"""

import numpy as np
from skimage import io, filters, measure
from skimage.feature import hog
from PIL import Image
import tifffile
import imagecodecs
import cv2
from scipy.ndimage import gaussian_filter
from pathlib import Path
import os
import shutil
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from tqdm import tqdm
import pandas as pd
from collections import Counter
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, accuracy_score
import tensorflow
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten, Conv3D, MaxPooling3D
from tensorflow.keras.models import load_model
from tensorflow.keras.callbacks import EarlyStopping
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn import metrics
import pickle
print(sys.getrecursionlimit())

# load model SVM
with open(Path(r'...\trained_svm_Oct-26-2023.pkl'), 'rb') as f:
    trained_svm_clf = pickle.load(f)
# load model RF
with open(Path(r'...\trained_rf_Oct-26-2023.pkl'), 'rb') as f:
    trained_rf_clf = pickle.load(f)
trained_rf_clf.verbose = False
# load model CNN
trained_cnn = load_model(r'...\trained_cnn_Oct-26-2023.h5')

global trained_svm_clf
global trained_rf_clf
global trained_cnn

### Note: can apply more image proecssing steps before classification###
def viz_img(tif_np):
    """visualize image"""
    plt.imshow(tif_np, cmap='gray')
    plt.show()

def plot_histogram(tif_np):
    """histogram distribution of pixel intensities"""
    n, bins, patches = plt.hist(tif_np.flatten(), bins='auto', color='#0504aa', alpha=0.7)
    plt.xlabel('Normalized pixel intensity')
    plt.ylabel('Frequency')
    plt.title('Distribution of normalized pixel intensity')
    plt.show()

def crop_around_point(image, center, crop_size):
    """ generate crops around soma"""
    # top and left corner
    top_left = (center[0] - crop_size // 2, center[1] - crop_size // 2)
    # bottom and right corner
    bottom_right = (center[0] + crop_size // 2, center[1] + crop_size // 2)
    # crop the image
    cropped = image[top_left[1]:bottom_right[1], top_left[0]:bottom_right[0]]
    
    return cropped

def preprocess_img(img, soma_x, soma_y, r):
    """ perform cropping, histogram normalization, thresholding-based segmentation, min-max scaling of image"""
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
    
    # min-max scaling
    ch1_cropped_normalized_r_scaled = (ch1_cropped_normalized_r - np.min(ch1_cropped_normalized_r)) / (np.max(ch1_cropped_normalized_r) - np.min(ch1_cropped_normalized_r))
    
    return ch1_cropped_normalized_r_scaled

def find_neighbors(mask, start_point):
    """find connected components of a given point"""
    
    stack = [start_point]
    visited = set()

    while stack:
        x, y = stack.pop()
        visited.add((x, y))

        # Check the neighboring cells
        neighbors = [(x-1, y), (x+1, y), (x, y-1), (x, y+1)]
        for nx, ny in neighbors:
            if 0 <= nx < len(mask) and 0 <= ny < len(mask[0]):
                if mask[nx][ny] == 1 and (nx, ny) not in visited:
                    stack.append((nx, ny))
                    visited.add((nx, ny))
    return visited

def morph_ops(mask, kernel_size, open_iter, dilate_iter):
    """apply morphological operations to further remove noise on the thresholded mask image"""
    kernel = np.ones (kernel_size, np.uint8)
    opening = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel, iterations=open_iter)
    sure_bg = cv2.dilate(opening, kernel, iterations=dilate_iter)
    return sure_bg

# function to scan throuhg enhanced image crops to classify D1 and D2
def batch_scan_tif_and_predict(brain_name, enhanced_crop_dir, r_around_soma_svm, r_around_soma_rf, r_around_soma_cnn):
    
    """
    enhanced_crop_dir: directory to the enhanced crops 
    r_around_soma: r in #of pixels to keep when cropping out the soma region
    
    returns: a dictionary with key=neuron name and value=voted predicted D1D2 labels
    """

    neuron_prediction = {}

    for neuron in [f for f in os.listdir(enhanced_crop_dir) if '.png' in f]:
        print(neuron)
      
        neuron_name = neuron.replace('.png', '')
        
        neuron_img = np.array(Image.open(enhanced_crop_dir/neuron))
        neuron_img = cv2.cvtColor(neuron_img, cv2.COLOR_BGR2GRAY)
        # print(neuron_img.shape)
        
        # process image
        # 1. svm image
        ch1_tif_cropped_svm = (neuron_img - np.min(neuron_img)) / (np.max(neuron_img) - np.min(neuron_img))
        ch1_tif_cropped_svm =  gaussian_filter(input = ch1_tif_cropped_svm, sigma = (3,3), mode='wrap')
        ch1_tif_cropped_svm = crop_around_point(ch1_tif_cropped_svm, (250, 250), r_around_soma_svm)
        # viz_img(ch1_tif_cropped_svm)

        # 2. rf image
        ch1_tif_cropped_rf = (neuron_img - np.min(neuron_img)) / (np.max(neuron_img) - np.min(neuron_img))
        ch1_tif_cropped_rf =  gaussian_filter(input = ch1_tif_cropped_rf, sigma = (3,3), mode='wrap')
        ch1_tif_cropped_rf = crop_around_point(ch1_tif_cropped_rf, (250, 250), r_around_soma_rf)
        
        # 3. cnn image
        ch1_tif_cropped_cnn = (neuron_img - np.min(neuron_img)) / (np.max(neuron_img) - np.min(neuron_img))
        ch1_tif_cropped_cnn = crop_around_point(ch1_tif_cropped_cnn, (250, 250), r_around_soma_cnn)
        ch1_tif_cropped_cnn = np.stack([ch1_tif_cropped_cnn] * 7, axis=-1)
        # ch1_tif_cropped_cnn.shape
        
        # make predcition
        neuron_prediction[neuron_name] = {}
        ch1_tif_cropped_svm_flatted = ch1_tif_cropped_svm.reshape(1, -1)
        ch1_tif_cropped_svm_flatted.shape
        
        ch1_tif_cropped_rf_flatted = ch1_tif_cropped_rf.reshape(1, -1)
        # 1) SVM model
        neuron_prediction[neuron_name]['SVM_Probability_D2'] = trained_svm_clf.predict_proba(ch1_tif_cropped_svm_flatted)[0][1]
        # 2) RF model
        neuron_prediction[neuron_name]['RF_Probability_D2'] = trained_rf_clf.predict_proba(ch1_tif_cropped_rf_flatted)[0][1]
        # 3) CNN model
        cnn_pred = trained_cnn.predict(ch1_tif_cropped_cnn.reshape(-1, ch1_tif_cropped_cnn.shape[0], ch1_tif_cropped_cnn.shape[1], ch1_tif_cropped_cnn.shape[2]))
        neuron_prediction[neuron_name]['CNN_Probability_D2'] = cnn_pred[0][0]
    
        # apply averaged probability to output the final probability for D2 neuron prediction
        sum_p = 0
        for values in neuron_prediction[neuron_name].values():
            sum_p += values
        neuron_prediction[neuron_name]['Average_Probability_D2'] = sum_p / len(neuron_prediction[neuron_name])
        neuron_prediction[neuron_name]['Average_Probability_D1'] = 1 - neuron_prediction[neuron_name]['Average_Probability_D2']
    
    # convert to dataframe and convert to human readable labels
    neuron_prediction_df = pd.DataFrame.from_dict(neuron_prediction, orient='index').reset_index()
    neuron_prediction_df.columns = ['file_path', 'SVM_Probability_D2', 'RF_Probability_D2', 'CNN_Probability_D2', 'Average_Probability_D2', 'Average_Probability_D1']
    neuron_prediction_df['Predicted_Label'] = np.where(neuron_prediction_df['Average_Probability_D2'] >= 0.5, 'D2', 'D1')
    neuron_prediction_df['difference'] = np.abs(neuron_prediction_df['Average_Probability_D2'] - neuron_prediction_df['Average_Probability_D1'])
    neuron_prediction_df['Certainty'] = np.where(neuron_prediction_df['difference'] <= 0.6, 'low', 'high')
    neuron_prediction_df = neuron_prediction_df[['file_path', 'SVM_Probability_D2', 'RF_Probability_D2', 'CNN_Probability_D2', 'Average_Probability_D1', 'Average_Probability_D2', 'Predicted_Label', 'Certainty']]
    
    neuron_prediction_df.to_csv(enhanced_crop_dir/rf'{brain}_model_D1D2_classification_1.csv', index=None)
    
    return neuron_prediction_df


#######################################################################################
############################# MAKE PREDICTIONS ########################################
#######################################################################################

r_around_soma_svm = 60
r_around_soma_rf = 70
r_around_soma_cnn = 70

failed_lst = []
global failed_lst

# define brain name of the neurons to predict
brain = '...'
# define the directory that has the enhanced 2D crop around the soma
enhanced_crop_dir = Path(r'...')

batch_scan_tif_and_predict(brain, enhanced_crop_dir, r_around_soma_svm, r_around_soma_rf, r_around_soma_cnn)    
    
    
    
    
    
    
