# -*- coding: utf-8 -*-
"""
Train 3 models (cnn, svm, rf) using human human-created gold standard dataset for D1- and D2-MSN classification based on channel images
final result is based on weighted pooled probabilities
The models will output probabilities instead of hard labels
images: dragonfly images 30x, D1-td tomatal labeling channel
"""

import numpy as np
from skimage import io, filters, measure
from skimage.feature import hog
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
from sklearn.utils import class_weight
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, accuracy_score
import tensorflow
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten, Conv3D, MaxPooling3D, Dropout, BatchNormalization, Conv2D, MaxPooling2D
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from tensorflow.keras.callbacks import EarlyStopping
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn import metrics
import pickle
print(sys.getrecursionlimit())

# define work directory
os.chdir(Path(r'...'))
print(os.getcwd())
from augmented import generator

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
    cropped = image[top_left[1]:bottom_right[1], top_left[0]:bottom_right[0],:]
    
    return cropped

def preprocess_img(img, soma_x, soma_y, r):
    """ perform cropping, histogram normalization, threshold-based segmentation, min-max scaling of image"""
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
    ch1_cropped_normalized_r = (ch1_cropped_normalized_r - np.min(ch1_cropped_normalized_r)) / (np.max(ch1_cropped_normalized_r) - np.min(ch1_cropped_normalized_r))
    
    return ch1_cropped_normalized_r

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

def posthoc(prediction, labels, model_name):
    
    """plot confusion matrix, grab incorrect classification to a different folder, output a csv summary file
       prediction: list of predicted labels
       labels: list of true labels
    """
    # confusion matrix and visualization
    cm = confusion_matrix(labels, prediction)

    ax= plt.subplot()
    sns.heatmap(cm, annot=True, fmt='g', ax=ax);  #annot=True to annotate cells, ftm='g' to disable scientific notation
    # labels, title and ticks
    ax.set_xlabel('Predicted labels')
    ax.set_ylabel('True labels'); 
    ax.set_title(rf'Confusion Matrix of D1D2 classification {model_name}'); 
    ax.xaxis.set_ticklabels(['D1', 'D2'])
    ax.yaxis.set_ticklabels(['D1', 'D2']);

    # calculate metrix
    accuracy = accuracy_score(labels, prediction)
    print(f'Accuracy: {accuracy}')

    tp, fp, tn, fn = cm[0,0], cm[1,0], cm[1,1], cm[0,1]
    sensitivity = tp/(tp+fp)
    specificity = tn/(tn+fn)
    print(f'Sensitivity: {sensitivity}\nSpecificity: {specificity}')
    

def morph_ops(mask, kernel_size, open_iter, dilate_iter):
    """apply morphological operations to further remove noise on the thresholded mask image"""
    kernel = np.ones (kernel_size, np.uint8)
    opening = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel, iterations=open_iter)
    sure_bg = cv2.dilate(opening, kernel, iterations=dilate_iter)
    return sure_bg

def plot_learning_curve(history_obj):
    # Extract loss and accuracy values
    train_loss = history_obj.history['loss']
    val_loss = history_obj.history['val_loss']
    train_accuracy = history_obj.history['accuracy']
    val_accuracy = history_obj.history['val_accuracy']
    
    # Create subplots for loss and accuracy
    plt.figure(figsize=(12, 5))
    
    # Plot training and validation loss
    plt.subplot(1, 2, 1)
    plt.plot(train_loss, label='Training Loss')
    plt.plot(val_loss, label='Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Training and Validation Loss')
    plt.legend()
    
    # Plot training and validation accuracy
    plt.subplot(1, 2, 2)
    plt.plot(train_accuracy, label='Training Accuracy')
    plt.plot(val_accuracy, label='Validation Accuracy')
    plt.xlabel('Epoch')
    plt.ylabel('Accuracy')
    plt.title('Training and Validation Accuracy')
    plt.legend()
    
    plt.tight_layout()
    plt.show()

#######################################################################################
####################### train a model based on gold-standard human assignment##############
######################################################################################
train_data_dir = Path(r'...')

# generate train images
neuron_list = []
img_list_svm = []
img_list_rf = []
img_list_cnn = []

r_around_soma_svm = 60
r_around_soma_rf = 70
r_around_soma_cnn = 70

for folder in tqdm([f for f in os.listdir(train_data_dir) if not f.startswith('x')]):
    
    ch1_tif_file = [f for f in os.listdir(train_data_dir/folder) if '.tif' in f][0]
    # get the soma location in local space
    QCd_swc = [f for f in os.listdir(train_data_dir/folder) if '_Final.swc' in f][0]
    neuron_name = folder
    with open(train_data_dir/folder/QCd_swc) as f:
        l = [l for l in f.readlines() if '#' not in l][0]
        l_split = l.split(' ')
        x, y, z, radii = int(float(l_split[2])),  int(float(l_split[3])),  int(float(l_split[4])), float(l_split[5])
        soma_coord = (x,y,z, radii)
    f.close()
        
    # only load particular z ranges (+- n slices around the soma)
    with tifffile.TiffFile(train_data_dir/folder/ch1_tif_file) as tif:
        ch1_tif = tif.asarray(key = slice(z-3, z+3+1))
        # swap axies
        ch1_tif = np.swapaxes(ch1_tif, 0, 2) # height, width, depth
    
    neuron_list.append(neuron_name)
    # process image for svm, blurring + specific size
    ch1_tif_cropped_svm = preprocess_img(ch1_tif, soma_coord[1], soma_coord[0], r_around_soma_svm)
    ch1_tif_cropped_svm_blurred = gaussian_filter(input = ch1_tif_cropped_svm, sigma = (3,3,3), mode='wrap')
    img_list_svm.append(ch1_tif_cropped_svm_blurred)
    
    # process image for rf, blurring + specific size
    ch1_tif_cropped_rf = preprocess_img(ch1_tif, soma_coord[1], soma_coord[0], r_around_soma_rf)
    ch1_tif_cropped_rf_blurred = gaussian_filter(input = ch1_tif_cropped_rf, sigma = (3,3,3), mode='wrap')
    img_list_rf.append(ch1_tif_cropped_rf_blurred)
    
    # process image for cnn, no blurring + specific size
    ch1_tif_cropped_cnn = preprocess_img(ch1_tif, soma_coord[1], soma_coord[0], r_around_soma_cnn)
    img_list_cnn.append(ch1_tif_cropped_cnn)
    
      
print(len(neuron_list))
print(len(img_list_svm))
print(len(img_list_rf))
print(len(img_list_cnn))

img_array_svm = np.array(img_list_svm)
img_array_rf = np.array(img_list_rf)
img_array_cnn = np.array(img_list_cnn)

# Read human-labeld cell type and create the labels for training dataset
CP_labeled_dir = Path(r'...')
CP_labeled = pd.read_csv('...')

label_list = []
for i, neuron in enumerate(neuron_list):
    try:    
        label = CP_labeled.loc[CP_labeled['file_name'] == neuron, 'CP_second_round'].values.item()
        label_list.append(label)
    except:
        print(neuron)

label_list = np.array(label_list)
label_list_numeric = LabelEncoder().fit_transform(label_list) # D1: 0, D2: 1

class_weights = class_weight.compute_class_weight('balanced', classes=np.unique(label_list_numeric), y=label_list_numeric)
class_weights_dict = {0: class_weights[0], 1: class_weights[1]}
print(class_weights_dict)

####################################################################################
########################### 1. Training SVM model ######################################
# only keep the middle slice for simplicity
img_array_svm_flat = img_array_svm[:,:,:, 3].reshape(img_array_svm.shape[0], -1)
img_array_svm_flat.shape
x_train_svm_flat, x_test_svm_flat, y_train, y_test = train_test_split(img_array_svm_flat, label_list_numeric, test_size = 0.3, random_state = 42)
print(x_train_svm_flat.shape, y_train.shape)
print(x_test_svm_flat.shape, y_test.shape)

svm_clf = svm.SVC(max_iter = -1, verbose = True, class_weight = class_weights_dict, probability = True)
params_grid = {'C': [1, 10, 100], 
               'gamma': [0.1, 0.01, 0.001,'scale'],
               'kernel': ['rbf', 'linear', 'poly']}
svm_grid_search = GridSearchCV(svm_clf, params_grid, scoring='accuracy', cv=3)
svm_grid_search.fit(x_train_svm_flat, y_train)

print('best params: ', svm_grid_search.best_params_)
print('best score: ', svm_grid_search.best_score_)

best_svm = svm_grid_search.best_estimator_ 

print('test accuracy: ', metrics.accuracy_score(y_test, best_svm.predict(x_test_svm_flat))) 

posthoc(best_svm.predict(x_test_svm_flat), y_test, 'SVM')

# save model to defined directory
with open(Path(r'...\trained_svm_Oct-26-2023.pkl'), 'wb') as f:
    pickle.dump(best_svm,f)

# find out the correctly vs. incorrectly classified samples
test_indices = []
for t in x_test_svm_flat:
    ind = np.where(np.all(img_array_svm_flat == t, axis=1))[0][0]
    test_indices.append(ind)
len(test_indices)

# correctly classified vs. incorrectly classified
svm_correct = []
svm_correct_true_label = []
svm_incorrect = []
svm_incorrect_true_label = []
svm_incorrect_neuron = []

for i,t in enumerate(x_test_svm_flat):
    reshaped_input = x_test_svm_flat[i].reshape(1,-1)
    if y_test[i] == best_svm.predict(reshaped_input):
        svm_correct.append(img_array_svm[test_indices[i]])
        svm_correct_true_label.append(label_list[test_indices[i]])
    elif y_test[i] != best_svm.predict(reshaped_input):
        svm_incorrect.append(img_array_svm[test_indices[i]])
        svm_incorrect_true_label.append(label_list[test_indices[i]])
        svm_incorrect_neuron.append(neuron_list[test_indices[i]])


i=1
viz_img(svm_correct[i][:,:,3])
print(rf"correctly classified as {svm_correct_true_label[i]}")

i=4
viz_img(svm_incorrect[i][:,:,3])
print(rf"{svm_correct_true_label[i]} cell incorrectly classified")

####################################################################################
########################### 2. Training RF model ######################################
img_array_rf_flat = img_array_rf[:,:,:, 3].reshape(img_array_rf.shape[0], -1)
img_array_rf_flat.shape
x_train_rf_flat, x_test_rf_flat, y_train, y_test = train_test_split(img_array_rf_flat, label_list_numeric, test_size = 0.3, random_state = 42)
print(x_train_rf_flat.shape, y_train.shape)
print(x_test_rf_flat.shape, y_test.shape)

rf_clf = RandomForestClassifier(verbose = True, class_weight = class_weights_dict)
params_grid = {'n_estimators': [100,200,300],
               'max_depth': [None, 10, 20,50],
               'min_samples_split': [2,5,10]}
rf_grid_search = GridSearchCV(rf_clf, params_grid, scoring='accuracy', cv=3)
rf_grid_search.fit(x_train_rf_flat, y_train)

print('best params: ', rf_grid_search.best_params_)
print('best score: ', rf_grid_search.best_score_) 

best_rf = rf_grid_search.best_estimator_ 

print('test accuracy: ', metrics.accuracy_score(y_test, best_rf.predict(x_test_rf_flat))) 

posthoc(best_rf.predict(x_test_rf_flat), y_test, 'Random Forest')

# correctly classified vs. incorrectly classified
rf_correct = []
rf_correct_true_label = []
rf_incorrect = []
rf_incorrect_true_label = []
rf_incorrect_neuron = []
for i,t in enumerate(x_test_rf_flat):
    reshaped_input = x_test_rf_flat[i].reshape(1,-1)
    if y_test[i] == best_rf.predict(reshaped_input):
        rf_correct.append(img_array_rf[test_indices[i]])
        rf_correct_true_label.append(label_list[test_indices[i]])
    elif y_test[i] != best_rf.predict(reshaped_input):
        rf_incorrect.append(img_array_rf[test_indices[i]])
        rf_incorrect_true_label.append(label_list[test_indices[i]])
        rf_incorrect_neuron.append(neuron_list[test_indices[i]])

i=8
viz_img(rf_correct[i][:,:,3])
print(rf"correctly classified as {rf_incorrect_true_label[i]}")

i=7
viz_img(rf_incorrect[i][:,:,3])
print(rf"{rf_incorrect_true_label[i]} cell incorrectly classified")

# save model to defined directory
with open(Path(r'...\trained_rf_Oct-25-2023.pkl'), 'wb') as f:
    pickle.dump(best_rf,f)

rf_incorrect_neuron
svm_incorrect_neuron
list(set(rf_incorrect_neuron) - set(svm_incorrect_neuron))

####################################################################################
######################## Training 3. 3D CNN model ######################################
# split into train and validation
x_train_cnn, x_test_cnn, y_train_cnn, y_test_cnn = train_test_split(img_array_cnn, label_list_numeric, test_size = 0.2, random_state = 42)
print(x_train_cnn.shape, y_train_cnn.shape) # 254 in training
print(x_test_cnn.shape, y_test_cnn.shape) # 64 in testing

# add data augmentation
x_train_cnn, x_val_cnn, y_train_cnn, y_val_cnn = train_test_split(x_train_cnn, y_train_cnn, test_size=0.2, random_state=42)
print(x_train_cnn.shape, y_train_cnn.shape) # 203 in training
print(x_val_cnn.shape, y_val_cnn.shape) # 51 in val

BATCH_SIZE = 10

train_aug = generator.customImageDataGenerator(rotation_range=90, 
                                               # width_shift_range=0.2,
                                               # height_shift_range=0.2,
                                               # zoom_range=0.2,
                                               horizontal_flip=True,
                                               vertical_flip=True)

train_datagen = train_aug.flow(np.expand_dims(x_train_cnn, axis=-1), y_train_cnn, batch_size=BATCH_SIZE, seed=42)
val_datagen = train_aug.flow(np.expand_dims(x_val_cnn, axis=-1), y_val_cnn, batch_size=BATCH_SIZE, seed=42)

model_cnn = Sequential()
model_cnn.add(Conv3D(32, kernel_size=(3, 3, 1), activation='relu', input_shape=(x_train_cnn.shape[1], x_train_cnn.shape[2], x_train_cnn.shape[3], 1)))
model_cnn.add(MaxPooling3D(pool_size=(2, 2, 1)))
model_cnn.add(Conv3D(64, kernel_size=(5, 5, 1), activation='relu', padding='same'))
model_cnn.add(MaxPooling3D(pool_size=(2, 2, 1)))
model_cnn.add(Conv3D(128, kernel_size=(7, 7, 1), activation='relu', padding='same'))
model_cnn.add(MaxPooling3D(pool_size=(2, 2, 1)))
model_cnn.add(Flatten())
model_cnn.add(Dense(128, activation='relu'))
model_cnn.add(Dense(1, activation='sigmoid'))

model_cnn.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model_cnn.summary()

early_stopping = EarlyStopping(monitor='val_loss', patience=7, restore_best_weights=True)

# history_cnn = model_cnn.fit(x_train, y_train, validation_split = 0.2, epochs=100, batch_size=BATCH_SIZE, callbacks=[early_stopping])
history_cnn = model_cnn.fit(train_datagen, epochs=100, validation_data = val_datagen, 
                            steps_per_epoch = len(y_train_cnn) // BATCH_SIZE,
                            validation_steps = len(y_val_cnn) // BATCH_SIZE,
                            callbacks = [early_stopping], verbose=1, 
                            class_weight = class_weights_dict)

predicted = (model_cnn.predict(x_test_cnn) > 0.5).astype(int)
print('test accuracy: ', accuracy_score(y_test_cnn, predicted)) # 90.63 %

plot_learning_curve(history_cnn)
posthoc(predicted, y_test_cnn, 'CNN')

# save model to the defined directory
model_cnn.save(Path(r'...\trained_cnn_Oct-25-2023.h5'))
