import difussionFilter
from skimage import exposure, io
import numpy as np
from skimage.morphology import opening, disk, square
from skimage.morphology import closing
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops
import cv2 as cv

def cropImage(img):
    if img.shape[2]==3:
        img = cv.cvtColor(img, cv.COLOR_BGR2GRAY)
    else:
        img  = img
    # Assuming 'img' is the grayscale image loaded earlier with OpenCV
    DIM = img.shape  # Get the dimensions of the image
    # Calculate crop dimensions
    X = int(DIM[1] * 0.1)  # 10% of the columns on the right
    Y = int(DIM[0] * 0.15)  # 15% of the rows from the top
    X2 = int(DIM[1] - DIM[1] * 0.1)  # Remaining 90% of the columns, to crop 10% from the left

    # Crop the image
    cropped_img = img[Y-1:DIM[0], X:X2]
    return cropped_img

def ApplyFilter(img, niter=20, kappa=30, gamma=0.25, option=2):
    return difussionFilter.anisodiff(img,niter=niter,kappa=kappa,gamma=gamma,option=option)

def strechedHistogram(img):
    p2, p98 = np.percentile(img, (2, 98))
    img_rescaled = exposure.rescale_intensity(img, in_range=(p2, p98))
    return img_rescaled

def morphologicalOpening(img, strelSize, type):
    if type=='disk':
        selem = disk(strelSize)
    elif type=='square':
        selem = square(strelSize)
    img_open = opening(img, selem)
    return img_open

def morphologicalClosing(img, strelSize):
    selem = disk(strelSize)
    imsubclos = closing(img, selem)
    return imsubclos

def ApplyOtsu(img):
    dim = img.shape
    X = int(dim[0] * 0.5)
    Y = int(dim[1] * 0.5)
    crop1 = img[:X, :Y]        # First section
    crop2 = img[:X, Y:]        # Second section
    crop3 = img[X:, :Y]        # Third section
    crop4 = img[X:, Y:]        # Fourth section

    # Calculate the Otsu threshold for the whole image and each section
    T = threshold_otsu(img)
    T1 = threshold_otsu(crop1)
    T2 = threshold_otsu(crop2)
    T3 = threshold_otsu(crop3)
    T4 = threshold_otsu(crop4)

    # Collect all thresholds and select the maximum
    thresholds = [T, T1, T2, T3, T4]
    max_threshold = max(thresholds)
    binary_image = img > max_threshold
    return binary_image




