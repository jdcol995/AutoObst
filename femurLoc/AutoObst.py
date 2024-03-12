import pydicom
from matplotlib import pyplot as plt
import numpy as np
from skimage.measure import label, regionprops
from Utils import cropImage, ApplyFilter, ApplyOtsu, morphologicalClosing, morphologicalOpening, strechedHistogram
import math
from matplotlib.pyplot import hist
from skimage.morphology import skeletonize

class AutoFemurLoc:
    """
    AutoFemurLoc is a class that, given a 2D ultrasound (US) image where the fetal femur is located,
    identifies the femur and mark the endpoints. If the metadata of the image includes the spatial resolution,
    it provides a measurement in millimeters. This class is based on the method presented in 
    "Fully automatic segmentation and measurement of the fetal femur," published in the Proceedings 
    Volume 10975 of the 14th International Symposium on Medical Information Processing and Analysis. 
    DOI: 10.1117/12.2511534.
    """

    def __init__(self, strel_1 = 25, strel_2 = 4, strel_3 = 9 ):
        """
        The class constructors need three parameters, the size of the structuring element for the
        three diferent morphological operators, this method is based on morphological operator to identify
        the bone region of the femur
        """
        self.strel_1 = strel_1
        self.strel_2 = strel_2
        self.strel_3 = strel_3


    def findFemur(self, img):
        """
        The function findFemur is the collarbone of this class, it has all the steps described in the 
        article, it returns three diferent set of points to find the extrema. the first one is
        the extrema of the bounding box of the region identified as femur; the second is the end points
        of the skeletonized region; and the final one is a compensated points that add the differences 
        between the original region and the skeletonized along the x-axis.
        """

        # Automatic crop of the original image
        croped_image = cropImage(img)
        # Anisotropic difussion filtering
        filtered_image = ApplyFilter(croped_image,20,30,0.25,2)
        # Histogram streching for contrast enhancing
        streched_image = strechedHistogram(filtered_image)
        # grayscale opening and subsequently subtraction
        open_grayscaled_image = morphologicalOpening(streched_image,self.strel_1, 'disk')
        subtracted_image = streched_image - open_grayscaled_image
        # Grayscale closing and threshold of the image
        closed_grayscaled_image = morphologicalClosing(subtracted_image,self.strel_2)
        thresholed_image = ApplyOtsu(closed_grayscaled_image)
        # Final opening on the binary image
        open_thresholed_image = morphologicalOpening(thresholed_image,self.strel_3,'square')
        # Labeling of the remaning regions
        labels = label(open_thresholed_image)

        # Removing smaller regions
        dim = open_thresholed_image.shape
        max_L = math.sqrt((dim[0]-1)**2 + (dim[1]-1)**2)

        lengths = []
        X1, X2, Y1, Y2 = [], [], [], []
        small_regions = []
        candidate_regions = []

        for region in regionprops(labels):
            minr, minc, maxr, maxc = region.bbox
            X1.append(minr)
            Y1.append(minc)
            X2.append(maxr)
            Y2.append(maxc)
            dist = np.sqrt((minr - maxr) ** 2 + (minc - maxc) ** 2)
            
            if dist > (max_L * 0.1):
                lengths.append(dist)
                candidate_regions.append(region.label)
            else:
                small_regions.append(region.label)

        for label_i in small_regions:
            labels[labels == label_i] = 0
        props = regionprops(labels, intensity_image=streched_image)    
        # Finding the Intensity mean and position of the centroid of each region
        I_region = [prop.intensity_mean for prop in props]
        centroid_region = [prop.centroid for prop in props]
        # Calculating and saving the entropy of each regionn as an indicator for texture
        entropies = []
        texture = []
        for region_label in candidate_regions:
            mask = labels == region_label
            region_image = np.where(mask, streched_image, 0)
            region_image = region_image / np.max(region_image)
            histogram, _ = np.histogram(region_image[region_image > 0], bins=256, range=(1, np.max(region_image)))
            p = histogram / np.sum(histogram)
            p = p[p > 0]
            H = -np.sum(p * np.log2(p))
            entropies.append(H)
        texture = [1 / H for H in entropies]
        # Calculating scores
        centroid_region = np.array(centroid_region)
        scores = []
        for k in range(len(lengths)):
            score = (1/4) * ((I_region[k] / max(I_region)) +
                            (lengths[k] / max(lengths)) +
                            (centroid_region[k, 0] / max(centroid_region[:, 0])) +
                            (texture[k] / max(texture)))
            scores.append(score)
        # Identify the region with the highest score
        femur_loc_index = np.argmax(scores)  # Index of the femur location in the candidate_regions list
        femur_loc = candidate_regions[femur_loc_index]  # The label of the femur location

        # Merging similar regions
        # Intensity limit
        intensity_limit = I_region[femur_loc_index] - 10
        # Centroid position limit
        centroid_limit = np.abs(centroid_region[:, 0] - centroid_region[femur_loc_index, 0])

        # Find regions with similar characteristics
        regions_i = [i for i in range(len(I_region)) if I_region[i] > intensity_limit]
        centroids_i = [i for i in range(len(centroid_region)) if centroid_limit[i] < 10]

        # If there are regions that meet both criteria, merge them
        similar_regions = set(regions_i).intersection(set(centroids_i))
        for region_index in similar_regions:
            region_label = candidate_regions[region_index]
            labels[labels == region_label] = femur_loc  # Merge by setting their labels to the femur location label

        mask = labels == femur_loc
        # Find the coordinates of the mask
        y, x = np.where(mask)
        # Find extrema points
        min_x, max_x = np.min(x), np.max(x)
        min_y, max_y = np.min(y), np.max(y)

        original_P1 = min_x, y[np.argmin(x)]
        original_P2 = max_x, y[np.argmax(x)]

        # Thining of the mask image
        skeleton = skeletonize(mask)
        # Find coordinates of the thinned (skeletonized) region
        r, c = np.where(skeleton)
        # Original extremal points
        uncompensated_P2 =  c[0], r[0]
        uncompensated_P1 = c[-1],r[-1]
        compP1, compP2 = abs(uncompensated_P1[0] - original_P1[0]), abs(uncompensated_P2[0] - original_P2[0])
        compensated_P1 = (uncompensated_P1[0]-compP1,uncompensated_P1[1])
        compensated_P2 = (uncompensated_P2[0]+compP2,uncompensated_P2[1])

        # return the three set of poits as an array, first set original bounding box end points,
        # second, uncompensated thined end ponits,
        # third, compensated thined end points
        return np.array([[original_P1,original_P2],[uncompensated_P1,uncompensated_P2],[compensated_P1,compensated_P2]])

    def set_markers(self,img,end_points):
        """ 
        This function marks over the original cropped image, it only takes
        a 2x2 array, it has to be given only a pair of points
        """
        croped_image = cropImage(img)
        plt.imshow(croped_image, cmap='gray')
        plt.scatter(end_points[0][0], end_points[0][1], color='red')
        plt.scatter(end_points[1][0], end_points[1][1], color='red')

    def meassure_length(self, pixel_spacing, end_points):
        """
        This function takes the spatial spacing an the endpoints to perform the meassure in milimeters,
        returns this value
        """
        length = np.sqrt((end_points[0,0] - end_points[1,0]) ** 2 + (end_points[0,1] - end_points[1,1]) ** 2)
        return length * pixel_spacing











