import cv2
import numpy as np

import sharedmatting


if __name__ == '__main__':
    im = cv2.imread('input.png')
    trimap = cv2.imread('trimap.png', cv2.IMREAD_GRAYSCALE)
    result = np.zeros(trimap.shape, dtype='uint8')
    sharedmatting.solve_alpha(im, trimap, result)
    cv2.imshow('a.png', result)
    cv2.waitKey()
