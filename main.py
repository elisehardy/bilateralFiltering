import cv2 as cv2
import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import math
from tkinter import *
size_window_gaussien = 5

#gausiien Ip = valeur de l'image I a la position p(x, y)
#Iq intensite pixel q

#F[I] = sortie du filtre F applique à l'image I
# q voisin du point p

#Y' = 0.2989 R + 0.5870 G + 0.1140 B


def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])


def gaussienKernel2D(val, sigma):
    return 1.0/(2*math.pi*(sigma**2)) * math.exp(-(val**2)/(2*(sigma**2)))


#def norme(a, b):
#    return math.sqrt((b[0]-a[0])**2 + (b[1]-a[1])**2)

def give_neighborhood(diametre_square, x, y, width, height):
    res = []
    for i in range(1, (diametre_square//2)):
        if x + i < width:
            res.append((x + i, y))
            if y + i < height:
                res.append((x + i, y + i))
            if y - i > height:
                res.append((x + i, y - 1))
        if x - i > width:
            res.append((x - i, y))
            if y - i > height:
                res.append((x - i, y - i))
            if y + i < height:
                res.append((x - 1, y + i))
        if y -i > height:
            res.append((x, y - i))
        if y + i < height:
            res.append((x, y + i))

    return res

def distance(x, y, i, j):
    return np.sqrt((x-i)**2 + (y-j)**2)


def bilateralFilter(image_source, image_filtred, signam_s, sigma_i):
    width = len(image_source)
    heigth = len(image_source[0])
    dim = 5
    for x in range(width):
        for y in range(heigth):
            wp = 0
            sum_gaus = 0
            #neig = give_neighborhood(dim, x,y, width, heigth)
            #for n in neig:
            i=0
            while i < 5:
                j = 0
                while j < 5:
                    neighbour_x = x - (5/2 - i)
                    neighbour_y = y - (5/2 - j)
                    if neighbour_x >= len(image_source):
                        neighbour_x -= len(image_source)
                    if neighbour_y >= len(image_source[0]):
                        neighbour_y -= len(image_source[0])
                    p = np.array((x, y))
                    q = np.array((neighbour_x, neighbour_y))
                    gaussien_pixel_s = gaussienKernel2D(np.linalg.norm(p-q), signam_s)
                    #gaussien_pixel_s = gaussienKernel2D(distance(neighbour_x,neighbour_y, x, y), signam_s)

                    gaussien_pixel_i = gaussienKernel2D(image_source[x][y]-image_source[abs(int(neighbour_x))][abs(int(neighbour_y))], sigma_i)
                    tmp = gaussien_pixel_s*gaussien_pixel_i
                    wp += tmp
                    sum_gaus += tmp*image_source[abs(int(neighbour_x))][abs(int(neighbour_y))]

                    j+=1
                i+=1
            image_filtred[x][y] = int(round(sum_gaus/wp))
    return image_filtred



def main():
    # image = mpimg.imread("chien2.png")
    # image = rgb2gray(image)
    # plt.imshow(image, cmap=plt.get_cmap('gray'), vmin=0, vmax=1)
    # plt.show()
    # image_filtred = np.zeros(image.shape)
    # image_filtred = bilateralFilter(image, image_filtred, 12.0, 16.0)
    # plt.imshow(image_filtred)
    # plt.imshow(image_filtred, cmap=plt.get_cmap('gray'), vmin=0, vmax=1)
    # plt.show()
    src = cv2.imread("chien2.png", 0)
    filtered_image_OpenCV = cv2.bilateralFilter(src, 5, 12.0, 16.0)
    cv2.imwrite("original_image_grayscale.png", src)
    cv2.imwrite("filtered_image_OpenCV.png", filtered_image_OpenCV)
    image_filtred = np.zeros(src.shape)

    filtered_image_own = bilateralFilter(src,image_filtred , 12.0, 16.0)
    cv2.imwrite("filtered_image_own.png", filtered_image_own)
    plt.imshow(filtered_image_own, cmap=plt.get_cmap('gray'), vmin=0, vmax=1)


if __name__ == "__main__":
    # execute only if run as a script
    main()