import re
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import bicgstab
import Reader as rd
from glob import glob
import matplotlib.patches as patches
import pandas as pd
import math
from scipy.misc import imread


import os

# Data cleaning & Annotation


def init_data(clear_path):
    def button_press(event):
        press.clear()
        press.append(int(event.ydata))
        press.append(int(event.xdata))
        return press

    def button_release(event):
        release.clear()
        release.append(int(event.ydata))
        release.append(int(event.xdata))
        return release

    def key_press_renew(event):
        temp = again
        if event.key == ' ':
            plt.close()
        if event.key == 'a':
            temp[0] = True
            plt.close()

    def key_press(event):
        if event.key == 'c':
            isRedExist[0] = True
        if event.key == ' ':
            plt.close()

    if os.path.exists(clear_path):
        print('Path exists!!!')
        return
    os.mkdir(clear_path)
    dataInfo = []
    index = 0
    labels = pd.read_csv("../stage1_labels.csv").T
    labels.columns = labels.iloc[0]
    label_dict = labels.to_dict()
    path_list = glob("./stage1_aps/*.aps")
    path_list.sort()
    i1 = 0
    i2 = 100
    i2 = 200
    i3 = 300
    i4 = 400
    i5 = 500
    i6 = 600
    i7 = 700
    i8 = 800
    i9 = 900
    i10 = 1000
    i11 = len(path_list)
    while index < len(path_list[0:2]):
        path = path_list[index]
        press = [40, 40]
        release = [100, 100]
        again = [False]
        isRedExist = [False]
        sample_id = re.search(r"(?<=s/)[^.]*", path).group(0)
        isThreatInBody = (label_dict[sample_id + '_Zone5']['Probability'] == 1 or
         label_dict[sample_id + '_Zone6']['Probability'] == 1 or
         label_dict[sample_id + '_Zone7']['Probability'] == 1 or
         label_dict[sample_id + '_Zone8']['Probability'] == 1 or
         label_dict[sample_id + '_Zone9']['Probability'] == 1 or
         label_dict[sample_id + '_Zone10']['Probability'] == 1)
        print(sample_id)
        if isThreatInBody:
            if label_dict[sample_id + '_Zone5']['Probability'] == 1:
                print('Zone5')
            if label_dict[sample_id + '_Zone6']['Probability'] == 1:
                print('Zone6')
            if label_dict[sample_id + '_Zone7']['Probability'] == 1:
                print('Zone7')
            if label_dict[sample_id + '_Zone8']['Probability'] == 1:
                print('Zone8')
            if label_dict[sample_id + '_Zone9']['Probability'] == 1:
                print('Zone9')
            if label_dict[sample_id + '_Zone10']['Probability'] == 1:
                print('Zone10')
        else:
            print('No threat in that portion')
        # Draw body detection
        data = rd.read_data(path)
        # we only get images from one certain angle
        img = np.rot90(data[:, :, ANGLE])
        # Create figure and axes
        fig = plt.figure(figsize=(img.shape[1]/100, img.shape[0]/100), dpi=100)
        ax = fig.add_axes([0., 0., 1., 1.])
        ax.imshow(img)
        # Firstly, blue rect is assumed to be
        fig.canvas.mpl_connect('button_press_event', button_press)
        fig.canvas.mpl_connect('button_release_event', button_release)
        fig.canvas.mpl_connect('key_press_event', key_press)
        plt.show()
        if release[0] < press[0]:
            temp = press
            press = release
            release = temp
        if release[1] < press[1]:
            temp = press[1]
            press[1] = release[1]
            release[1] = temp
        line = {'ID': sample_id, 'BlueStartX': press.copy()[0], 'BlueStartY': press.copy()[1],
                'BlueEndX': release.copy()[0], 'BlueEndY': release.copy()[1], 'RedStartX': [], 'RedStartY': [],
                'RedEndX': [], 'RedEndY': []}

        # Draw threat detection if exists
        if isThreatInBody and isRedExist[0]:
            while(isRedExist[0]):
                isRedExist[0] = False
                fig = plt.figure(figsize=(img.shape[1] / 100, img.shape[0] / 100), dpi=100)
                ax = fig.add_axes([0., 0., 1., 1.])
                ax.imshow(img)
                rect = patches.Rectangle((line['BlueStartY'], line['BlueStartX']),
                                          line['BlueEndY'] - line['BlueStartY'],
                                          line['BlueEndX'] - line['BlueStartX'],
                                          linewidth=1, edgecolor='b', facecolor='none')
                ax.add_patch(rect)
                fig.canvas.mpl_connect('button_press_event', button_press)
                fig.canvas.mpl_connect('button_release_event', button_release)
                fig.canvas.mpl_connect('key_press_event', key_press)
                plt.show()
                if release[0] < press[0]:
                    temp = press
                    press = release
                    release = temp
                if release[1] < press[1]:
                    temp = press[1]
                    press[1] = release[1]
                    release[1] = temp
                line['RedStartX'].append(press.copy()[0])
                line['RedStartY'].append(press.copy()[1])
                line['RedEndX'].append(release.copy()[0])
                line['RedEndY'].append(release.copy()[1])
        else:
            line['RedStartX'].append(-1)
            line['RedStartY'].append(-1)
            line['RedEndX'].append(-1)
            line['RedEndY'].append(-1)

        # Draw resulting figure
        fig = plt.figure(figsize=(img.shape[1]/100, img.shape[0]/100), dpi=100)
        ax = fig.add_axes([0., 0., 1., 1.])
        ax.imshow(img)
        if isThreatInBody:
            rect1 = patches.Rectangle((line['BlueStartY'], line['BlueStartX']),
                                          line['BlueEndY'] - line['BlueStartY'],
                                          line['BlueEndX'] - line['BlueStartX'],
                                     linewidth=1, edgecolor='b', facecolor='none')
            ax.add_patch(rect1)
            for i in range(len(line['RedStartX'])):
                ax.add_patch(patches.Rectangle((line['RedStartY'][i], line['RedStartX'][i]),
                                 line['RedEndY'][i] - line['RedStartY'][i],
                                 line['RedEndX'][i] - line['RedStartX'][i],
                                 linewidth=1, edgecolor='r', facecolor='none'))
        else:
            rect2 = patches.Rectangle((line['BlueStartY'], line['BlueStartX']),
                                      line['BlueEndY'] - line['BlueStartY'],
                                      line['BlueEndX'] - line['BlueStartX'],
                                 linewidth=3, edgecolor='b', facecolor='none')
            ax.add_patch(rect2)
        fig.canvas.mpl_connect('key_press_event', key_press_renew)
        plt.show()
        if not again[0]:
            index = index + 1
            dataInfo.append(line)

        # Check if the pixels are correct or not
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        press = [line['BlueStartX'], line['BlueStartY']]
        release = [line['BlueEndX'], line['BlueEndY']]
        for i in range(press[0], release[0], 1):
            for j in range(press[1] - 5, press[1] + 5, 1):
                img[i][j] = 0.00089142355
        for i in range(press[0], release[0], 1):
            for j in range(release[1] - 5, release[1] + 5, 1):
                img[i][j] = 0.00089142355
        for i in range(press[1], release[1], 1):
            for j in range(press[0] - 5, press[0] + 5, 1):
                img[j][i] = 0.00089142355
        for i in range(press[1], release[1], 1):
            for j in range(release[0] - 5, release[0] + 5, 1):
                img[j][i] = 0.00089142355
        if isThreatInBody:
            press = [line['RedStartX'], line['RedStartY']]
            release = [line['RedEndX'], line['RedEndY']]
            for i in range(press[0], release[0], 1):
                for j in range(press[1] - 5, press[1] + 5, 1):
                    img[i][j] = 0.00089142355
            for i in range(press[0], release[0], 1):
                for j in range(release[1] - 5, release[1] + 5, 1):
                    img[i][j] = 0.00089142355
            for i in range(press[1], release[1], 1):
                for j in range(press[0] - 5, press[0] + 5, 1):
                    img[j][i] = 0.00089142355
            for i in range(press[1], release[1], 1):
                for j in range(release[0] - 5, release[0] + 5, 1):
                    img[j][i] = 0.00089142355
        ax.imshow(img)
        fig.canvas.mpl_connect('key_press_event', key_press_renew)
        plt.show()
        '''
    columnList = ['ID', 'BlueStartX', 'BlueStartY',  'BlueEndX', 'BlueEndY', 'RedStartX', 'RedStartY', 'RedEndX', 'RedEndY']
    (pd.DataFrame(dataInfo)[columnList]).to_csv(clear_path + '/NewLabels.csv', mode='w', index=False)
# Function call
# init_data('initDataBatch1')


def e_neighbour_graph(img, e):
    array = np.array(img)
    rowN = np.shape(array)[0]
    colN = np.shape(array)[1]
    result = np.zeros((rowN * colN, rowN * colN))
    for row in range(rowN):
        for col in range(colN):
            for rowDiff in range(-e, e + 1, 1):
                for colDiff in range(-e, e + 1, 1):
                    if(row + rowDiff < rowN and -1 < row + rowDiff and col + colDiff < colN and -1 < col + colDiff):
                        if(rowDiff * rowDiff + colDiff * colDiff <= e * e):
                            result[row * colN + col][(row + rowDiff) * colN + (col + colDiff)] = 1
    return result


def fully_connected_graph(img, sigma):
    array = np.array(img)
    rowN = np.shape(array)[0]
    colN = np.shape(array)[1]
    result = np.zeros((rowN * colN, rowN * colN))
    for row in range(rowN):
        for col in range(colN):
            for row2 in range(rowN):
                for col2 in range(colN):
                    weight = math.exp(- math.pow(img[row][col] - img[row2][col2], 2) / (2 * math.pow(sigma, 2)))
                    result[row * colN + col][row2 * colN + col2] = weight
    return result


def fully_connected_e_neighbour_graph(img, sigma, e):
    array = np.array(img)
    rowN = np.shape(array)[0]
    colN = np.shape(array)[1]
    result = np.zeros((rowN * colN, rowN * colN))
    for row in range(rowN):
        for col in range(colN):
            for rowDiff in range(-e, e + 1, 1):
                for colDiff in range(-e, e + 1, 1):
                    if(row + rowDiff < rowN and -1 < row + rowDiff and col + colDiff < colN and -1 < col + colDiff):
                        if(rowDiff * rowDiff + colDiff * colDiff <= e * e):
                            row2 = row + rowDiff
                            col2 = col + colDiff
                            weight = math.exp(- math.pow(img[row][col] - img[row2][col2], 2) / (2 * math.pow(sigma, 2)))
                            result[row * colN + col][row2 * colN + col2] = weight
    return result
'''
# Check for graph implementations
ANGLE = 0 # 0-15, 0 is front, 8 is back
path_list = glob("./stage1_aps/*.aps")
path = path_list[0]
data = rd.read_data(path_list[0])
# we only get images from one certain angle
# img = np.rot90(data[:,:, ANGLE])
# imgplt = plt.imshow(img)
plt.show()


img2 = imread("./sample.jpg")
img2 = np.asarray(img2,np.dtype(np.int32))
# plt.imshow(img2)
# plt.show()
A = np.array([[1,2,3], [0,5,6]])
# A = img2
print("Image")
print(np.shape(A))
print(A)
print("e neighbour graph")
print(e_neighbour_graph(A, 1))
print("Fully connected & e neighbour graph")
print(np.round(fully_connected_e_neighbour_graph(A, 1, 1), decimals=2))
print("Fully connected graph")
print(np.round(fully_connected_graph(A, 1), decimals=2))
'''


def laplacian(graph):
    diag = np.diag(sum(graph))
    lap = diag - graph
    return [lap, diag]


def orthogonalize(diag, s):
    nodeNum = len(diag)
    p = np.ones((nodeNum, 1)) / pow(sum(np.diag(np.diag(sum(diag)))), .5)
    projection = np.dot(np.transpose(s), np.dot(diag, p))
    v = s - projection[0][0] * p
    # return [v, projection]
    return v


def PRwalk(graph, alpha, s):
    if alpha < 0:
        print('Error: alpha must be positive!')
        return -1
    nodeNum = len(graph)
    [L, D] = laplacian(graph)
    s = orthogonalize(D, s)
    v, flag = bicgstab((L + alpha * D), np.dot(D, s),tol=1e-3, maxiter=nodeNum)

    v = np.reshape(v, (len(v), 1))

    v = orthogonalize(D, v)
    return v

# Check for method implementations
# A = np.array([[1, 2, 3], [1, 5, 6]])
# G = np.round(fully_connected_e_neighbour_graph(A, 10, 1), decimals=1)
# print("Laplacian Check")
# [L, D] = laplacian(G)
# print(L)
# print(D)
# print("Orthogonalize check")
# A = np.array([2,3,3])
# D = np.diag(A)
# s = np.array([[1], [2], [3]])
# t = orthogonalize(D, s)
# print(t)
# print("PRwalk check")
# G = np.array([[0, 1], [1, 0]])
# s = np.array([[1], [2]])
# print(PRwalk(G, 1, s)) # output should be [[-0.167], [ 0.167]]

def normalize(D, s):
    norm = np.dot(np.dot(np.transpose(s), D), s)
    v = s/math.sqrt(norm)
    return v

myImg = np.asarray(imread("1a10297e6c3101af33003e6c4f846f47/1a10297e6c3101af33003e6c4f846f47.aps.jpg"), dtype=np.int)
# print(myImg)
# print(myImg.shape) # 148 row, 176 columns
graph = fully_connected_e_neighbour_graph(myImg, 5, 3)
seed = np.zeros((myImg.shape[0] * myImg.shape[1], 1))
seed[40 * myImg.shape[1] + 90][0] = 1 # Example seed
'''
myImg = np.asarray(imread("3a8696b99b2d1b28be62389d48d697be/3a8696b99b2d1b28be62389d48d697be.aps.jpg"), dtype=np.int)
graph = fully_connected_e_neighbour_graph(myImg, 5, 3)
# print(np.round(graph, decimals=3))
seed = np.zeros((myImg.shape[0] * myImg.shape[1], 1))
seed[20 * myImg.shape[1] + 65][0] = 1 # Example seed
'''

[L, D] = laplacian(graph)
s = orthogonalize(D, seed)
s = normalize(D, s)
alpha = .001
v = PRwalk(graph, alpha, s)
normalization = (np.dot(np.dot(np.transpose(v), D), v))[0][0]
print(normalization)
correlation = (np.dot(np.dot(np.dot(np.transpose(v), D), s), np.dot(np.dot(np.transpose(v), D), s)) / normalization)[0][0]
print(correlation)

while correlation < 1 / len(graph):
    alpha = 2 * alpha
    v = PRwalk(graph, alpha, s)
    normalization = (np.dot(np.dot(np.transpose(v), D), v))[0][0]
    correlation = \
        (np.dot(np.dot(np.dot(np.transpose(v), D), s), np.dot(np.dot(np.transpose(v), D), s)) / normalization)[0][0]
    print(correlation)

heat = np.zeros(myImg.shape)
for i in range(myImg.shape[0]):
    for j in range(myImg.shape[1]):
        heat[i][j] = v[i * myImg.shape[1] + j]

'''
# Saving on top of the cropped image
plt.imshow(myImg)
plt.imshow(heat, cmap='hot', interpolation='nearest', alpha=.5)
# plt.savefig('1a10297e6c3101af33003e6c4f846f47/RW_cropped_inThreat.eps')
'''

'''
# Saving on top of full image
# data = rd.read_data("1a10297e6c3101af33003e6c4f846f47/1a10297e6c3101af33003e6c4f846f47.aps")
img = np.rot90(data[:, :, 0])
myImgFull = np.asarray(img, dtype=np.double)
heat = np.zeros(myImgFull.shape)
for i in range(myImg.shape[0]):
    for j in range(myImg.shape[1]):
        # heat[204 + i][168 + j] = v[i * myImg.shape[1] + j] # 204, 168 are from annotated data
print(myImg.shape)
print(heat.shape)
plt.show()
plt.imshow(myImgFull)
plt.imshow(heat, cmap='hot', interpolation='nearest',alpha = .5)
# plt.savefig('1a10297e6c3101af33003e6c4f846f47/RW_FullImage_outsideofThreat_MiddlePointSeed.eps')
'''
Ps = np.copy(v[:,0])
Ps = Ps - min(Ps)
Ps = Ps / sum(Ps)
# Flattening image
imgArray = np.zeros(myImg.shape[0] * myImg.shape[1])
for i in range(myImg.shape[0]):
    for j in range(myImg.shape[1]):
        imgArray[i * myImg.shape[1] + j] = myImg[i][j]
LValues = np.argsort(-Ps)

plt.plot(Ps[list(LValues)])
plt.axis([0, 3500, 0, .0005])
plt.show()
# plt.savefig('ScoreCalculation/SeedOutThreadProbabilitySortingZoomed.eps')
# 1a10297e6c3101af33003e6c4f846f47
# L is chosen to be 2500 for seed in  threat
# L is chosen to be 3500 for seed out threat
# 3a8696b99b2d1b28be62389d48d697be
# L is chosen to be 2500 for seed in  threat
# L is chosen to be 3500 for seed out threat

numerator = 0
denominator = 0
for i in range(3500):
    numerator += Ps[LValues[i]] * imgArray[LValues[i]]
    denominator += Ps[LValues[i]]
score = numerator / math.sqrt(denominator)
# 1a10297e6c3101af33003e6c4f846f47
# Score for seed in  threat is 19.78153
# Score for seed out threat is 67.78561
# 3a8696b99b2d1b28be62389d48d697be
# Score for seed in  threat is 41.96559
# Score for seed out threat is 84.67606
