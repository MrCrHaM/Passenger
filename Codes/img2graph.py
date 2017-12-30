import numpy

def e_neighbour_graph(img, e):
    array = numpy.array(img)
    rowN = numpy.shape(array)[0]
    colN = numpy.shape(array)[1]
    result = numpy.zeros((rowN * colN, rowN * colN))
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
