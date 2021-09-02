import numpy as np
#初始：地址
filepath = "E:/User/Desktop/datalab/u-wind1.csv"

file = open(filepath, "rb")
filedata=np.loadtxt(file, delimiter=";")
file.close()
filearray=np.array(filedata)
#初始：行数，列数
xrow =15
ycol = 15

dataname = np.zeros(shape=(xrow,ycol))
for i in range(xrow):
    for j in range(ycol):
        dataname[i,j]=filearray[i, j]
print(dataname[1,1])