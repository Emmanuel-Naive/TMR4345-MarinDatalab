import numpy as np
#初始：地址
def ReadCsvData(filepath,xrow,ycol):
    file = open(filepath, "rb")
    filedata=np.loadtxt(file, delimiter=";")
    file.close()
    filearray=np.array(filedata)

    dataname = np.zeros(shape=(xrow,ycol))
    for i in range(xrow):
        for j in range(ycol):
            dataname[i,j]=filearray[i, j]
    return dataname

# filepath = "E:/User/Desktop/datalab/u-wind1.csv"
# wind = ReadCsvData(filepath,10,10)
# print(wind)