import numpy as np
import gdal
import math
from PIL import Image
#0.读取parameters存储在变量中
filePath = 'data/'
fileName = 'parameter.txt'
txtData = [line for line in open(filePath+fileName)]
P = []
#起始五个参考点坐标
listX =[] 
listY =[]
listZ =[]
def calSateVec(coeX, coeY, coeZ,t):
    x = coeX[0]+coeX[1]*t+coeX[2]*pow(t,2)+coeX[3]*pow(t,3)
    y = coeY[0]+coeY[1]*t+coeY[2]*pow(t,2)+coeY[3]*pow(t,3)
    z = coeZ[0]+coeZ[1]*t+coeZ[2]*pow(t,2)+coeZ[3]*pow(t,3)
    vx = coeX[1]+2*coeX[2]*t+3*coeX[3]*pow(t,2)
    vy = coeY[1]+2*coeY[2]*t+3*coeY[3]*pow(t,2)
    vz = coeZ[1]+2*coeZ[2]*t+3*coeZ[3]*pow(t,2)
    ax = 2*coeX[2]+6*coeX[3]*t
    ay = 2*coeY[2]+6*coeY[3]*t
    az = 2*coeZ[2]+6*coeZ[3]*t

    return [float(x),float(y),float(z),float(vx),float(vy),float(vz),float(ax),float(ay),float(az)]

for i in range(5):#0-4 rows
    P.append([])
    [tempX,tempY,tempZ] = txtData[i].strip().split(' ') 
    # P[i].append(tempX)
    # P[i].append(tempY)
    # P[i].append(tempZ)
    listX.append(float(tempX))
    listY.append(float(tempY))
    listZ.append(float(tempZ))
matX = np.mat(listX).reshape((5,1))
matY = np.mat(listY).reshape((5,1))
matZ = np.mat(listZ).reshape((5,1))
#第一点成像时间和时间间隔
[t_initial,t_delta] = txtData[6].strip().split(' ')
#四个角点经纬度
x_corner = []
y_corner = []
for i in range(8,12):
    [tempX,tempY] = txtData[i].strip().split(' ')
    x_corner.append(tempX)
    y_corner.append(tempY)
#采样频率
[t_firstRow,t_freRow] = txtData[13].strip().split(' ')
[t_firstCol,t_freCol] = txtData[14].strip().split(' ')
#投影坐标系长短半轴
[AxisLengthA, AxisLengthB] = txtData[-1].strip().split(' ')
 
#1. 根据五个参考点坐标和对应成像时间，计算卫星位置矢量的函数参数。
#attention
listT = []
for i in range(0,5):
    listT.append(0+float(t_delta)*i)#attention: if there 0 -- float(t_initial),cal 真实值，而不是相对值。0的话后面输入t_old就是t_old-t_initial
#arrayT is a 1*20 list, the same as T in word
arrayT = []
for i in range(0,4*5):
    if i%4 == 0 :
        arrayT.append(1)
        continue
    arrayT.append(pow(listT[int(i/4)],i%4))
matT = np.mat(arrayT).reshape((5,4))
#coeX,Y,Z save coes of formation. 4*1,a0,a1,a2,a3 matches t^0 - t^3
coeX = ((matT.T*matT).I*(matT.T))*matX
coeY = (matT.T*matT).I*(matT.T)*matY
coeZ = (matT.T*matT).I*(matT.T)*matZ


#2. Dem内插，90m-- 30m；经纬度在具体地像素处理时计算

from gdalconst import *
fileDEM = 'dem.img'#attention: read img directly here
dataset = gdal.Open(filePath+fileDEM,GA_ReadOnly)
geoTransform = dataset.GetGeoTransform()
heightDEM = dataset.RasterYSize
widthDEM = dataset.RasterXSize
band = dataset.GetRasterBand(1)


fileSAR = 'sar.001'
imgSAR = gdal.Open(filePath+fileSAR,GA_ReadOnly)
bandSAR = imgSAR.GetRasterBand(1)

heightSAR = imgSAR.RasterYSize
widthSAR = imgSAR.RasterXSize


#3. 多普勒t公式，对于dem块（利用角点画出minX maxX，minY，maxY）的矩形框，并且遍历所有像元，得到DEM相应的SAR坐标
#判断DEM上的角点坐标 (minx,miny)  (maxX,maxY)，可以大大减少计算量。


#attention:椭圆的长短半轴怎么用（难道要经纬度转化为大地坐标？）

pi = 3.1415926
#遍历dem.img(加密采样后)的所有像元，每个像元存储着经纬度坐标（理想的情形）;然后画框。减少计算量


maxX = int((float(min(x_corner))-geoTransform[3])/geoTransform[5])
minY = int((float(min(y_corner))-geoTransform[0])/geoTransform[1])
minX = int((float(max(x_corner))-geoTransform[3])/geoTransform[5])
maxY = int((float(max(y_corner))-geoTransform[0])/geoTransform[1])

#for i in range(widthDEM):
#    for j in range(heightDEM):
#        latitudeDEM = (geoTransform[3] + geoTransform[5]*i)*pi/180
#        longtitudeDEM = (geoTransform[0]+j*geoTransform[1])*pi/180
#        if latitudeDEM == min(x_corner) and longtitudeDEM == min(y_corner):
#            [minX,minY] = [i,j]
#        if latitudeDEM == max(x_corner) and longtitudeDEM == max(y_corner):
#            [maxX,maxY] = [i,j]

#在DEM上画完框后，开始计算方位向和距离向的时间;并且内插出每个DEM像元点（i，j）对应的修正辐射值。
#新的图像的坐标关系为：(i,j)  ---  (i-minX,j-minY)

newSAR = Image.new('L',(maxX-minX+1,maxY-minY+1))
#attention:有可能找不到绝对一样的经纬度点，可能直接DEMx；DEMy的全图像了。
AxisLengthA = float(AxisLengthA)*1000
AxisLengthB = float(AxisLengthB)*1000

#minX += 2000
#minY += 2500
#maxX = minX + 10
#maxY = minY + 10
minX = 0
minY = 0
for i in range(minX,maxX+1):
    for j in range(minY,maxY+1):

        #获取灰度值，计算高程,这里是每行来读取的grey数据，因为整幅图读出来没问题，只是Python需要解码，这个二进制码太长了。。。会崩
        scanLine = band.ReadRaster(0,i,band.XSize,1,band.XSize,1,GDT_Float32)
        import struct
        greyDEM = struct.unpack('f'*band.XSize,scanLine)
        #因为需要双线性插值，要两行才行。
        greySAR1 = bandSAR.ReadRaster(0,i,bandSAR.XSize,1,bandSAR.XSize,1,GDT_Float32)
        greySAR1 = struct.unpack('f'*bandSAR.XSize,greySAR1)

        greySAR2 = bandSAR.ReadRaster(0,i+1,bandSAR.XSize,1,bandSAR.XSize,1,GDT_Float32)
        greySAR2 = struct.unpack('f'*bandSAR.XSize,greySAR2)

        latitudeDEM = (geoTransform[3] + geoTransform[5]*i)*pi/180
        longtitudeDEM = (geoTransform[0]+j*geoTransform[1])*pi/180
        hDEM = greyDEM[j] #attention:
        e = math.sqrt(float(AxisLengthA)**2-float(AxisLengthB)**2)/float(AxisLengthA)
        N = float(AxisLengthA)/math.sqrt(1-(e*math.sin(latitudeDEM))**2)
        xDEM = (N+hDEM)* math.cos(latitudeDEM)*math.cos(longtitudeDEM)
        yDEM = (N+hDEM)* math.cos(latitudeDEM)*math.sin(longtitudeDEM)
        zDEM = (N*(1-e**2)+hDEM)*math.sin(latitudeDEM)

        #计算初始迭代时刻，中间行成像时间
        t_mid = float(t_firstRow) + int(heightSAR/2)*(1/float(t_freRow))
        t_old = t_mid
        #attention:
        t_old = 2*float(t_delta)
        deltaT = 1
        num_iteration = 0
        #开始迭代，求出实际此DEM坐标对应的方位向t
        while (num_iteration < 10 and abs(deltaT) >= 1.0e-15 ):
            [x,y,z,vx,vy,vz,ax,ay,az] = calSateVec(coeX,coeY,coeZ,t_old)
            #attention, delta mostly has pro.
            deltaT = (vx*(x-xDEM)+vy*(y-yDEM)+vz*(z-zDEM))/((ax*(x-xDEM)+ay*(y-yDEM)+az*(z-zDEM))+vx*vx+vy*vy+vz*vz)
            t_old -= deltaT
            num_iteration += 1
        #根据距离向和方位向时间计算SAR的行列号
        i_SAR = (t_old + float(t_initial)- float(t_firstRow))*(float(t_freRow))#attention:ori: t_freRow*2
        trueSateVars = []
        trueSateVars = calSateVec(coeX,coeY,coeZ,t_old)
        R = (trueSateVars[0] - xDEM)**2 + (trueSateVars[1] - yDEM)**2 +(trueSateVars[2] - zDEM)**2
        
        delta_tRange = (float(t_freCol)*1.0e+6)#attention:MHZ:1.0e-6;rect to 1.0e-3
        j_SAR = ((2*math.sqrt(R)/299792458.458)- float(t_firstCol)*1.0e-3)*delta_tRange

        #4.内插此坐标，得到灰度，生成新影像
        #双线性内插出DEM像元（i，j）对应的实际辐射值
        if (int(i_SAR) or int (i_SAR)+1) not in range(0,heightSAR):
            newPixel = 0
            #print('ding!')
            newSAR.putpixel((i-minX,j-minY),newPixel)
            continue
        if (int(j_SAR) or int (j_SAR)+1) not in range(0,widthSAR):
            newPixel = 0
            #print('ding!')
            newSAR.putpixel((i-minX,j-minY),newPixel)
            continue
        ratioX = i_SAR - int(i_SAR)
        ratioY = j_SAR - int(j_SAR)
        #attention:这里的像素值，可能需要复数计算一下。
        if(int(j_SAR) >= len(greySAR1) or int(j_SAR) >= len(greySAR2) or int(j_SAR)+1 >=len(greySAR1) or  int(j_SAR)+1 >=len(greySAR2)):
            newPixel = 0
            newSAR.putpixel((i-minX,j-minY),newPixel)
            continue

        #pixelLeftUp = greySAR1[int(j_SAR)]
        #pixelRightUp = greySAR2[int(j_SAR)]
        #pixelLeftDown = greySAR1[int(j_SAR)+1]
        #pixelRightDown = greySAR2[int(j_SAR)+1]

        pixelLeftUp = greySAR1[2*int(j_SAR)]**2+greySAR1[2*int(j_SAR)+1]**2
        pixelRightUp = greySAR2[2*int(j_SAR)]**2+greySAR2[2*int(j_SAR)+1]**2
        pixelLeftDown = greySAR1[2*(int(j_SAR)+1)]**2+greySAR1[2*(int(j_SAR)+1)+1]**2
        pixelRightDown = greySAR2[2*(int(j_SAR)+1)]**2+greySAR2[2*(int(j_SAR)+1)+1]**2

        newPixel = ratioY * (ratioX*pixelLeftUp+(1-ratioX)*pixelRightUp) + (1-ratioY)*(ratioX*pixelLeftDown+(1-ratioX)*pixelRightDown)
        #attention 可能这个newPixel需要int一下。
        newSAR.putpixel((i-minX,j-minY),int(newPixel))
        
        #print('ding!ding!ding!')
fileNewSAR = 'newSAR'
newSAR.save(filePath+fileNewSAR+'.bmp','bmp')
newSAR.show()
print('see result %s.bmp' % fileNewSAR )

