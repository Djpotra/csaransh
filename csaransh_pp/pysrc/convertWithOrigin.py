import os


class ConvertWithOrigin():
    
    def __init__(self):

        self.latConst = 2.87  # 晶格常数
        self.edge = 80  # box大小

    def __readAtoms(self, filename):
        '''
        读取xyz坐标，并存到列表中
        :param filename:
        :return:
        '''
        atoms_xyz = []

        # 读取xyz文件，并将数据写入到列表中存储
        with open(filename, 'r') as f:
            f.readline()  # 读掉前两行
            f.readline()
            line = f.readline()
            line = line.strip('\n')  # 去除每次读到的换行符

            while line != "":
                each_line = line.split()  # 将读到的字符串以空格隔开
                atom = [float(each_line[1]), float(each_line[2]), float(each_line[3])]  # 只取三维坐标
                atoms_xyz.append(atom)  # 将该字典加入列表
                line = f.readline()
                line = line.strip('\n')

        return atoms_xyz

    def convertWithOrigin(self, stepFile, originFile,outFile):

        stepXYZ = self.__readAtoms(stepFile)
        originXYZ = self.__readAtoms(originFile)

        sizeX = 8
        sizeY = 10
        timeX = int(self.edge / sizeX)
        timeY = int(self.edge / sizeY)
        index = 0

        with open(outFile,"w") as f:
            f.write("1024000\n")
            f.write("CCCCCCCCCCCCCCCCCCCCCC\n")
            for k in range(timeX):
                for t in range(timeY):
                    for z in range(80):
                        for y in range(t * sizeY, t * sizeY + sizeY):
                            for x in range(k * sizeX, k * sizeX + sizeX):
                                print("converting  "+str(index)+" in 1024000",end="\r")
                                vertex = []
                                vertex.append(stepXYZ[index][0] - originXYZ[index][0] + x*self.latConst)
                                vertex.append(stepXYZ[index][1] - originXYZ[index][1] + y*self.latConst)
                                vertex.append(stepXYZ[index][2] - originXYZ[index][2] + z*self.latConst)
                                f.write("Fe  "+str(round(vertex[0],6))+"  "+str(round(vertex[1],6))+"  "+str(round(vertex[2],6))+"\n")
                                index += 1
                                
                                center = []
                                center.append(stepXYZ[index][0] - originXYZ[index][0] + (x+0.5)*self.latConst)
                                center.append(stepXYZ[index][1] - originXYZ[index][1] + (y+0.5)*self.latConst)
                                center.append(stepXYZ[index][2] - originXYZ[index][2] + (z+0.5)*self.latConst)
                                f.write("Fe  "+str(round(center[0],6))+"  "+str(round(center[1],6))+"  "+str(round(center[2],6))+"\n")
                                index += 1
            f.close()
            print("\ndone\n")

if __name__ == "__main__":
    cwo = ConvertWithOrigin()
    stepDir = "/home/zxy/csaransh/csaransh_pp/data/beforeConvert"
    originPath = "/home/zxy/csaransh/csaransh_pp/data/origin.xyz"

    for dirpath, _, filenames in os.walk(stepDir):
        for filename in filenames:
            stepPath = "/home/zxy/csaransh/csaransh_pp/data/beforeConvert/"+filename
            outPath = "/home/zxy/csaransh/csaransh_pp/data/generic/out_"+filename
            print("converting "+stepPath)
            cwo.convertWithOrigin(stepPath,originPath,outPath)