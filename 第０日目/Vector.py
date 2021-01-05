class Vector3:
    def __init__(self,x=0,y=0,z=0):
        self.__x=x
        self.__y=y
        self.__z=z

    def getx(self):
        return(self.__x)

    def gety(self):
        return(self.__y)

    def getz(self):
        return(self.__z)

    def givex(self,x):
        self.__x=x

    def givey(self,y):
        self.__y=y

    def givez(self,z):
        self.__z=z

    def getv(self):
        return([self.__x,self.__y,self.__z])
