import Vector
import numpy as np
import math

class RK4:
	def __init__(self,dt=0.01,r=Vector.Vector3(0,0,0),v=Vector.Vector3(0,0,0),a=Vector.Vector3(0,0,0)):
		self.__r=r
		self.__v=v
		self.__dt=dt
		self.__a=a

	def rk4(self,dir):
		if dir==0:
			xp=self.__r.getx()
			vp=self.__v.getx()
			a=self.__a.getx()
		elif dir==1:
			xp=self.__r.gety()
			vp=self.__v.gety()
			a=self.__a.gety()
		else:
			xp=self.__r.getz()
			vp=self.__v.getz()
			a=self.__a.getz()

		k1=vp
		j1=a(xp,vp)
		k2=vp+0.5*j1*self.__dt
		j2=a(xp+0.5*k1*self.__dt,vp+0.5*j1*self.__dt)
		k3=vp+0.5*j2*self.__dt
		j3=a(xp+0.5*k2*self.__dt,vp+0.5*j2*self.__dt)
		k4=vp+j3*self.__dt
		j4=a(xp+k3*self.__dt,vp+j3*self.__dt)
		k=(k1+2*k2+2*k3+k4)/6
		j=(j1+2*j2+2*j3+j4)/6

		xp=xp+k*self.__dt
		vp=vp+j*self.__dt
		return(xp,vp)

	def rk4x(self):
		x,v=self.rk4(0)
		self.__r.givex(x)
		self.__v.givex(v)
		return(x,v)

	def rk4y(self):
		x,v=self.rk4(1)
		self.__r.givey(x)
		self.__v.givey(v)
		return(x,v)

	def rk4z(self):
		x,v=self.rk4(2)
		self.__r.givez(x)
		self.__v.givez(v)
		return(x,v)
