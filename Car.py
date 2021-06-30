import numpy as np
import math
from matplotlib import pyplot as plt
import copy

class Car():
	def __init__(self):
		self.load_mass_kg = 2000
		self.full_mass_kg = 1800
		self.all_mass_kg = self.load_mass_kg+self.full_mass_kg
		self.wheel_redis_m = 0.367
		self.efficient = 0.85
		self.rolling_risistance = 0.013
		self.air_area_m2 = 2.77
		self.i0 = 5.83
		self.If_kg_m2 = 0.218
		self.Iw1_kg_m2 = 1.798
		self.Iw2_kg_m2 = 3.598
		self.ig = np.array([5.56,2.769,1.644,1.00,0.793])
		self.L_m = 3.2
		self.a_m = 1.947
		self.b_m = self.L_m-self.a_m
		self.hg_m = 0.9
		self.n_min_r_min = 600
		self.n_max_r_min = 4000
		self.n = np.round(np.arange(self.n_min_r_min,self.n_max_r_min,0.01),2)
		self.ua = np.array([self.get_ua(i) for i in [0,1,2,3,4]]).T
		self.__ax = [-19.313,295.27,-165.44,40.874,-3.8445]
		self.__num = 5
		self.Tq_n_m = np.array([self.__tqn_n(t,self.__num,self.__ax) for t in self.n])
		self.__init_SS()

	def __init_SS(self):
		self.SS = 1+(1/self.full_mass_kg)*((self.Iw1_kg_m2+self.Iw2_kg_m2)+(self.If_kg_m2*(self.i0**2)*(self.ig**2)*self.efficient))/(self.wheel_redis_m**2)

	def __tqn_n(self,n,num,ax):
		i = 0
		Tq = 0.000
		while i<num:
			Tq += ax[i]*((n/1000))**i
			i+=1
		return Tq

	def get_ua(self,level):
		return 0.377*self.wheel_redis_m*self.n/(self.i0*self.ig[level])

	def get_n(self,level,ua):
		return ua*(self.i0*self.ig[level])/(0.377*self.wheel_redis_m)
	
	def get_Ft(self,level):
		return  self.Tq_n_m*self.i0*self.ig[level]*self.efficient/self.wheel_redis_m

	def get_Ff_Fw(self,ua):
		return 9.8*self.rolling_risistance*self.full_mass_kg+self.air_area_m2*(ua**2)/21.15

	def caculate_Ft_Ff_Fw_diagramm(self):
		fig_Ft_Ff_Fw = plt.figure()
		ax = fig_Ft_Ff_Fw.add_subplot(111)
		self.Ft = np.array([self.get_Ft(i) for i in [0,1,2,3,4]]).T
		ax.plot(self.ua,self.Ft)
		ua = np.linspace(0,self.ua[-1][-1],len(self.n))
		ax.plot(ua,self.get_Ff_Fw(ua),linestyle='--')
		self.Umax = copy.deepcopy(self.caculate_MaxUa())
		ax.plot([self.Umax[0],self.Umax[0]],[0,self.Umax[1]*2],linestyle='--')
		plt.annotate("Umax=%.2f"%self.Umax[0]+"km/h",xy=(self.Umax[0],self.Umax[1]),xytext=(-150,60),textcoords='offset points',fontsize=16,
             arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=.2'))
		self.imax = self.caculate_Maxi()
		ax.plot([self.imax[1],self.imax[1]],[0,self.imax[2]*1.3],linestyle='--')
		plt.annotate("imax=%.2f"%self.imax[0],xy=(self.imax[1],self.imax[2]),xytext=(50,0),textcoords='offset points',fontsize=16,
             arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=.2'))
		plt.annotate("Cmax=%.2f"%self.imax[3],xy=(self.imax[1],self.imax[2]),xytext=(50,30),textcoords='offset points',fontsize=16,
             arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=.2'))
		plt.show()

	def caculate_MaxUa(self):
		Ft5 = self.Ft[:,4]
		ua = self.ua[:,4]
		Ff_Fw = self.get_Ff_Fw(ua)
		i = 0
		Umax = ua[0]
		Ft_max = Ft5[0]
		while i<len(self.n):
			if Ft5[i]>=Ff_Fw[i]:
				Umax = copy.deepcopy(ua[i])
				Ft_max = copy.deepcopy(Ft5[i])
				i+=1
			elif Ft5[i]<Ff_Fw[i]:
				break
		return [Umax,Ft_max,i]

	def caculate_Maxi(self):
		Ft1 = self.Ft[:,0]
		ua = self.ua[:,0]
		Ff_Fw = self.get_Ff_Fw(ua)
		Fmax = Ft1-Ff_Fw
		fmax = np.max(Fmax)
		m = np.where(Fmax==fmax)
		amax = math.asin(fmax/(9.8*self.full_mass_kg))
		imax = math.tan(amax)
		return [imax,ua[m],Ft1[m],fmax/(9.8*self.full_mass_kg)]

	def caculate_a(self):
		Ff_Fw = self.get_Ff_Fw(self.ua[:,0:4])
		#print(self.Ft.shape,Ff_Fw.shape)
		a = (self.Ft[:,0:4]-Ff_Fw[:,0:4])/(self.SS[0:4]*self.full_mass_kg)
		Ft5 = self.Ft[:,4]
		ua = self.ua[:self.Umax[2],4]
		Ff_Fw_5 = self.get_Ff_Fw(ua)
		a_2 = (Ft5[:self.Umax[2]]-Ff_Fw_5)/(self.SS[4]*self.full_mass_kg)
		#print(a_2.shape)
		a_2_minux = 1/a_2
		a_minux = 1/a
		self.a_1 = copy.deepcopy([a_minux,a_2_minux,ua])
		return [a_minux,a_2_minux,ua]

	def diagramm_a(self):
		fig_a = plt.figure()
		ax = fig_a.add_subplot(111)
		a_1 = self.caculate_a()
		ax.plot(self.ua[:,:4],a_1[0])
		ax.plot(a_1[2],a_1[1])
		plt.ylim(0,a_1[0][-1][-1]*1.5)
		plt.show()

	def caculate_accleration_time(self,max_u):
		level = 0
		i = 0
		time = 0
		u_now = self.ua[i][level]
		while(u_now<=max_u):
			if i==len(self.n)-1:
				u_now = self.ua[i][level]
				n = round(self.get_n(level+1,u_now),2)
				i = np.where(self.n==n)
				i = i[0][0]
				level+=1
			elif i<len(self.n)-1:
				u_now = self.ua[i][level]
				time+=(self.ua[i+1][level]-u_now)*self.a_1[0][i][level]/3.6
				i+=1
		print(time)

	def caculate_changeLevel_ua(self):
		return 0