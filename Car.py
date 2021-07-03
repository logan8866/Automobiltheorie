import numpy as np
import math
from matplotlib import pyplot as plt
import copy
import time
from scipy import interpolate
import mpl_toolkits.mplot3d

class Car():
	def __init__(self):
		self.load_mass_kg = 2000
		self.full_mass_kg = 1800
		self.all_mass_kg = self.load_mass_kg+self.full_mass_kg
		self.wheel_redis_m = 0.367
		self.efficient = 0.85
		self.rolling_risistance = 0.013
		self.air_area_m2 = 2.77
		self.__init_verschiedene_i0()
		self.__init_car_i0()
		#self.i0 = 5.83
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
		self.__init_caculate()
		self.__init_break()

	def __init_caculate(self):
		self.n = np.round(np.arange(self.n_min_r_min,self.n_max_r_min,0.01),2)
		self.ua = np.array([self.__init_ua(i) for i in [0,1,2,3,4]]).T
		self.__ax = np.array([-19.313,295.27,-165.44,40.874,-3.8445])
		self.__num = 5
		#self.Tq_n_m = np.array([self.__fitting_n(t/1000,self.__num,self.__ax) for t in self.n])
		#self.Tq_n_m = self.__fitting_n(self.n/1000,self.__num,self.__ax)
		self.Tq_n_m = np.poly1d(np.flip(self.__ax))(self.n/1000)
		self.__init_SS()
		self.Ft = np.array([self.__init_Ft(i) for i in [0,1,2,3,4]]).T
		self.Ff_Fw = self.__init_Ff_Fw(self.ua)
		self.caculate_MaxUa()
		self.caculate_Maxi()
		self.caculate_a()
		self.__init_feulkonsum_charater()

	def __init_SS(self):
		self.SS = 1+(1/self.full_mass_kg)*((self.Iw1_kg_m2+self.Iw2_kg_m2)+(self.If_kg_m2*(self.i0**2)*(self.ig**2)*self.efficient))/(self.wheel_redis_m**2)

	def __init_ua(self,level):
		return 0.377*self.wheel_redis_m*self.n/(self.i0*self.ig[level])

	def get_u(self,level,n):
		return 0.377*self.wheel_redis_m*n/(self.i0*self.ig[level])

	def get_n(self,level,ua):
		return ua*(self.i0*self.ig[level])/(0.377*self.wheel_redis_m)
	
	def __init_Ft(self,level):
		return  self.Tq_n_m*self.i0*self.ig[level]*self.efficient/self.wheel_redis_m

	def __init_Ff_Fw(self,ua):
		return 9.8*self.rolling_risistance*self.full_mass_kg+self.air_area_m2*(ua**2)/21.15

	def caculate_Ft_Ff_Fw_diagramm(self):
		fig_Ft_Ff_Fw = plt.figure()
		ax = fig_Ft_Ff_Fw.add_subplot(111)
		ax.plot(self.ua,self.Ft)
		ua = np.linspace(0,self.ua[-1][-1],len(self.n))
		ax.plot(ua,self.__init_Ff_Fw(ua),linestyle='--')
		#self.Umax = copy.deepcopy(self.caculate_MaxUa())
		ax.plot([self.Umax[0],self.Umax[0]],[0,self.Umax[1]*2],linestyle='--')
		plt.annotate("Umax=%.2f"%self.Umax[0]+"km/h",xy=(self.Umax[0],self.Umax[1]),xytext=(-150,60),textcoords='offset points',fontsize=16,
             arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=.2'))
		#self.imax = self.caculate_Maxi()
		ax.plot([self.imax[1],self.imax[1]],[0,self.imax[2]*1.3],linestyle='--')
		plt.annotate("imax=%.2f"%self.imax[0],xy=(self.imax[1],self.imax[2]),xytext=(50,0),textcoords='offset points',fontsize=16,
             arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=.2'))
		plt.annotate("Cmax=%.2f"%self.imax[3],xy=(self.imax[1],self.imax[2]),xytext=(50,30),textcoords='offset points',fontsize=16,
             arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=.2'))
		plt.show()

	def caculate_MaxUa(self):
		Ft5 = self.Ft[:,4]
		ua = self.ua[:,4]
		Ff_Fw = self.Ff_Fw[:,4]
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
		self.Umax = [Umax,Ft_max,i]
		#return [Umax,Ft_max,i]

	def caculate_Maxi(self):
		Ft1 = self.Ft[:,0]
		ua = self.ua[:,0]
		Ff_Fw = self.Ff_Fw[:,0]
		Fmax = Ft1-Ff_Fw
		fmax = np.max(Fmax)
		m = np.where(Fmax==fmax)
		amax = math.asin(fmax/(9.8*self.full_mass_kg))
		imax = math.tan(amax)
		self.imax = [imax,ua[m],Ft1[m],fmax/(9.8*self.full_mass_kg)]
		#return [imax,ua[m],Ft1[m],fmax/(9.8*self.full_mass_kg)]

	def caculate_a(self):
		Ff_Fw = self.Ff_Fw[:,0:4]
		#print(self.Ft.shape,Ff_Fw.shape)
		a = (self.Ft[:,0:4]-Ff_Fw[:,0:4])/(self.SS[0:4]*self.full_mass_kg)
		Ft5 = self.Ft[:,4]
		ua = self.ua[:self.Umax[2],4]
		Ff_Fw_5 = self.Ff_Fw[:self.Umax[2],4]
		a_2 = (Ft5[:self.Umax[2]]-Ff_Fw_5)/(self.SS[4]*self.full_mass_kg)
		#print(a_2.shape)
		a_2_minux = 1/a_2
		a_minux = 1/a
		self.a_1 = copy.deepcopy([a_minux,a_2_minux,ua])
		#return [a_minux,a_2_minux,ua]

	def diagramm_a(self):
		fig_a = plt.figure()
		ax = fig_a.add_subplot(111)
		#a_1 = self.caculate_a()
		ax.plot(self.ua[:,:4],self.a_1[0])
		ax.plot(self.a_1[2],self.a_1[1])
		plt.ylim(0,self.a_1[0][-1][-1]*1.5)
		plt.show()

	def caculate_accleration_time(self,max_u):
		level = 0
		i = 0
		time = 0
		u_now = self.ua[i][level]
		if max_u>self.ua[-1][3]:
			print("too high speed!")
			return 0
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
		print("accleration time is: ",time)
		self.temp_level = level
		return time

	def caculate_changeLevel_ua(self):
		return 0

	def __init_feulkonsum_charater(self):
		self.n_manche = np.array([815,1207,1614,2012,2603,3006,3403,3804])
		self.__b_params = np.array([
		1326.8,-416.46,72.379,-5.8629,0.17768,
		1354.7,-303.98,36.657,-2.0553,0.043072,
		1284.4,-189.75,14.524,-0.51184,0.0068164,
		1122.9,-121.59,7.0035,-0.18517,0.0018555,
		1141.0,-98.893,4.4763,-0.091077,0.00068906,
		1051.2,-73.714,2.8593,-0.05138,0.00035032,
		1233.9,-84.478,2.9877,-0.047449,0.00028230,
		1129.7,-45.291,0.71113,-0.00075215,-0.000038568,]).reshape(8,5)
		self.__six_situation_params = [
		50,7.2,7.2,[25,25],0,
		200,16.7,23.9,[25,40],0.25,
		450,22.5,46.4,[40,40],0,
		625,14.0,60.4,[40,50],0.2,
		875,18.0,78.4,[50,50],0,
		1075,19.3,97.7,[50,25],-0.36]
		self.__init_pe()
		self.__init_bx()
		self.__init_b_pe_n_func()
		self.caculate_6S_Q()

	def __init_bx(self):
		self.b_n_x = []
		for i in range(len(self.n_manche)):
			self.b_n_x.append(np.poly1d(np.flip(self.__b_params[i,:])))

	def __init_b_pe_n_func(self):
		b_z = []
		Tq = np.arange(0,np.max(self.Tq_n_m),0.1)
		for i in range(len(self.n_manche)):
			b_z.append(self.b_n_x[i](Tq*self.n_manche[i]/9550))
		b_z = np.array(b_z)
		x,y = np.meshgrid(self.n_manche,Tq)
		self.func=interpolate.interp2d(x,y,b_z.T,kind='cubic')
		#look n-Tq grafik
		# ax=plt.subplot(111,projection='3d')
		# ax.plot_surface(x,y,self.func(self.n_manche,Tq))
		# ax.scatter(x,y,self.func(self.n_manche,Tq),s=.5,c='g')
		# plt.show()

	def caculate_Q_100km_h(self,level):
		Q = []
		for i in range(0,len(self.n),100):
			b = self.func(self.n[i],self.pe_b[i,level]*9550/self.n[i])[0]
			#Q.append(b)
			#print(b,self.n[i],self.pe[i,level])
			Q.append(self.pe_b[i,level]*b/(1.02*7.05*self.ua[i,level]))
		return np.array(Q)

	def from_u_get_level(self,u):
		index = np.min(np.where(np.round(self.ua-u,1)==0)[1])
		return index
		
	def caculate_Q(self,u,s):
		level = self.from_u_get_level(u)
		n = self.get_n(level,u)
		u_index = np.where(self.n==round(n,2))[0]
		b = self.func(n,self.pe_b[u_index,level])
		Q = b*(self.pe_b[u_index,level])*s/(102*7.07*u)
		return Q

	def caculate_aQ(self,u1,u2,a):
		t = 1/(3.6*a)
		Q = 0
		level1 = self.from_u_get_level(u1)
		level2 = self.from_u_get_level(u2)
		n = self.get_n(level1,u1)
		n_end = self.get_n(level2,u2)
		u_index = np.where(self.n==round(n,2))[0]
		while(n<=n_end):
			p = self.pe[u_index,level1]
			b = self.func(n,p)
			Q+=b*p/(367.1*7.02)
			u1+=1
			level1 = self.from_u_get_level(u1)
			n = self.get_n(level1,u1)
			u_index = np.where(self.n==round(n,2))[0]
		return Q

	def caculate_6S_Q(self):
		self.Q = 0
		self.Q+=self.caculate_Q(25,50)
		self.Q+=self.caculate_Q(40,200)
		self.Q+=self.caculate_Q(50,250)
		self.Q+=self.caculate_aQ(25,40,0.25)
		self.Q+=self.caculate_aQ(40,50,0.2)
		self.Q_6 = (self.Q/1075)*100
		print("Q:",self.Q,self.Q_6)


	def __init_pe(self):
		self.i = 0
		self.pf = (self.full_mass_kg*9.8*self.rolling_risistance*self.ua/3600)/self.efficient
		self.pw = (self.air_area_m2*self.ua**3/76140)/self.efficient
		self.pj = (self.ua*(self.Ft-self.Ff_Fw)/3600)/self.efficient
		self.pi = (self.full_mass_kg*9.8*self.i*self.ua/3600)/self.efficient
		self.pe = self.pf+self.pw+self.pj+self.pi
		self.pe_b = self.pf+self.pw
		self.pemax = np.mean(np.max(self.pe,axis=0))

	def diagramm_pe(self):
		fig_pe = plt.figure()
		ax = fig_pe.add_subplot(111)
		ax.plot(self.ua,self.pe)
		ax.plot([0,self.ua[-1][-1]],[self.pemax,self.pemax],linestyle="--")
		ua = ua = np.linspace(0,self.ua[-1][-1],len(self.n))
		pf = (self.full_mass_kg*9.8*self.rolling_risistance*ua/3600)/self.efficient
		pw = (self.air_area_m2*ua**3/76140)/self.efficient
		ax.plot(ua,pf+pw)
		plt.show()

	def diagramm_Q(self):
		fig_pe = plt.figure()
		ax = fig_pe.add_subplot(111)
		Q4 = self.caculate_Q_100km_h(3)
		Q5 = self.caculate_Q_100km_h(4)
		ax.plot(self.ua[::100,3],Q4)
		ax.plot(self.ua[::100,4],Q5)
		# ax.plot(Q5)
		# ax.plot(Q4)
		plt.show()

	def __init_verschiedene_i0(self):
		self.verschiedene_i0 = [5.17,5.43,5.83,6.17,6.33]

	def __init_car_i0(self,i=2):
		self.i0 = self.verschiedene_i0[i]

	def __reinit_car(self,i):
		self.__init_car_i0(i)
		self.__init_caculate()

	def __from_pe_n_get_b(self,n,pe):
		return self.func(n,pe)[0]

	def caculate_accleration_feulkonsum(self,u):
		u = self.get_u(0,self.n_min_r_min)
		u_end = u
		Q = 0
		while(u<=u_end):
			n = self.get_n(self.from_u_get_level(u),u)
			n_index = np.where(self.n==round(n,2))[0]
			b = self.__from_pe_n_get_b(n,self.pe_b[n_index,self.from_u_get_level(u)])
			Q += b*self.pe_b[n_index,self.from_u_get_level(u)]/(1.02*7.05*u*u_end)
			u+=1
		#print(Q[0])
		return Q[0]

	def diagramm_wirschaft_beschleunigen(self):
		u = 80
		t = []
		Q = []
		for i in range(len(self.verschiedene_i0)):
			self.__reinit_car(i)
			t.append(self.caculate_accleration_time(u))
			Q.append(copy.deepcopy(self.Q_6))
		fig_pe = plt.figure()
		ax = fig_pe.add_subplot(111)
		ax.plot(Q,t)
		ax.scatter(Q,t,c="g")
		plt.show()
	
	def __init_break(self):
		self.L_b = 3.95
		self.b_b = 0.38
		self.__init_null_break()
		self.__init_full_break()
		self.__init_fei_E()

	def __init_null_break(self):
		self.break_null_mass = 4080
		self.null_hg = 0.845
		self.null_a = 2.1
		self.null_b = self.L_b-self.null_a

	def __init_full_break(self):
		self.break_full_mass = 9290
		self.full_hg = 1.17
		self.full_a = 2.95
		self.full_b = self.L_b-self.full_a

	def __init_fei_E(self):
		z = np.arange(0,1,0.01)
		self.null_fei_f = self.b_b*z*self.L_b/(self.null_b+z*self.null_hg)
		self.full_fei_f = self.b_b*z*self.L_b/(self.full_b+z*self.full_hg)
		self.null_fei_r = (1-self.b_b)*z*self.L_b/(self.null_a-z*self.null_hg)
		self.full_fei_r = (1-self.b_b)*z*self.L_b/(self.full_a-z*self.full_hg)
		fei = np.arange(0,1,0.01)
		self.null_E_f = self.null_b/(self.b_b*self.L_b-fei*self.null_hg)
		self.null_E_r = self.null_a/((1-self.b_b)*self.L_b+fei*self.null_hg)
		self.full_E_f = self.full_b/(self.b_b*self.L_b-fei*self.full_hg)
		self.full_E_r = self.full_a/((1-self.b_b)*self.L_b+fei*self.full_hg)

	def diagramm_fei_E(self):
		z = np.arange(0,1,0.01)
		fei = np.arange(0,1,0.01)
		fei_index = np.where(np.round(self.full_E_r-self.full_E_f,2)==0)[0][0]
		fig_1 = plt.figure()
		ax1 = fig_1.add_subplot(111)
		ax1.plot(z,self.null_fei_f,c="g")
		ax1.plot(z,self.full_fei_f)
		ax1.plot(z,self.null_fei_r)
		ax1.plot(z,self.full_fei_r)
		ax1.text(z[50],self.null_fei_f[50],'null f')
		ax1.text(z[50],self.full_fei_f[50],'full f')
		ax1.text(z[70],self.null_fei_r[70],'null r')
		ax1.text(z[70],self.full_fei_r[70],'full r')
		plt.show()
		fig_2 = plt.figure()
		ax = fig_2.add_subplot(111)
		#ax.plot(fei,self.null_E_f)
		ax.plot(fei,self.null_E_r)
		ax.plot(fei[:fei_index+1],self.full_E_f[:fei_index+1])
		ax.plot(fei[fei_index:],self.full_E_r[fei_index:])
		#ax.text(fei[50],self.null_E_f[50],'null Ef')
		ax.text(fei[50],self.full_E_f[50],'full Ef')
		ax.text(fei[70],self.null_E_r[70],'null Er')
		ax.text(fei[70],self.full_E_r[70],'full Er')
		ax.plot([fei[fei_index],fei[fei_index]],[0,1.2],linestyle="--")
		ax.text(fei[fei_index],0.2,'fei='+str(fei[fei_index]))
		plt.ylim(0,1.2)
		plt.show()

	def caculate_break_s(self,u,fei,t21,t22):
		fei_1 = np.arange(0,1,0.01)
		fei_index = np.where(fei_1==fei)[0][0]
		null_E = min(self.null_E_f[fei_index],self.null_E_r[fei_index])
		s = (t21+t22/2)*u/3.6+u**2/(25.92*9.8*null_E*fei)
		print("null mass s is:"+str(s))
		full_E = min(self.full_E_f[fei_index],self.full_E_r[fei_index])
		s = (t21+t22/2)*u/3.6+u**2/(25.92*9.8*full_E*fei)
		print("full mass s is:"+str(s))

	def caculate_break_s_when_kaputt(self,rof,u,fei,t21,t22):
		if rof == 'f':
			null_z = self.null_a*fei/(self.L_b+fei*self.null_hg)
			full_z = self.full_a*fei/(self.L_b+fei*self.full_hg)
		if rof == 'r':
			null_z = self.null_b*fei/(self.L_b-fei*self.null_hg)
			full_z = self.full_b*fei/(self.L_b-fei*self.full_hg)
		s = (t21+t22/2)*u/3.6+u**2/(25.92*9.8*null_z)
		print("when "+rof+"_break is down, null_break_s is: "+str(s))
		s = (t21+t22/2)*u/3.6+u**2/(25.92*9.8*full_z)
		print("when "+rof+"_break is down, full_break_s is: "+str(s))













	

