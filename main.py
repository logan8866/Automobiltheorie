from Car import Car
import numpy as np
from matplotlib import pyplot as plt

if __name__ == "__main__":
	car = Car()
	car.caculate_Ft_Ff_Fw_diagramm()
	car.diagramm_a()
	car.caculate_accleration_time(80)
	car.diagramm_pe()
	car.diagramm_Q()
	#car.from_u_get_level(50)
	car.diagramm_wirschaft_beschleunigen()
	car.caculate_6S_Q()
	car.diagramm_fei_E()
	car.caculate_break_s(30,0.8,0.02,0.02)
	car.caculate_break_s_when_kaputt('f',30,0.8,0.02,0.02)
	car.caculate_break_s_when_kaputt('r',30,0.8,0.02,0.02)