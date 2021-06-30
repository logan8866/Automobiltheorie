from Car import Car
import numpy as np
from matplotlib import pyplot as plt

if __name__ == "__main__":
	car = Car()
	#print(car.ua.shape)
	car.caculate_Ft_Ff_Fw_diagramm()
	car.diagramm_a()
	car.caculate_accleration_time(100)
	
