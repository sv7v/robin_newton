#!/usr/bin/env python
'''
Project 'Robin-Newton' (RN)
See:
1. Тестовое_задание_по_численным_методам.docx
2. Части "математическая постановка задачи" и "численный метод" по ссылке: https://docs.google.com/document/d/1yJ9DMjvC3wu9AfAU4ENC8eFZgz_kq7e430SLJj0Sjt8/edit?usp=sharing
3. Копия условия задачи: https://docs.google.com/document/d/1Y_SM1xSkq_DAGMlGw02bsXH4HyZxXd_Z/edit?usp=sharing&ouid=110237950124914475014&rtpof=true&sd=true
'''
from math       import pi
from itertools  import count

class Behavior:
	'''Поведение модели.'''
	def __init__(s, rad1,           # R₁  радиус
	                rad2,           # R₂
	                dtime_s,        # ΔT  время от старта до окончания в секундах
			dtime_int,      #                  ---//---        в квантах (штуках)
			mass_vol,       # ρ   плотность (масса на объём)     кг/м^3
			heatCap_mass,   # Cᵨ  удельная теплоёмкость          Дж/(кг*К)
# TO DO Забыт Кельвин в условии задачи.
			termCond,       # λ   теплопроводность               Вт/(м*К)
			n,              #     количество слоёв (ячеек сетки)
			alpha1,         # α₁
			alpha2,         # α₂
			temp1,          # T₁  температура на внутренней стенке трубы
			temp2):         # T₂
		'''Передача в модель всего её описания.'''
		if 0 < n: 'Ok'
		else:      raise Exception()
		s.n = n

		s.dr = (rad2 - rad1) / n                # толщина слоя (ячейки)
		s.rad1 = rad1
		s.rad2 = rad2

		per1 = 2*pi * ((rad1 + s.dr) / 2)       # обхват (периметр) слоя крайнего с края 1
		V1_div = s.dr * per1                    # объём этого слоя поделённый на длину трубы (далее -- ПНДТ или _div)
		m1_div = V1_div * mass_vol              # масса этого слоя ПНДТ
		s.heatCap1_div = heatCap_mass * m1_div  # теплоёмкость этого слоя ПНДТ

		s.dt = dtime_s / dtime_int              # величина кванта времени в секундах

		s.termCond = termCond
		s.alpha1 = alpha1
		s.alpha2 = alpha2
		s.t1 = temp1
		s.t2 = temp2
	#def __init__():

	def _NewtonRichmann_div(alpha, t0, t1, r):
		'''Мощность текущая между внешним миром и крайней
		ячейкой ПНДТ -- закон Ньютона-Рихмана.'''
		return alpha * (t0 - t1) * 2 * pi * r

	# Мощность ПОЛученная 0-й  ячейкой от внешней среды ПНДТ.
	# Мощность ОТДанная   n-ой ячейкой    внешней среде ПНДТ.
	def P1_NewtonRichmann_div(s, t): return Behavior._NewtonRichmann_div(s.alpha1, s.t1,   t , s.rad1)
	def P2_NewtonRichmann_div(s, t): return Behavior._NewtonRichmann_div(s.alpha2,   t , s.t2, s.rad2)

	def power_div(s, i,             #
	                 temp_i,        # температура  i-ой  ячейки
	                 temp_next):    #   ---//---  i+1-ой
		'''Мощность передаваемая от i-ой ячейке к следующей.
		Пропорциональна коэффициенту теплопроводности и площади.
		Обратно пропорциональна шагу сетки.'''

		# площадь торцов слоя
		#    i   |  -1  |  n-1
		# скобки | rad1 | rad2
		S_div = (s.rad1 + (s.rad2-s.rad1)*((i+1)/s.n)) * 2*pi

		return (temp_next - temp_i) * S_div * s.termCond / s.dr

	def dtemp(s, i,         # i-я ячейка (от 0 до n-1)
	             P_prev,    # мощность ПОЛучаемая   от предыдущей ячейки
	             P_next):   # ---//--- ПЕРЕдаваемая     следующей
		'''Изменение температуры ячейки.'''
		# Отношение радиусов i-ой и 0-й ячеек. Оно же равно
		# оношению объёмов и отношению масс этих ячеек.
		b = (s.rad1 + s.dr*(0.5+i)) / (s.rad1 + s.dr/2)
		return (P_prev - P_next) * s.dt / (s.heatCap1_div * b)
#class Behavior:

class pair_iterator():
	def __init__(self, i, buff, flag_curr, flag_my):
		self._i    = i
		self._buff = buff
		self._flag = flag_curr
		self._my   = flag_my

	def __iter__(self):
		return self

	def __next__(self):
		if   0 == len(self._buff):
			ret = next(self._i)
			self._buff.append(ret)
			self._flag[0] = not self._my
		elif self._flag[0] == self._my:
			ret = self._buff.pop(0)
		else:
			ret = next(self._i)
			self._buff.append(ret)
		return ret
	#def __next__(self):
#class pair_iterator():

def pair_iter(iterable):
	i    = iterable
	buff = []
	flag = [None]
	return (pair_iterator(i, buff, flag, True ),
	        pair_iterator(i, buff, flag, False))

def pair(x):
	ret = pair_iter(x)
	next(ret[1])
	return zip(ret[0], ret[1])

class State:
	'''Состояние. Время и температура.'''
	def __init__(self,
	             temp,     # температура ячеек
	             beh):     # поведение модели
		self._time = 0
		self._temp = temp
		self._beh  = beh

	def _power_div(self, t):
		'''Мощность получаемая ячейкой.'''
		t0, t1 = next(t)
		yield self._beh.P1_NewtonRichmann_div(t1)

		yield self._beh.power_div(0, t0, t1)

		for i,(t_p,t_n) in zip(count(1), t):
			yield self._beh.power_div(i, t_p, t_n)

		yield self._beh.P2_NewtonRichmann_div(t_p)

	def _dtemp(self, P):
		'''Изменение температуры.'''
		for i, (P_p, P_n) in zip(count(), pair(P)):
			yield self._beh.dtemp(i, P_p, P_n)

	def getNext(self):
		'''Возвращает состояние в следующий момент времени.'''
		t1, t2 = pair_iter(self._temp)

		ret = State(map(sum, zip(t2,
		                         self._dtemp(self._power_div(pair(t1))))),
		            self._beh)
		ret._time = self._time + 1

		return ret
	def get(self):
		return self._temp
#class State:

def run(rad1,           # R₁
        rad2,           # R₂
        dtime_s,        # ΔT
        dtime_int,      #
        mass_vol,       # ρ
        heatCap_mass,   # Cᵨ
        termCond,       # λ
        n,
        temp0,          # T₀ᵢ
        alpha1,         # α₁
        alpha2,         # α₂
        temp1,          # T₁
        temp2):         # T₂

	temp0 = iter(temp0)
	b = Behavior(rad1         = rad1,
	             rad2         = rad2,
	             dtime_s      = dtime_s,
	             dtime_int    = dtime_int,
	             mass_vol     = mass_vol,
	             heatCap_mass = heatCap_mass,
	             termCond     = termCond,
	             n            = n,
	             alpha1       = alpha1,
	             alpha2       = alpha2,
	             temp1        = temp1,
	             temp2        = temp2)

	s = State(map(lambda i:sum(i)/2,
	              pair(temp0)),
	          b)
	for i in range(dtime_int):
		s = s.getNext()

	yield temp1

	for t in map(lambda i:sum(i)/2,
	             pair(s.get())):
		yield t

	yield temp2
#def run():

def main():
	from sys        import argv
	from json       import load
	from csv        import reader

	from sys        import getrecursionlimit
	from sys        import setrecursionlimit

	name_param = argv[1]
	data_name  = argv[2]

	with open(name_param, 'r') as fd_param:
		param = load(fd_param)
		with open(data_name, 'r') as fd_data:
			temp0 = reader(fd_data, delimiter=';')
			temp0 = next(temp0)
			temp0 = map(float,temp0)

			ret = run(param['R_1'],
			          param['R_2'],
			          param['Deltat'],
				  param['dtime_int'],
			          param['ro'],
			          param['C_ro'],
			          param['lambda'],
				  param['n'],
				  temp0,
			          param['alpha_1'],
			          param['alpha_2'],
			          param['T_1'],
			          param['T_2'])

			rl = getrecursionlimit()
			if rl/2 < param['dtime_int']:
				setrecursionlimit(2 * param['dtime_int'])

			print(';'.join(map(str, ret)))

if __name__ == "__main__":
	main()
