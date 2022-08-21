from re import findall
from subprocess import getoutput
from numpy import array,sqrt,sum

global p_i,p_s,p_e,volume
unit = 1.60217662e3

def getnum(s):
	return findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def read_para_OUTCAR():
	global p_i,p_e,p_s,volume
	volume = float(getnum(getoutput("grep 'volume' OUTCAR|tail -1 "))[0])
	numbers = getnum(getoutput("grep 'dipole moment' OUTCAR"))
	p_i = [float(numbers[0]),float(numbers[1]),float(numbers[2])]
	p_s = [float(numbers[4]),float(numbers[5]),float(numbers[6])]
	p_e = [float(numbers[7]),float(numbers[8]),float(numbers[9])]

def print_dip():
	global p_i,p_e,volume
	p_tot = array(p_i)+array(p_e)
	p = sqrt(sum(p_tot**2))/volume
	print(p*unit)

read_para_OUTCAR()
print_dip()
