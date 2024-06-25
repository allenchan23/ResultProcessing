import math
from rocketcea.cea_obj import CEA_Obj,add_new_fuel
import numpy as np

card_str = """
fuel C2H4(G) C 2.0 H 4.0 wt%=100.0
h,cal=12532.24962   t(k)=298.15  rho=0.00118
"""
add_new_fuel('C2H4', card_str)
C = CEA_Obj( oxName='GOX', fuelName='C2H4')
C.get_Cstar()

MR = C.getMRforER(ERphi = 1.0)
print(MR)

print(C.get_Cstar(Pc=86*0.1450, MR=MR)*0.3048)
cstar = C.get_Cstar(Pc= 86*0.1450, MR=MR)*0.3048
mdot = 0.02
pc = 86*10**3
A_t = cstar* mdot/pc
print(A_t)
print(np.sqrt(A_t/np.pi)*1000*2)
