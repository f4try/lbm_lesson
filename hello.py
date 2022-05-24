import numpy as np
import matplotlib.pyplot as plt 
import math
x=np.linspace(-2,2.001,100)
print(x)
y=x**3
y2=x**2
y3=x
plt.scatter(x,y,"r--")
plt.plot(x,y2,"b--")
plt.plot(x,y3,"g--")
plt.grid(True)
plt.show()