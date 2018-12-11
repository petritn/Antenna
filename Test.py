from numpy import *  # for outer and a range
import pylab as p  # for figure
from mpl_toolkits.mplot3d import Axes3D

fig = p.figure()
ax = Axes3D(fig)

theta = arange(0,pi,pi/10)
phi = arange(0,2*pi,pi/10)
r = 2 * pow(math.e, -((theta**4)/(0.25**2)))  # need to distort the radius by some function

x = r*outer(cos(phi), sin(theta))
y = r*outer(sin(phi), sin(theta))
z = r*outer(ones(phi.shape), cos(theta))


ax.plot_wireframe(x,y,z)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

p.show()