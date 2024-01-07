import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
from numpy import trapz
from scipy.integrate import cumtrapz
import math

#Updates the x and y values of the SFD for edgecase where a train wheel is directly on the support
def edgecase(x,y):
  #support A
  if x.count(0)>2:
    x.remove(0)
    y_points.pop(1)
    y_points[1]-=P_train[1]
  #support B
  elif x.count(1200)>=2:
    x.pop()
    y_points.pop(-2)
    y_points[-2]-=P_train[-2]
  return x,y

#Cross section properties
height = 140
ybar = 100.4
ybot = ybar
ytop = height - ybar
pi = math.pi
I = 577619
Qcent = 15335
Qglue = 9743.4
base = 2.54
E = 4000
mu = 0.2
#Max stresses and strains from material properties
S_tens = 30
S_comp = 6
T_max = 4
T_gluemax = 2

# Initial values
total_area = 813 * 1016 * 1.27

n = 1200 #analyzed for each mm along the length
P = 400  # weight of train

x_train = [52, 228, 392, 568, 732, 908] #distance along the train
#P_train = [400 / 6] * 6  #Load case 1
P_train=[90,90,66.65,66.65,66.65,66.65] #Load case 2

#Values to store maximum forces
max_shear = 0
max_moment = 0
shear_distance = 0
moment_distance = 0
shearatm = 0 #shear at max moment
momentats = 0 #moment at max shear

# Create plot
plt.figure(figsize=(8, 6))
plt.axhline(y=0, color='red', linestyle='--', label='y = 0') #Draw an x-axis

#Ay and By at max moment
Aym = 0
Bym = 0
#Ay and By at max shear
Ays = 0
Bys = 0

#edge cases where a wheel of the train is directly on one of the supports
badx = [52,228,392,568,732,908,1148,972,808,632,468,292,1252,1428,1592,1768,1932,2108]
#Which graph to display
bmd = False
sfd = False

while not bmd and not sfd:
  graph = input("Which graph would you like? (BMD or SFD?): ")
  if graph.lower() == 'bmd':
    bmd = True
  elif graph.lower()=='sfd':
    sfd = True
shear = []
# For every position on the bridge starting at 52
for i in range(52, 2109):
  sum_momentsA = 0
  force_y = 0
  #graph points
  x_points = []
  y_points = []

  # For every wheel force
  for a in range(6):
    #if the wheel is on the bridge
    distance = i-x_train[a]
    if distance >=0 and distance<=1200:
      force_y -= P_train[a] #add to sum of forces in y-direction
      sum_momentsA -= P_train[a] * (i - x_train[a]) #add to sum of the moments at A
      x_points.append(i - x_train[a]) #add distance to x-points

  #calculate By and Ay
  B_y = abs(sum_momentsA / 1200)
  A_y = abs(force_y + B_y)
  y = A_y

  #Check if any wheels are directly on the supports
  # For every wheel force
  for a in range(6):
    #if the wheel is on the bridge
    distance = i-x_train[a]
    if distance >=0 and distance<=1200:
      y -= P_train[-a - 1]  # update the current y position
      y_points.append(y)  # add to y-points
      #y -= P_train[-a - 1] #update the current y position
      #y_points.append(y) #add to y-points

  #Get the graph points
  num_points = len(y_points)
  x_points.reverse()
  #Insert start point
  x_points.insert(0, 0)
  y_points.insert(0, A_y)
  #Insert end point
  x_points.append(1200)
  y_points.append(y_points[-1] + B_y)
  #Graph starts at 0
  x_points.insert(0,0)
  y_points.insert(0,0)

  #Calculate area under SFD for BMD
  areas = [0] * (len(x_points) - 1) #Create array for areas
  for j in range(0, len(x_points) - 2):
    x_dis = x_points[j + 1] - x_points[j] #calculate distance between two forces
    y_dis = max(y_points[j + 1], y_points[j]) #take the max of the two shear
    #Calculate the area and add it to array
    if j > 0:
      areas[j] = areas[j - 1] + (-x_dis * y_dis)
    else:
      areas[0] = -x_dis * y_dis

  #Check for maxes
  potential_max = abs(max(y_points, key=abs))
  potential_moment = abs(max(areas, key=abs))
  if potential_max > max_shear: #Check if this is the maximum shear so far
    #Update values
    max_shear = potential_max
    shear_distance = i
    Ays = A_y
    Bys = B_y
    momentats = potential_moment
  if potential_moment > max_moment: #Check if this is the maximum moment so far
    #Update values
    max_moment = potential_moment
    moment_distance = i
    shearatm = potential_max
    Aym = A_y
    Bym = B_y
  if bmd:
    x_points.remove(0)
    #Connect the points with horizontal lines
    plt.plot(x_points,areas,linestyle='-',linewidth=1)
  else:
    if i in badx: #Checks for edge case
      x_points, y_points = edgecase(x_points, y_points)
    #Connect the points with slopes
    plt.step(x_points, y_points, linestyle='-', linewidth=1, where='post')

#Display Calculated Values
print(f"Ay at max moment = {Aym}")
print(f"By at max moment= {Bym}")
print(f"Ay at max shear = {Ays}")
print(f"By at max shear = {Bys}")
print("\nMax Forces: ")
print(f"Max shear: {max_shear} at {shear_distance}")
print(f"Max Moment: {max_moment} at {moment_distance}")
print(f"Shear at max moment: {shearatm}")
print(f"Moment at max shear:{momentats}")

#SFD Graph
# Customize the plot
if bmd:
  plt.title('BMD')
  plt.xlabel('Distance')
  plt.ylabel('Moment')
else:
  plt.title('SFD')
  plt.xlabel('Distance')
  plt.ylabel('Force')
#plt.step(x_points, y_points, linestyle='-', linewidth=1, label='SFD')  # Add this line
plt.legend()
plt.grid(True)

# Show the plot
plt.show()


#-------------------------------------------------------------------------
#Calculations

#Calculate max compresive and tensile stress
S_top = max_moment * ytop / I
S_bot = max_moment * ybot / I
#Calculate max shear
T_cent = (max_shear * Qcent) / (I * base)
base = 20
T_glue = (max_shear * Qglue) / (I * base)

#Check buckling
#Case 1
t = 2.54
b = 85
S_buck1 = ((4 * (pi**2) * E) / (12 * (1 - mu**2))) * (t / b)**2
#Case 2
t = 2.54
b = 22.54
S_buck2 = (0.425 * (pi**2) * E) * (t / b)**2 / (12 * (1 - mu**2))
#Case 3
t = 1.27
b = 37.06
S_buck3 = (6 * (pi**2) * E) * (t / b)**2 / (12 * (1 - mu**2))
#Case 4
t = 1.27
height = 37.06
V_buck = (5 * (pi**2) * E) * ((t / height)**2) / (12 - (1 - mu)**2)

#Shear buckling in web
t = 1.27
height = 134.92
shearBuck = (5 * (pi**2) * E) * ((t / height)**2) / (12 - (1 - mu)**2)

print("\nFOS:")
#Calculate factor of safety
FOS_tens = S_tens / S_bot
FOS_comp = S_comp / S_top
FOS_shear = T_max / T_cent
FOS_glue = T_gluemax / T_glue
FOS_buck1 = S_buck1 / S_top
FOS_buck2 = S_buck2 / S_top
FOS_buck3 = S_buck3 / S_top
FOS_buckV = V_buck / T_cent
#Dictionary of all FOS
fos = {'tension':FOS_tens,
      'compression':FOS_comp,
      'shear': FOS_shear,
      'glue': FOS_glue,
      'buck1': FOS_buck1,
      'buck2':FOS_buck2,
      'buck3':FOS_buck3,
      'buckV':FOS_buckV}
#Get minimum FOS
minFOS_key = min(fos, key=fos.get)
minFOS_value = fos[minFOS_key]
print(
    f"tens: {FOS_tens}\ncomp:{FOS_comp}\nshear:{FOS_shear}\nglue:{FOS_glue}\nbuck1:{FOS_buck1}\nbuck2:{FOS_buck2}\nbuck3:{FOS_buck3}\nshear buck:{FOS_buckV}"
)

print(f"\nMinimum FOS: {minFOS_key} with value of {minFOS_value}")

#Calculate max forces
Mf_tens = FOS_tens * max_moment
Mf_comp = FOS_comp * max_moment
Vf_shear = FOS_shear * max_shear
Vf_glue = FOS_glue * max_shear
Mf_buck1 = FOS_buck1 * max_moment
Mf_buck2 = FOS_buck2 * max_moment
Mf_buck3 = FOS_buck3 * max_moment
Vf_buckV = FOS_buckV * max_shear
print("\nMax forces:")
print(
    f"Max tension: {Mf_tens}\nMax compression: {Mf_comp}\nMax shear: {Vf_shear}\nMax glued: {Vf_glue}\nMax buckling1:{Mf_buck1}\nMax buckling2:{Mf_buck2}\nMax buckling3:{Mf_buck3}\nMax Shear Buckling:{Vf_buckV}"
)
