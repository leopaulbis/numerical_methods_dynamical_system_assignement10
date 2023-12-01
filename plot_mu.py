import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


file_path_xp_1 = 'results_I1_xp.csv'  # Replace 'your_file.csv' with your CSV file's path
data_xp_1 = pd.read_csv(file_path_xp_1,header=None)
xp_1=data_xp_1.iloc[:,0]

file_path_xp_2="results_I2_xp.csv"
data_xp_2=pd.read_csv(file_path_xp_2,header=None)
xp_2=data_xp_2.iloc[:,0]

file_path_xp_3="results_I3_xp.csv"
data_xp_3=pd.read_csv(file_path_xp_3,header=None)
xp_3=data_xp_3.iloc[:,0]


mu1= np.linspace(0.001, 0.015, round((0.015 - 0.001) / 0.00001) + 1)
mu2= np.linspace(0.015, 0.05, round((0.05 - 0.015) / 0.0001) + 1)
mu3= np.linspace(0.05, 0.49, round((0.49 - 0.05) / 0.001) + 1)

mu=np.concatenate((mu1,mu2),axis=0)
mu=np.concatenate((mu,mu3),axis=0)
xp=np.concatenate((xp_1,xp_2),axis=0)
xp=np.concatenate((xp,xp_3),axis=0)
print(np.max(xp))

print(np.shape(mu))
plt.scatter(mu,xp,s=1)
plt.show()


mu_bis=[]
for i in range(len(mu)):
    if mu[i]<0.1:
        mu_bis.append(mu[i])

size=len(mu_bis)
xp_bis=xp[:size]

plt.scatter(mu_bis,xp_bis,s=2)
plt.show()
