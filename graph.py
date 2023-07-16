import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

df = pd.read_csv("scatterplot.csv")

kde = stats.gaussian_kde(df['x'], bw_method='scott')

x = df['x'].values
xmin, xmax = x.min(), x.max()
xgrid = np.linspace(xmin, xmax, 100)
plt.plot(xgrid, kde(xgrid))

plt.xlabel('x')
plt.ylabel('Density')
plt.title('Kernel Density Estimation')
plt.show()
