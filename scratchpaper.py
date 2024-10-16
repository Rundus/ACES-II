import matplotlib.pyplot as plt

fig, ax = plt.subplots(2,sharex=True)

xData = [1,2,3]
yData = [4,5,6]

ax[0].plot(xData,yData)
# ax[1].axis('off')

plt.subplots_adjust(hspace=0)
plt.show()