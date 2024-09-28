"""
Post process for the results from the excel file
(1) The strain energy curves for plate and total;
(2) The strain energy ratio for plate to the total;
(3) The X configuration for the specified frame;
(4) The Y Configuration for the specified frame;
(5) The Crease spring Moment results for the specified frame
(6) The Average Crease Moment during the unfolding process.
(7) The last configuration of Crease when the unfolding
"""
import pandas as pd
from math import *
import matplotlib.pyplot as plt

jobnumber = 14
fpath = "SingleCrease"+str(jobnumber)+".xlsx"
df =pd.read_excel(fpath)


"""
(1) The strain energy curves for plate and total;
(2) The strain energy ratio for plate to the total;
"""
time = df.loc[:, "time"]
PlateALLSE = df.loc[:, "PlateALLSE"]
TotalALLSE = df.loc[:, "TotalALLSE"]
CreaseALLSE = TotalALLSE - PlateALLSE
fig, ax = plt.subplots()
plt.xlabel('Time')
plt.ylabel('ALLSE')
line1, = ax.plot(time, PlateALLSE, linewidth=2.0)
line2, = ax.plot(time, TotalALLSE, linewidth=2.0)
line3, = ax.plot(time, CreaseALLSE, linewidth=2.0)
CreaseToPlate = CreaseALLSE /TotalALLSE
print("SingleCrease%s, The ratio of crease ALLSE to total ALLSE is %f" %(jobnumber, CreaseToPlate.mean()))
ax2 = ax.twinx()
line4, = ax2.plot(time, CreaseToPlate, color="red", linewidth=2.0)
plt.legend([line1, line2, line3, line4], ["PlateALLSE", "TotalALLSE", "CreaseALLSE", "CreaseToPlate"])
# ax2.legend()
# ax2.set(xlim=(0, 1), ylim=(, 1))
figname = "SingleCrease"+str(jobnumber)+"ALLSE"
fig.suptitle(figname)
plt.ylabel('Ratio')
fig.savefig(figname+'.jpg')
plt.show()

"""
(3) The X configuration for the specified frame;
Specified frame list Frames[] and the simulation process percent Prcess[]
"""
step = ceil(len(df)/5)
Frames = [x for x in range(0, len(df), step)]
Frames.append(len(df)-1)
Process = [float(x)/len(df) for x in Frames]

fig, ax = plt.subplots()
figname = "SingleCrease"+str(jobnumber)+"XCoord"
fig.suptitle(figname)
plt.xlabel('Coordinate X')
plt.ylabel('Coordinate Z')
for i in range(len(Frames)):
    XListX = list(filter(lambda x: "XCoordX" in x, df.columns.values.tolist()))
    XListZ = list(filter(lambda x: "XCoordZ" in x, df.columns.values.tolist()))
    dfXListX = df[XListX].loc[Frames[i], :]
    dfXListZ = df[XListZ].loc[Frames[i], :]
    label = str(floor(Process[i]*100))+"%"
    ax.plot(dfXListX, dfXListZ, linewidth=2.0, label=label)
ax.legend()
fig.savefig(figname+'.jpg')
plt.show()


"""
(4) The Y configuration for the specified frame;
"""
fig, ax = plt.subplots()
figname = "SingleCrease"+str(jobnumber)+"YCoord"
fig.suptitle(figname)
plt.xlabel('Coordinate Y')
plt.ylabel('Coordinate relative Z')
for i in range(len(Frames)):
    YListY = list(filter(lambda x: "YCoordY" in x, df.columns.values.tolist()))
    YListZ = list(filter(lambda x: "YCoordZ" in x, df.columns.values.tolist()))
    dfYListY = df[YListY].loc[Frames[i], :]
    dfYListZ = df[YListZ].loc[Frames[i], :]
    miniZ = dfYListZ.min()
    dfYListZ = dfYListZ -miniZ
    label = str(floor(Process[i]*100))+"%"
    ax.plot(dfYListY, dfYListZ, linewidth=2.0, label=label)
ax.legend()
fig.savefig(figname+'.jpg')
plt.show()


"""
(5) The Crease Spring Moment configuration for the specified frame;
"""
fig, ax = plt.subplots()
figname = "SingleCrease"+str(jobnumber)+"CreaseMoment"
fig.suptitle(figname)
plt.xlabel('Coordinate Y')
plt.ylabel('Crease relative Moment percent')
for i in range(len(Frames)):
    YListY = list(filter(lambda x: "YCoordY" in x, df.columns.values.tolist()))
    CMoment = list(filter(lambda x: "Element ASSEMBLY" in x, df.columns.values.tolist()))
    dfYListY = df[YListY].loc[Frames[i], :]
    dfCMoment = df[CMoment].loc[Frames[i], :]
    miniM = dfCMoment.min()
    dfCMoment = (dfCMoment - miniM)/abs(miniM)
    label = str(floor(Process[i]*100))+"%"
    ax.plot(dfYListY, dfCMoment, linewidth=2.0, label=label)
ax.legend()
fig.savefig(figname+'.jpg')
plt.show()

"""
(6) The Average Crease Moment during the unfolding process.
"""
time = df.loc[:, "time"]
CMoment = list(filter(lambda x: "Element ASSEMBLY" in x, df.columns.values.tolist()))
dfCMomentM = df[CMoment].mean(axis=1)
fig, ax = plt.subplots()
plt.xlabel('Time')
plt.ylabel('Average Crease Moment')
ax.plot(time, dfCMomentM, linewidth=2.0)
figname = "SingleCrease"+str(jobnumber)+"Average_Crease_Moment"
fig.suptitle(figname)
fig.savefig(figname+'.jpg')
plt.show()

"""
(7) The Average Crease Moment during the unfolding process.
"""

YListY = list(filter(lambda x: "YCoordY" in x, df.columns.values.tolist()))
YListZ = list(filter(lambda x: "YCoordZ" in x, df.columns.values.tolist()))
dfYListY = df[YListY].loc[df.index[-1], :]
dfYListZ = df[YListZ].loc[df.index[-1], :]
miniZ = dfYListZ.min()
dfYListZ = dfYListZ - miniZ
fig, ax = plt.subplots()
figname = "SingleCrease"+str(jobnumber)+"YCoordLast"
fig.suptitle(figname)
plt.xlabel('Coordinate Y')
plt.ylabel('Coordinate relative Z')
ax.plot(dfYListY, dfYListZ, linewidth=2.0)
fig.savefig(figname+'.jpg')
plt.show()