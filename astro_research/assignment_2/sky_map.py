"""
Created by: Kristian Ivanov
URN: 6534278

The script only works if the MiklyWaySats.txt is in the same directory as the
script.

(Folder directory example)


sky_map.py       MilkyWaySats.txt


Note:
For some reason this seems to work better when executed with >>>python sky_map.py
"""
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Define custom column names
colnames=['Name', 'RA', 'Dec', 'Mv','vel']

#Extract the txt data (separation is tab (\t))
sky_ext = pd.read_csv('MilkyWaySats.txt',header=None,sep='\t',skiprows=1,\
                        names=colnames)
#Magnitude limit
mv_lim = -7.5


#Make astropy skycoords with the columns of data
sky_data = SkyCoord(ra=np.array(sky_ext['RA']),dec=np.array(sky_ext['Dec']),unit='deg')
sky_data_galact = sky_data.galactic #galactic conversion


#Plotting
fig = plt.figure()
equat = plt.subplot(211, projection="aitoff")
plt.grid(True)
plt.title("Equatorial Coordinates")
plt.xlabel('Right Ascension [degrees]')
plt.ylabel('Declination [degrees]')
galact = plt.subplot(212, projection="aitoff")
plt.grid(True)
plt.title('Galactic Coordinates')
plt.xlabel('Longtitude')
plt.ylabel('Latitude')

#logic gate for color coding
for row in range(len(sky_ext)):
    if   -500<=sky_ext['vel'][row]<-250:
        color='red'
    elif -250<=sky_ext['vel'][row]<0:
        color='orange'
    elif    0<=sky_ext['vel'][row]<250:
        color='green'
    elif  250<=sky_ext['vel'][row]<500:
        color='blue'
    if sky_ext['Mv'][row] < mv_lim:
        equat.scatter(sky_data.ra.value[row],sky_data.dec.value[row],marker='*',c=color)
        galact.scatter(sky_data_galact.l.value[row],sky_data_galact.b.value[row],marker='*',c=color)
    elif sky_ext['Mv'][row] > mv_lim:
        equat.scatter(sky_data.ra.value[row],sky_data.dec.value[row],marker='^',c=color)
        galact.scatter(sky_data_galact.l.value[row],sky_data_galact.b.value[row],marker='^',c=color)

#A rather clunky but still functional way to make a custom legend and colorcode points
legend_mk1=plt.scatter(0,0,marker='^',label='Ultra-faint',c='k')
legend_mk2=plt.scatter(0,0,marker='*',label='Classic',c='k')
legend_mk3=plt.scatter(0,0,marker='o',label='-500 km/s < vel < -250 km/s',c='red')
legend_mk4=plt.scatter(0,0,marker='o',label='-250 km/s < vel < 0 km/s',c='orange')
legend_mk5=plt.scatter(0,0,marker='o',label='0 km/s < vel < 250 km/s',c='green')
legend_mk6=plt.scatter(0,0,marker='o',label='250 km/s < vel < 500 km/s',c='blue')
fig.legend()
legend_mk1.remove()
legend_mk2.remove()
legend_mk3.remove()
legend_mk4.remove()
legend_mk5.remove()
legend_mk6.remove()

plt.show()
