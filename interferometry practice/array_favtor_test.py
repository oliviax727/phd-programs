# pylint: skip-file
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:03:28 2022

@author: Aishwarya Selvaraj

"""
import matplotlib.colors as colors
import matplotlib.pylab as plt
# from scipy import interpolate
from matplotlib import cm
import numpy as np
# import xlrd

#%%

theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
# d = lamda/2
l = np.sin(theta)*np.cos(phi)
m = np.sin(theta)*np.sin(phi)

gain_ns = np.sqrt(1- l*l)
gain_ew = np.sqrt(1- m*m)

# Considering the element pattern: 
h = 0.3/lamda
ground_plane_effect = 2*np.sin(2*h*np.pi*np.cos(theta))

n = 8
distance_dipoles = 1.1
drange = ((n-1)*distance_dipoles)/2.
d = np.linspace(-drange, drange, n)
#%% With all dipoles active voltage response for MWA. 

mwa_af = np.zeros((np.shape(theta)), dtype = complex)
plt.figure()
count = 0 
for p in range(n):
    for q in range(n):
        mwa_af+= np.exp(1j*beta*(l*d[p]+m*d[q]))
        plt.plot(d[q], d[p], '*r')
        plt.text(d[q], d[p], count, bbox=dict(facecolor='red', alpha=0.5), fontweight ='bold',fontsize=12 )
        count+= 1
        # print(d[q], d[p])
res = (mwa_af*ground_plane_effect)  
res_xx = gain_ew*res
res_yy = gain_ns*res

power_yy = np.real(res_yy*np.conj(res_yy))
power_xx = np.real(res_xx*np.conj(res_xx))

xx_norm = power_xx/np.max(power_xx)
yy_norm = power_yy/np.max(power_yy)

# if n ==8:
#     plt.title("Tile layout for CRAM")
# elif n == 4:
#     plt.title("Tile layout for MWA")
# else:
#     pass

# plt.figure()
# cs = plt.pcolormesh(l, m, mwa_active, norm=colors.SymLogNorm(linthresh=0.3, linscale= 0.3 ,vmin=np.nanmin(mwa_active), vmax= np.nanmax(mwa_active), base=10), cmap='jet_r', shading='auto') 
# bar = plt.colorbar(cs, ticks=[0, 1, 10, 100, 10**2], orientation='vertical', extend='both', label=r'$\log T\ [\mathrm{K}]$')
# plt.title('voltage response of MWA')
# plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
# plt.ylabel('directional cosine m', fontsize=12, fontweight ='bold')



# 3D plot of voltage response
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(l, m, yy_norm, rstride=1, cstride=1, antialiased=True)
plt.show()
print(np.sum(yy_norm))
#%% With dead dipoles, voltage response for MWA. 
mwa_af = np.zeros((np.shape(theta)), dtype = complex)
 
dead_dipole = [53, 44, 57, 33, 16, 26]

count = 0 
plt.figure()
for p in range(n):
    for q in range(n):
        condition = count in dead_dipole
        if condition == True:
            # mwa_af+= np.ones((np.shape(mwa_af)))*(-1)
            # mwa_af-= np.exp(1j*beta*(l*d[q]+m*d[p]))
            print(count)
            plt.plot(d[q], d[p], '*k')
            plt.text(d[q], d[p], count, bbox=dict(facecolor='black', alpha=0.5), fontweight ='bold',fontsize=12 )
           
        else:
            mwa_af+= np.exp(1j*beta*(l*d[q]+m*d[p]))
            plt.plot(d[q], d[p], '*r')
            plt.text(d[q], d[p], count, bbox=dict(facecolor='red', alpha=0.5), fontweight ='bold',fontsize=12 )
        count += 1

res = (mwa_af*ground_plane_effect)  
res_xx = gain_ew*res
res_yy = gain_ns*res

power_yy = np.real(res_yy*np.conj(res_yy))
power_xx = np.real(res_xx*np.conj(res_xx))

xx_norm = power_xx/np.max(power_xx)
yy_norm = power_yy/np.max(power_yy)

# plt.figure()
# cs = plt.pcolormesh(l, m, mwa_dead, norm=colors.SymLogNorm(linthresh=0.3, linscale= 0.3 ,vmin=np.nanmin(mwa_dead), vmax= np.nanmax(mwa_dead), base=10), cmap='jet_r', shading='auto') 
# bar = plt.colorbar(cs, ticks=[0, 1, 10, 100, 10**2], orientation='vertical', extend='both', label=r'$\log T\ [\mathrm{K}]$')
# plt.title('voltage response of MWA')
# plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
# plt.ylabel('directional cosine m', fontsize=12, fontweight ='bold')



# 3D plot of voltage response without the ground plane effect. 
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(l, m, yy_norm, rstride=1, cstride=1, antialiased=True)
plt.show()
print(np.sum(yy_norm))



# print((np.sum(mwa_active) - np.sum(mwa_dead))/np.sum(mwa_active))
#%%
# # voltage response for CRAM. 
# n = 8
# afx = 0; afy = 0
# drange = ((n-1)*(lamda/2))/2.
# d = np.linspace(-drange, drange, n)
# # d = np.linspace(0, (n-1)*lamda/2, n)
# cram_af = 0
# for p in range(n):
#     for q in range(n):
#         cram_af = cram_af + np.exp(1j*beta*(l*d[q]+m*d[p]))
# cram_res = abs(cram_af*ground_plane_effect)



# plt.figure()
# cs = plt.pcolormesh(l, m, cram_res, norm=colors.SymLogNorm(linthresh=0.3, linscale=0.3 ,vmin=np.nanmin(cram_res), vmax= np.nanmax(cram_res), base=10), cmap='jet_r', shading='auto') 
# bar = plt.colorbar(cs, ticks=[0, 1, 10, 100, 10**2], orientation='vertical', extend='both')
# # plt.title('voltage response of CRAM')
# plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
# plt.ylabel('directional cosine m',  fontsize=12, fontweight ='bold')
# plt.show()


#%%
dead_dipole = [57, 49, 41, 33, 16, 26, 25, 60, 45, 53]
dead_dipole = np.arange(5)

biast_x = []
biast_y = []
mwa_afx = np.zeros((np.shape(theta)), dtype = complex)
mwa_afy = np.zeros((np.shape(theta)), dtype = complex)
count = 0 
for p in range(n):
    for q in range(n):
        if count in dead_dipole:
            mwa_afx+= 0
            mwa_afy+= 0            
        elif count in biast_x:
            mwa_afx+= 0
            mwa_afy+= np.exp(1j*beta*(l*d[q]+m*d[p]))
        elif count in biast_y:
            mwa_afx+= np.exp(1j*beta*(l*d[q]+m*d[p]))
            mwa_afy+= 0
            
        else:
            mwa_afx+= np.exp(1j*beta*(l*d[q]+m*d[p]))
            mwa_afy+= np.exp(1j*beta*(l*d[q]+m*d[p]))
        count += 1


res_xx = gain_ew*(mwa_afx*ground_plane_effect)  
res_yy = gain_ns*(mwa_afy*ground_plane_effect)  

power_yy = np.real(res_yy*np.conj(res_yy))
power_xx = np.real(res_xx*np.conj(res_xx))

xx_norm = power_xx/np.max(power_xx)
yy_norm = power_yy/np.max(power_yy)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(l, m, yy_norm, rstride=1, cstride=1, antialiased=True)
plt.show()
print(np.sum(yy_norm))

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(l, m, xx_norm, rstride=1, cstride=1, antialiased=True)
# plt.show()
