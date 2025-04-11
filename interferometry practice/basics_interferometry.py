# pylint: skip-file
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:03:28 2022

@author: Aishwarya Selvaraj
This is the script written to simulate primary beams for MWA/CRAM tiles. 
It includes 1D array factor. 
"""

from __future__ import print_function, division
import os
import matplotlib.colors as colors
import matplotlib.pylab as plt
from scipy import interpolate
from matplotlib import cm
import numpy as np
import xlrd



save_results_to = './'
#%% equi-distant antenna sources 
'''
This code generates the array factor for n number of equidistant sources 
which are placed at a distance of lamda/2 from each other. These isotropic
antennas are assumed to be placed on the z-axis .
'''
delta = 0 # phase difference of one source with respect to another. 
angle = np.arange(0, 180) # its covers the horizon. 
theta = np.deg2rad(angle) 
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
d = lamda/2
l = np.cos(theta)
chi = beta*d*l  + delta # total phase difference.
n = 8
af = (np.sin(n*chi/2))/(n*np.sin(chi/2)) #normalised array factor. 
# plt.figure()
plt.plot(angle, af)
plt.xlabel('angle in degrees', 
               fontweight ='bold')
plt.ylabel('normalised amplitude', 
               fontweight ='bold')

# plt.figure()
# plt.polar(theta, af)


#%%
f = (150) 
lamda = ((300)/f)
beta = (2*np.pi/ lamda)

angle = np.linspace(-90, 90+1, 513)
theta = np.deg2rad(angle)
h = 0.3/lamda #standard value in MWA papers. 
element_pattern = np.sin(2*h*np.pi*np.cos(theta))#ground plane effect



angle = np.linspace(0, 180+1, 513) # its covers the horizon. 
theta = np.deg2rad(angle) 




n = 4
drange = ((n-1)*1.1)/2.
d = np.linspace(-drange, drange, n)
l = np.cos(theta)
mwa_af = 0;
for k in range(n):
    mwa_af += np.exp((1j*beta*l*d[k]))
    
mwa_res = mwa_af*element_pattern
    

# plt.plot(angle, mwa_res)
# plt.xlabel('angle in degrees', 
#                fontweight ='bold')
# plt.ylabel('normalised amplitude', 
#                fontweight ='bold')


n = 8
drange = ((n-1)*1.1)/2.
d = np.linspace(-drange, drange, n)
l = np.cos(theta)
cram_af = 0;
for k in range(n):
    cram_af += np.exp((1j*beta*l*d[k]))
    
cram_res = cram_af*element_pattern
    

# plt.plot(angle, cram_res)
# plt.xlabel('angle in degrees', 
#                fontweight ='bold')
# plt.ylabel('normalised amplitude', 
#                fontweight ='bold')

plt.figure()
power1 = np.real(mwa_res*np.conj(mwa_res))
power1/= np.max(power1)
plt.plot(angle, power1, '-k')
plt.xlabel('angle in degrees', 
               fontweight ='bold')
plt.ylabel('normalised amplitude', 
               fontweight ='bold')


power2 = np.real(mwa_res*np.conj(cram_res))
power2/= np.max(power2)
plt.plot(angle, power2, 'r')
plt.xlabel('angle in degrees', 
               fontweight ='bold')
plt.ylabel('normalised amplitude', 
               fontweight ='bold')


power3 = np.real(cram_res*np.conj(cram_res))
power3/= np.max(power3)
plt.plot(angle, power3, 'b')
plt.xlabel('angle in degrees', 
               fontweight ='bold')
plt.ylabel('normalised amplitude', 
               fontweight ='bold')


plt.legend(['MWA-MWA', 'MWA-CRAM', 'CRAM-CRAM'])


from scipy.interpolate import UnivariateSpline
from scipy import integrate
import pylab as pl
signal = power1
spline = UnivariateSpline(angle, signal - np.max(signal)/2, s =0 )
r1, r2 = spline.roots()
FWHM = (r2-r1)
sigma = np.deg2rad(FWHM)/2.3548
print(FWHM, 1/(sigma))
area = integrate.simps(signal, angle )
print(np.deg2rad(area))
plt.plot(angle, signal)
pl.axvspan(r1, r2, facecolor = 'g', alpha = 0.5)
pl.show()

# plt.figure()
# fft_cram_mwa = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(power2)))
# plt.plot(np.real(fft_cram_mwa))
#%% for different values on n.
'''
This is section, 1D array factor is calculated for different values of n, i.e. 
for different number of isotropic antenna sources. But, they are considered 
to be kept equidistant from each other. 
'''
delta = 0 # phase difference of one source with respect to another. 
angle = np.arange(0, 180) # its covers the horizon. 
theta = np.deg2rad(angle) 
f = (150*10**6) 
lamda = (3*10**8)/f
beta = (2*np.pi/ lamda)
d = lamda/2
chi = beta*d*np.cos(theta) + delta # total phase difference.
fig, axs = plt.subplots(2, 4)
plt.xlim(angle[0], angle[-1]+1)
N = [2, 4, 8, 10, 16, 25, 50, 100]
count = 0
for row in range(2):
    for col in range(4):
        n = N[count]
        af = (np.sin(n*chi/2))/(n*np.sin(chi/2))
        axs[row, col].plot(angle, af)
        axs[row, col].legend(["$" + str(n) + "$"], loc ="best")
        count = count+1
        axs[row, col].set_xlabel('angle in degrees', 
               fontweight ='bold')
        axs[row, col].set_ylabel('normalised amplitude', 
               fontweight ='bold')




#%% When the antennas are placed at different lengths from each other.
'''
In this section, 1D array factor is calulated for n = 5 sources which are placed
at different distance from each other. 
''' 
delta = 0 
angle = np.arange(0, 180)
theta = np.deg2rad(angle)
f = (150*10**6)
lamda = (3*10**8)/f
beta = (2*np.pi/ lamda)
n = 5

d1 = lamda/2
chi1 = beta*d1*np.cos(theta) + delta 

d2 = lamda/4
chi2 = beta*d2*np.cos(theta) + delta 

d3 = lamda/2
chi3 = beta*d3*np.cos(theta) + delta 

d4 = lamda/4
chi4 = beta*d4*np.cos(theta) + delta 


af = (1 + np.exp(1j*chi1) + np.exp(2j*chi2) + np.exp(3j*chi3) + np.exp(4j*chi4))/n
plt.figure(2)
plt.plot(angle, np.real(af))
plt.plot(angle, np.imag(af))
plt.legend(['real terms', 'imaginary terms'],  loc ="top right")
plt.xlabel('angle in degrees', 
               fontweight ='bold')
plt.ylabel('normalised amplitude', 
               fontweight ='bold')
#%%
'''
This is the same code as above but is automated by writing inside a for loop. 
'''
delta = 0
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
# theta = np.deg2rad(x) 
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
d = lamda/2 
n = 4

l = np.sin(theta)*np.cos(phi)
m = np.sin(phi)*np.sin(theta)

chix = beta*d*l + delta # total phase difference.
chiy = beta*d*m  + delta # total phase difference.

afx = 0; afy = 0

for k in range(n):
    afx = afx + np.exp((1j*chix*k))
    afy = afy + np.exp((1j*chiy*k))
    
af = afx*afy


plt.imshow(np.imag(afx))
# plt.xlabel('angle in degrees', 
#                fontweight ='bold')
# plt.ylabel('normalised amplitude', 
#                fontweight ='bold')

#%% Ground plane effect. 
'''
This sections consists of the code that plots the element pattern of a dipole
that is nothing but the ground plane effect. 
'''
angle = np.arange(-90, 90+1) 
theta = np.deg2rad(angle)
f = [80*10**6, 90*10**6, 111*10**6, 150*10**6, 180*10**6, 220*10**6]
# f = [100*10**6, 150*10**6, 200*10**6, 300*10**6, 350*10**6]
plt.figure()
for k in range(len(f)):
    lamda = ((3*10**8)/f[k])
    h = 0.3/lamda
    l = np.cos(theta)
    pf = np.sin(2*h*np.pi*l)
    # total_voltage = af*pf
    plt.plot(angle, pf, label = str(f[k]/10**6) + "MHz")
    plt.legend()

plt.xlabel('angle in degrees', 
               fontweight ='bold')
plt.ylabel('normalised amplitude', 
               fontweight ='bold')                                                                                    


#%% 28th March
'''
 Considering the ground plane effect in 1D
'''



delta = 0 # phase difference of one source with respect to another. 
angle = np.arange(0, 180+1) 
theta = np.deg2rad(angle) 
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
d = lamda/2
l = np.cos(theta) 
chi = beta*d*l + delta # total phase difference.
n = 8
af = (np.sin(n*chi/2))/(n*np.sin(chi/2)) #normalised array factor. 


angle = np.arange(-90, 90+1)
theta = np.deg2rad(angle)
h = 0.3/lamda #standard value in MWA papers. 
element_pattern = np.sin(2*h*np.pi*np.cos(theta))#ground plane effect
# plt.figure(1)
# plt.plot(angle, element_pattern)


total_voltage = af*element_pattern

plt.figure(2)
plt.plot(l, element_pattern, 'r')
plt.plot(l, af, 'b') 
plt.plot(l, total_voltage, 'k')
plt.xlabel('cos\u03B8', fontweight ='bold',  fontsize=12)
plt.ylabel('normalised amplitude', fontweight ='bold',  fontsize=12)
plt.legend(['ground plane effect', 'array factor', 'total voltage'], loc ="best")




# def lin_interp(x, y, i, half):
#     return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

# def half_max_x(x, y):
#     half = max(y)/2.0
#     signs = np.sign(np.add(y, -half))
#     zero_crossings = (signs[0:-2] != signs[1:-1])
#     zero_crossings_i = np.where(zero_crossings)[0]
#     return [lin_interp(x, y, zero_crossings_i[0], half),
#             lin_interp(x, y, zero_crossings_i[1], half)]

# # make some fake data


# # find the two crossing points
# hmx = half_max_x(l, total_voltage)

# # print the answer
# fwhm = hmx[1] - hmx[0]
# print("FWHM:{:.3f}".format(fwhm))


#%% 2D representation of array factor
'''
When the linear array of isotropic antennas are placed on the xy-plane.
on xy-plane
theta = [0, 90] degrees
phi = [0, 360] degrees
d = lamda/2
n = 8
considering ground plane effect, where h = 0.3/lamda
'''

delta = 0 # phase difference of one source with respect to another. 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
# theta = np.deg2rad(x) 
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
d = lamda/2
n = 4

l = np.sin(theta)*np.cos(phi)
m = np.sin(phi)*np.sin(theta)

chix = beta*d*l + delta # total phase difference.
chiy = beta*d*m  + delta # total phase difference.

afx = 0; afy = 0

for k in range(n):
    afx = afx + np.exp((1j*chix*k))
    afy = afy + np.exp((1j*chiy*k))

af = abs(afx*afy)
# Considering the element pattern: 
h = 0.3/lamda
ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))

total_voltage = af*ground_plane_effect


# plt.figure()
# plt.imshow(total_voltage,origin = 'lower',  vmin = np.min(total_voltage), vmax = np.max(total_voltage))


# 3D plot of voltage response with the ground plane effect. 
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(l, m, total_voltage, rstride=1, cstride=1, antialiased=True)
plt.title('voltage response with ground plane effect')
plt.show()


# 3D plot of voltage response without the ground plane effect. 
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(l, m, af, rstride=1, cstride=1, antialiased=True)
plt.title('voltage response without ground plane effect')
plt.show()

plt.figure()
plt.pcolormesh(l, m, total_voltage) 
plt.show()

#%% 2D representation of array factor
'''
When the linear array of isotropic antennas are placed on the xz-axis.

'''

delta = 0 # phase difference of one source with respect to another. 
th = np.deg2rad(np.arange(0, 90))
ph = np.deg2rad(np.arange(0, 360))
[theta, phi] = np.meshgrid(th, ph)
# theta = np.deg2rad(x) 
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
d = lamda/2
n = 8
l = np.sin(theta)*np.cos(phi)
m = np.cos(theta) 
chix = beta*d*l  + delta # total phase difference.\
chiz = beta*d*m + delta # total phase difference.
afx = 0; afz = 0
for k in range(n):
    afx = afx + np.exp((1j*chix*k))
    afz = afz + np.exp((1j*chiz*k))

af = abs(afx*afz)
h = 0.3
ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))
total_voltage = af*ground_plane_effect

plt.figure()
plt.imshow(total_voltage, vmin = np.min(total_voltage), vmax = np.max(total_voltage))

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(l, m, total_voltage, rstride=1, cstride=1, antialiased=True)
plt.show()



#%% 2D representation of array factor
'''
When the linear array of isotropic antennas are placed on the yz-axis.

'''

delta = 0 # phase difference of one source with respect to another. 
th = np.deg2rad(np.arange(0, 180))
ph = np.deg2rad(np.arange(-90, 90))
[theta, phi] = np.meshgrid(th, ph)
# theta = np.deg2rad(x) 
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
d = lamda/2
n = 8

chiy = beta*d*np.sin(phi)*np.sin(theta)  + delta # total phase difference.\
chiz = beta*d*np.cos(theta)  + delta # total phase difference.
afy = 0; afz = 0
for k in range(n):
    afy = afy + np.exp((1j*chiy*k))
    afz = afz + np.exp((1j*chiz*k))

af = abs(afy*afz)
plt.figure()
plt.imshow(af, vmin = np.min(af), vmax = np.max(af))

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(theta, phi, af, rstride=1, cstride=1, antialiased=True)
plt.show()





#%% 2D representation of array factor
'''
When the linear array of isotropic antennas are placed on the xy-plane.
 - on xy-plane
 - theta = [0, 90] degrees
 - phi = [0, 360] degrees
 - d = lamda/2
 - n = 4
 - considering ground plane effect, where h = 0.3

'''

delta = 0 # phase difference of one source with respect to another. 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
# theta = np.deg2rad(x) 
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
d = lamda/2
n = 4

l = np.sin(theta)*np.cos(phi)
m = np.sin(phi)*np.sin(theta)

chix = beta*d*l + delta # total phase difference.
chiy = beta*d*m  + delta # total phase difference.

afx = 0; afy = 0

for k in range(n):
    afx = afx + np.exp((1j*chix*k))
    afy = afy + np.exp((1j*chiy*k))

af = abs(afx*afy)

# Considering the element pattern: 
h = 0.3/lamda
ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))

total_voltage = af*ground_plane_effect


# plt.figure()
# plt.imshow(total_voltage,origin = 'lower',  vmin = np.min(total_voltage), vmax = np.max(total_voltage))


# 3D plot of voltage response with the ground plane effect. 
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(l, m, total_voltage, rstride=1, cstride=1, antialiased=True)
plt.show()


# 3D plot of voltage response without the ground plane effect. 
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(l, m, af, rstride=1, cstride=1, antialiased=True)
# plt.show()

plt.figure()
plt.pcolormesh(l, m, total_voltage) 
plt.title("Polar plot of voltage response: MWA")
plt.show()

#%% Primary beam response. 
'''
This code is written only to create a subplot for MWA and CRAM polar plot. 
When the linear array of isotropic antennas are placed on the xy-plane.
 - on xy-plane
 - theta = [0, 90] degrees
 - phi = [0, 360] degrees
 - d = lamda/2
 - n = 4, 8
 - considering ground plane effect, where h = 0.3/lamda
''' 


delta = 0 # phase difference of one source with respect to another. 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
d = lamda/2

l = np.sin(theta)*np.cos(phi)
m = np.sin(phi)*np.sin(theta)
chix = beta*d*l + delta # total phase difference.
chiy = beta*d*m  + delta # total phase difference.

# Considering the element pattern: 
h = 0.3/lamda
ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))


# voltage response for MWA. 
n = 4
afx = 0; afy = 0
for k in range(n):
    afx = afx + np.exp((1j*chix*k))
    afy = afy + np.exp((1j*chiy*k))

mwa_af = (afx*afy)
mwa_res = abs(mwa_af*ground_plane_effect)


# voltage response for CRAM. 
n = 8
afx = 0; afy = 0
for k in range(n):
    afx = afx + np.exp((1j*chix*k))
    afy = afy + np.exp((1j*chiy*k))

cram_af = (afx*afy)
cram_res = abs(cram_af*ground_plane_effect)

plt.figure()
cs = plt.pcolormesh(l, m, mwa_res, norm=colors.SymLogNorm(linthresh=0.3, linscale= 0.3 ,vmin=np.min(mwa_res), vmax= np.max(mwa_res), base=10), cmap='jet_r', shading='auto') 
bar = plt.colorbar(cs, orientation='vertical', extend='both')
# plt.title('voltage response of MWA')
plt.figure()
cs = plt.pcolormesh(l, m, cram_res, norm=colors.SymLogNorm(linthresh=0.3, linscale=0.3 ,vmin=np.min(cram_res), vmax= np.max(cram_res), base=10), cmap='jet_r', shading='auto') 
# plt.title('voltage response of CRAM')
plt.title('voltage response of MWA and CRAM')
plt.xlabel('l')
plt.ylabel('m')
plt.show()

#%%
''''
This code is written only to create a subplot for MWA and CRAM polar plot. 
When the linear array of isotropic antennas are placed on the xy-plane.
 - on xy-plane
 - theta = [0, 90] degrees
 - phi = [0, 360] degrees
 - d = lamda/2 centered at the centre of the tile
 - n = 4, 8
 - considering ground plane effect, where h = 0.3/lamda
''' 

delta = 0 # phase difference of one source with respect to another. 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
l = np.sin(theta)*np.cos(phi)
m = np.sin(theta)*np.sin(phi)

# Considering the element pattern: 
h = 0.3/lamda
ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))


# voltage response for MWA. 
n = 4
afx = 0; afy = 0
drange = ((n-1)*(lamda/2))/2.
d = np.linspace(-drange, drange, n)
mwa_af = 0
for p in range(n):
    for q in range(n):
        mwa_af = mwa_af + np.exp(1j*beta*(l*d[q]+m*d[p]))

# mwa_res = abs(mwa_af*ground_plane_effect)
mwa_res = abs(mwa_af*ground_plane_effect)

# voltage response for CRAM. 
n = 8
afx = 0; afy = 0
drange = ((n-1)*(lamda/2))/2.
d = np.linspace(-drange, drange, n)
cram_af = 0
for p in range(n):
    for q in range(n):
        cram_af = cram_af + np.exp(1j*beta*(l*d[q]+m*d[p]))
        
# cram_res = abs(cram_af*ground_plane_effect)
cram_res = abs(cram_af*ground_plane_effect)

fontsize = 20
markersize=10
plt.figure(figsize=(13, 12))
# mwa_res[:, 23:] = 0
cs = plt.pcolormesh(l, m, mwa_res, norm=colors.SymLogNorm(linthresh=0.3, linscale= 0.3 ,vmin=np.nanmin(mwa_res), vmax= np.nanmax(mwa_res), base=10), cmap='jet_r', shading='auto') 
# cs.set_clim(0, 10)
bar = plt.colorbar(cs, ticks=[0, 1, 10, 100, 10**2], orientation='vertical', extend='both')
bar.ax.tick_params(labelsize=fontsize)  
# plt.title('voltage response of MWA')
plt.xlabel('directional cosine l',  fontsize=fontsize)
plt.ylabel('directional cosine m', fontsize=fontsize)
plt.xticks(fontsize=fontsize)  # Increase the font size of the x-axis tick labels to 10
plt.yticks(fontsize=fontsize)  # Increase the font size of the y-axis tick labels to 10
os.chdir('D:\\University of Curtin\\EoR foreground mitigation with the CRAM\\Paper1\\Figures')
plt.savefig('primary_beam_MWA.pdf', format='pdf', dpi=600)  # Adjust DPI as needed

fontsize = 20
markersize=10
plt.figure(figsize=(13, 12))
# cram_res[:, 13:] = 0
cs = plt.pcolormesh(l, m, cram_res, norm=colors.SymLogNorm(linthresh=0.3, linscale=0.3 ,vmin=np.nanmin(cram_res), vmax= np.nanmax(cram_res), base=10), cmap='jet_r', shading='auto') 
bar = plt.colorbar(cs, ticks=[0, 1, 10, 100, 10**2], orientation='vertical', extend='both')
bar.ax.tick_params(labelsize=fontsize)  
# plt.title('voltage response of CRAM')
# cs.set_clim(0, 10)
plt.xlabel('directional cosine l',  fontsize=fontsize)
plt.ylabel('directional cosine m',  fontsize=fontsize)
plt.xticks(fontsize=fontsize)  # Increase the font size of the x-axis tick labels to 10
plt.yticks(fontsize=fontsize)  # Increase the font size of the y-axis tick labels to 10
os.chdir('D:\\University of Curtin\\EoR foreground mitigation with the CRAM\\Paper1\\Figures')
plt.savefig('primary_beam_CRAM.pdf', format='pdf', dpi=600)  # Adjust DPI as needed


print (np.sum(cram_res)/np.sum(mwa_res))


# def source_position(pos1, pos2, N_sources):
#     theta = np.linspace(pos1, pos2, N_sources)
#     for i in range(N_sources):     
#         S[i] = [theta[i], 90, 1]
#     return S

# pos1 = 0
# pos2 = 25
# N_sources = 10
# S = np.zeros((N_sources, 3))
# S = source_position(pos1, pos2, N_sources)
# l1 = np.sin(np.deg2rad(S[:, 0]))*np.cos(np.deg2rad(S[:, 1]))
# m1 = np.sin(np.deg2rad(S[:, 0]))*np.sin(np.deg2rad(S[:, 1]))
 
# for i in range(len(l1)): 
#     plt.plot(l1[i], m1[i], '-or')
    

# plt.savefig(save_results_to + 'power_res_MWA_CRAM.png',  bbox_inches="tight", dpi = 300)

#%% array factor for different frequencies: 
    
delta = 0 # phase difference of one source with respect to another. 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
l = np.sin(theta)*np.cos(phi)
m = np.sin(theta)*np.sin(phi)
f = np.arange(80, 330, 10)*10**6


for itr in range(len(f)):
    C = (3*10**8)
    lamda = (C/f[itr])
    beta = (2*np.pi/ lamda)
    # print(" The frequency: ", f[itr], "Wavelength: ", lamda)

    

    # voltage response for MWA. 
    n = 8
    drange = ((n-1)*1.1)/2.
    d = np.linspace(-drange, drange, n)
    mwa_af = 0
    
    for p in range(n):
        for q in range(n):
            mwa_af = mwa_af + np.exp(1j*beta*(l*d[p]+m*d[q]))
            # print(d[p], d[q])
    
       
    h = 0.3/lamda
    ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))

    
    gain_ns = np.sqrt(1 - l*l)
    gain_ew = np.sqrt(1 - m*m)

    res_yy = gain_ns*mwa_af*ground_plane_effect
    res_xx = gain_ew*mwa_af*ground_plane_effect

    p_yy = np.real(res_yy*np.conj(res_yy))
    p_xx =  np.real(res_xx*np.conj(res_xx))
    
    

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(l, m, gain_ns, rstride=1, cstride=1, antialiased=True)
    plt.title("NS Gain at: " +str(int(f[itr]/(10**6))) + "MHz")
    plt.xlabel('l')
    plt.ylabel('m')
    plt.show()

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(l, m, gain_ew, rstride=1, cstride=1, antialiased=True)
    plt.title("EW Gain at: " +str(int(f[itr]/(10**6))) + "MHz")
    plt.xlabel('l')
    plt.ylabel('m')
    plt.show()

    # print(np.min(gain_ns), np.max(gain_ns))
    # print(np.min(gain_ew), np.max(gain_ew))
    # print("--------------------------------------------------------")
    # plt.figure()
    # cs = plt.pcolormesh(l, m, p_yy, norm=colors.SymLogNorm(linthresh=0.3, linscale= 0.3 ,vmin=np.nanmin(p_yy), vmax= np.nanmax(p_yy), base=10), cmap='jet_r', shading='auto') 
    # bar = plt.colorbar(cs, ticks=[0, 1, 10, 100, 10**2], orientation='vertical', extend='both', label=r'$\log T\ [\mathrm{K}]$')
    # plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
    # plt.ylabel('directional cosine m', fontsize=12, fontweight ='bold')
    # plt.title("Power response at: " +str(int(f[itr]/(10**6))) + "MHz")

    # plt.figure()
    # cs = plt.pcolormesh(l, m, p_xx, norm=colors.SymLogNorm(linthresh=0.3, linscale= 0.3 ,vmin=np.nanmin(p_xx), vmax= np.nanmax(p_xx), base=10), cmap='jet_r', shading='auto') 
    # bar = plt.colorbar(cs, ticks=[0, 1, 10, 100, 10**2], orientation='vertical', extend='both', label=r'$\log T\ [\mathrm{K}]$')
    # plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
    # plt.ylabel('directional cosine m', fontsize=12, fontweight ='bold')
    # plt.title("Power response at: " +str(int(f[itr]/(10**6))) + "MHz")
#%%
f = [100*10**6, 150*10**6, 200*10**6, 300*10**6]
delta = 0 # phase difference of one source with respect to another. 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
d = 0.6

n = 8
l = np.sin(theta)*np.cos(phi)
m = np.sin(phi)*np.sin(theta)
for k in range(len(f)):
    lamda = ((3*10**8)/f[k])
    beta = (2*np.pi/ lamda)
    
  
    chix = beta*d*l + delta # total phase difference.
    chiy = beta*d*m + delta # total phase difference.
    
    afx = 0; afy = 0
    
    for i in range(n):
        afx = afx + np.exp((1j*chix*i))
        afy = afy + np.exp((1j*chiy*i))
    
    af = abs(afx*afy)
    
    # Considering the element pattern: 
    h = 0.3/lamda
    ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))
    
    total_voltage = af*ground_plane_effect
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(l, m, total_voltage, rstride=1, cstride=1, antialiased=True)
    plt.title(label = str(f[k]/10**6) + " MHz" +  " and n=8")
    plt.show()

    plt.figure()
    plt.pcolormesh(l, m, total_voltage) 
    plt.title(label = str(f[k]/10**6) + " MHz" +  " and n=8")
    plt.show()


#%% MWAxCRAM and MWAxMWA in 1D



f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
d = lamda/2

delta = 0 # phase difference of one source with respect to another. 
angle = np.arange(0, 180+1) # its covers the horizon. 
theta = np.deg2rad(angle) 
l = np.cos(theta)
chi = beta*d*l 

h = 0.3/lamda 
angle = np.arange(-90, 90+1)
theta = np.deg2rad(angle)
element_pattern = np.sin(2*h*np.pi*np.cos(theta))#ground plane effect


n = 4
mwa_af = ((np.sin(n*chi/2))/(n*np.sin(chi/2)))*element_pattern #normalised array factor. 

n = 8
cram_af = ((np.sin(n*chi/2))/(n*np.sin(chi/2)))*element_pattern #normalised array factor. 


plt.figure()
#MWA x CRAM
mwa_cram = mwa_af*cram_af
plt.plot(l, mwa_cram)

#MWA x MWA
mwa_mwa = mwa_af*mwa_af
plt.plot(l, mwa_mwa)

#CRAM-CRAM
cram_cram = cram_af*cram_af
plt.plot(l, cram_cram)


# plt.title("Power Response for " + str(f/10**6) + " MHz" )
plt.xlabel('directional cosine l', 
               fontweight ='bold')
plt.ylabel('normalised amplitude', 
               fontweight ='bold')
plt.legend(["MWA-CRAM", "MWA-MWA", 'CRAM-CRAM'], loc = 'best')



#%% MWA x CRAM in 2D

delta = 0 # phase difference of one source with respect to another. 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
# d = lamda/2




l = np.sin(theta)*np.cos(phi)
m = np.sin(phi)*np.sin(theta)
chix = beta*d*l + delta # total phase difference.
chiy = beta*d*m  + delta # total phase difference.

# Considering the element pattern: 
h = 0.3/lamda
ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))


# voltage response for MWA. 
n = 4
drange = ((n-1)*(1.1))/2.
d = np.linspace(-drange, drange, n)
afx = 0; afy = 0
for k in range(n):
    afx = afx + np.exp((1j*chix*k))
    afy = afy + np.exp((1j*chiy*k))

mwa_af = (afx*afy)
mwa_res = mwa_af*ground_plane_effect

# voltage response for CRAM. 
n = 8
afx = 0; afy = 0
for k in range(n):
    afx = afx + np.exp((1j*chix*k))
    afy = afy + np.exp((1j*chiy*k))

cram_af = (afx*afy)
cram_res = cram_af*ground_plane_effect


# MWA x MWA
power_1 = np.real(mwa_res*np.conj(mwa_res))

# MWA x CRAM
power_2 = abs(mwa_res*np.conj(cram_res))
# power_2 = power_2/np.max(power_2)


fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Power Response ')
ax1.imshow(power_1, origin = 'lower',  vmin = np.min(power_1), vmax = np.max(power_1))
ax1.set_title('Power response MWA-MWA')
ax2.imshow(power_2, origin = 'lower',  vmin = np.min(power_2), vmax = np.max(power_2))
ax2.set_title('Power response MWA-CRAM')



fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw=dict(projection='3d'))
fig.suptitle('3D plot of Power response ')
# ax1 = fig.gca(projection='3d')
surf = ax1.plot_surface(l, m, power_1, rstride=1, cstride=1, antialiased=True)
ax1.set_title('Power response MWA-MWA')
# ax1 = fig.gca(projection='3d')
surf = ax2.plot_surface(l, m, power_2, rstride=1, cstride=1, antialiased=True)
ax2.set_title('Power response MWA-CRAM')
plt.show()



# fig, (ax1, ax2) = plt.subplots(1, 2)
# fig.suptitle('Power Response in polar plot')
# ax1.pcolormesh(l, m, power_1) 
# ax1.set_title('Power response MWA-MWA')
# ax2.pcolormesh(l, m, power_2) 
# ax2.set_title('Power response MWA-CRAM')
# plt.show()



fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Power Response')
cs = ax1.pcolormesh(l, m, power_1, norm=colors.SymLogNorm(linthresh=0.3, linscale= 0.3 ,vmin=np.min(power_1), vmax= np.max(power_1), base=10), cmap='jet_r', shading='auto')
cs = ax2.pcolormesh(l, m, power_2, norm=colors.SymLogNorm(linthresh=0.3, linscale= 0.3 ,vmin=np.min(power_2), vmax= np.max(power_2), base=10), cmap='jet_r', shading='auto')
bar = plt.colorbar(cs, orientation='vertical', extend='both')
plt.show()


# cs = plt.pcolormesh(l, m, power_1, norm=colors.SymLogNorm(linthresh=0.3, linscale= 0.3 ,vmin=np.min(power_1), vmax= np.max(power_1), base=10), cmap='jet_r', shading='auto')
# bar = plt.colorbar(cs, orientation='vertical', extend='both')
# plt.title('Power Response of MWA x MWA')

#%% Power response of MWA x CRAM 

delta = 0 # phase difference of one source with respect to another. 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
# d = lamda/2
l = np.sin(theta)*np.cos(phi)
m = np.sin(theta)*np.sin(phi)
# chix = beta*d*l + delta # total phase difference.
# chiy = beta*d*m + delta # total phase difference.

# Considering the element pattern: 
h = 0.3/lamda
ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))


# voltage response for MWA. 
n = 4
afx = 0; afy = 0
drange = ((n-1)*(lamda/2))/2.
d = np.linspace(-drange, drange, n)


for idx in range(n):
    chix = beta*d*l + delta 
    afx = afx + np.exp(1j*(idx)*chix)


for idx in range(n):
    chiy = beta*d*m + delta 
    afy = afy+ np.exp(1j*(idx)*chiy)




        
mwa_af = (afx*afy)
mwa_res = mwa_af*ground_plane_effect


# voltage response for CRAM. 
n = 8
afx = 0; afy = 0
drange = ((n-1)*(lamda/2))/2.
d = np.linspace(-drange, drange, n)
for idx in range(n):
    chix = beta*d*l + delta 
    afx = afx + np.exp(1j*(idx)*chix)

for idx in range(n):
    chiy = beta*d*m + delta 
    afy = afy+ np.exp(1j*(idx)*chiy)

cram_af = (afx*afy)
cram_res = cram_af*ground_plane_effect


# MWA x MWA
power_1 = np.real(mwa_res*np.conj(mwa_res))

# MWA x CRAM
power_2 = np.real(mwa_res*np.conj(cram_res))


plt.figure()
cs = plt.pcolormesh(l, m, power_1, norm=colors.SymLogNorm(linthresh= 0.3, linscale= 0.3 ,vmin=np.nanmin(power_1), vmax= np.nanmax(power_1), base=10), cmap='jet_r', shading='auto')
bar = plt.colorbar(cs, orientation='vertical', extend='both')
plt.title('Power Response of MWA x CRAM')
plt.xlabel('l')
plt.ylabel('m')


plt.figure()
cs = plt.pcolormesh(l, m, power_2, norm=colors.SymLogNorm(linthresh= 0.3, linscale= 0.3 ,vmin=np.nanmin(power_2), vmax= np.nanmax(power_2), base=10), cmap='jet_r', shading='auto')
bar = plt.colorbar(cs, orientation='vertical', extend='both')
plt.title('Power Response of MWA x CRAM')
plt.xlabel('l')
plt.ylabel('m')


#%% Superimposition of power response of MWA x CRAM with simulated sources. 
#This.
delta = 0 # phase difference of one source with respect to another. 
theta = np.deg2rad(np.linspace(0, 90+1, 513))
phi = np.deg2rad(np.linspace(0, 360+1, 513))
[theta, phi] = np.meshgrid(theta, phi)
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
l = np.sin(theta)*np.cos(phi)
m = np.sin(theta)*np.sin(phi)

# Considering the element pattern: 
h = 0.3/lamda
ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))


# voltage response for MWA. 
n = 4
afx = 0; afy = 0
drange = ((n-1)*(1.1))/2.
d = np.linspace(-drange, drange, n)
# d = np.linspace(0, (n-1)*lamda/2, n)
mwa_af = 0
for p in range(n):
    for q in range(n):
        mwa_af = mwa_af + np.exp(1j*beta*(l*d[q]+m*d[p]))

mwa_res = mwa_af*ground_plane_effect


# voltage response for CRAM. 
n = 8
afx = 0; afy = 0
drange = ((n-1)*(1.1))/2.
d = np.linspace(-drange, drange, n)
# d = np.linspace(0, (n-1)*lamda/2, n)
cram_af = 0
for p in range(n):
    for q in range(n):
        cram_af = cram_af + np.exp(1j*beta*(l*d[q]+m*d[p]))



cram_res = cram_af*ground_plane_effect


# MWA x MWA
power_1 = np.real(mwa_res*np.conj(mwa_res))
power_1 = power_1/np.max(power_1)

# MWA x CRAM
power_2 = np.real(mwa_res*np.conj(cram_res))
power_2 = power_2/np.max(power_2)

angle = np.linspace(0, 181, 513)
plt.plot(angle, power_1[256, :], '-k')
plt.plot(angle, power_2[256, :], '-r')

# plt.figure()
# cs = plt.pcolormesh(l, m, power_1, norm=colors.SymLogNorm(linthresh= 0.3, linscale= 0.3 ,vmin=np.nanmin(power_1), vmax= np.nanmax(power_1), base=10), cmap='jet_r', shading='auto')
# bar = plt.colorbar(cs,ticks=[0, 1, 10, 100], orientation='vertical', extend='both')
# # plt.title('Power Response of MWA x MWA')
# plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
# plt.ylabel('directional cosine m',  fontsize=12, fontweight ='bold')



# plt.figure()
# cs = plt.pcolormesh(l, m, power_2, norm=colors.SymLogNorm(linthresh= 0.3, linscale= 0.3 ,vmin=np.nanmin(power_2), vmax=np.nanmax(power_2), base=10), cmap='jet_r', shading='auto')
# bar = plt.colorbar(cs,ticks=[-100,-10, 0, 10, 100], orientation='vertical', extend='both')
# # plt.title('Power Response of MWA x CRAM')
# #bar.set_label('X+Y')
# plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
# plt.ylabel('directional cosine m',  fontsize=12, fontweight ='bold')


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(l, m, power_1, rstride=1, cstride=1, antialiased=True)
# plt.title('Power Response of MWA x MWA')
# plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
# plt.ylabel('directional cosine m',  fontsize=12, fontweight ='bold')
# plt.show()

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(l, m, power_2, rstride=1, cstride=1, antialiased=True)
# plt.title('Power Response of MWA x CRAM')
# plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
# plt.ylabel('directional cosine m',  fontsize=12, fontweight ='bold')
# plt.show()


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(l, m, abs(mwa_res), rstride=1, cstride=1, antialiased=True)
# plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
# plt.ylabel('directional cosine m',  fontsize=12, fontweight ='bold')
# plt.show()


# # 3D plot of voltage response without the ground plane effect. 
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(l, m, abs(mwa_af), rstride=1, cstride=1, antialiased=True)
# plt.xlabel('directional cosine l',  fontsize=12, fontweight ='bold')
# plt.ylabel('directional cosine m',  fontsize=12, fontweight ='bold')
# plt.show()

#%%
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors
from matplotlib.ticker import LogFormatter 

cs = plt.pcolormesh(l, m, power_2, norm=colors.SymLogNorm(linthresh= 0.3, linscale= 0.3 ,vmin=-100, vmax=100, base=10), cmap='jet_r', shading='auto')

cb = plt.colorbar(cs, ticks=[-100,-10, 0, 10, 100], orientation='vertical', extend='both', format = LogFormatter)
#%%


# Sources kept between wedges. 
def source_position(pos1, pos2, N_sources):
    theta = np.linspace(pos1, pos2, N_sources)
    for i in range(N_sources):     
        S[i] = [theta[i], 60, 1]
    return S

pos1 = 0
pos2 = 10
N_sources = 10
S = np.zeros((N_sources, 3))
S = source_position(pos1, pos2, N_sources)
l1 = np.sin(np.deg2rad(S[:, 0]))*np.cos(np.deg2rad(S[:, 1]))
m1 = np.sin(np.deg2rad(S[:, 0]))*np.sin(np.deg2rad(S[:, 1]))

 
for i in range(len(l1)): 
    plt.plot(l1[i], m1[i], '-ok')
# plt.savefig(save_results_to + 'power_res_MWA_CRAM.png',  bbox_inches="tight", dpi = 300)


    
plt.figure()
th_polar = np.deg2rad(S[:, 0])
r = np.sin(np.deg2rad(S[:, 1]))
plt.subplot(111, projection='polar')
plt.plot(th_polar, r, '.k')
plt.title('polar plot of sources placed from zenith to horizon')
# plt.savefig(save_results_to + 'polar_plot_source.png',  bbox_inches="tight", dpi = 300)

    
#%% Voltage Response and Visibility function. 




def voltage_response(f, n, th, ph):
    delta = 0 # phase difference of one source with respect to another. 

    [theta, phi] = np.meshgrid(th, ph)  
    
    lamda = ((3*10**8)/f)
    beta = (2*np.pi/ lamda)
    
    l = np.sin(theta)*np.cos(phi)
    m = np.sin(theta)*np.sin(phi)

    # Considering the element pattern: 
    h = 0.3/lamda
    ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))


    drange = ((n-1)*(1.1))/2.
    d = np.linspace(-drange, drange, n)
    af = 0
    for p in range(n):
        for q in range(n):
            af = af + np.exp(1j*beta*(l*d[q]+m*d[p]))
            
    res = af*ground_plane_effect
        
    return res, l, m
  


def mwa_distance_calculation(filename, f):
    
    wb = xlrd.open_workbook(filename)
    sheet = wb.sheet_by_name("Hex and LB tiles")
    x = []; y = []
    for i in range(93, sheet.nrows):
        x.append(sheet.cell_value(i, 4))
        y.append(sheet.cell_value(i, 5))

    # plt.plot(x, y, '.k')
    
    n = len(x)
    nb = int(n*(n-1)/2.0)
    distance_map = np.zeros((int(nb), 3))
    lamda = ((3*10**8)/f)
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            distance_map[k, 0] = x[i] - x[j] 
            distance_map[k, 1] = y[i] - y[j] 
            k = k+1
    distance_map = distance_map/lamda
        
    # plt.figure(); plt.plot(distance_map[:, 0], distance_map[:, 1], '.k')

    return distance_map


def cram_distance_calculation(filename, f):
    
    wb = xlrd.open_workbook(filename)
    sheet = wb.sheet_by_name("Hex and LB tiles")
    x = []; y = []
    for i in range(93, sheet.nrows):
        x.append(sheet.cell_value(i, 4))
        y.append(sheet.cell_value(i, 5))

    # plt.figure();  plt.plot(x, y, '.k')
    
    nb = len(x)
    distance_map = np.zeros((int(nb), 3))
    lamda = ((3*10**8)/f)
    for i in range(nb):
        distance_map[i, 0] = x[i] - 11.019
        distance_map[i, 1] = y[i] - 122.058 
    
    distance_map = distance_map/lamda
        
    # plt.figure(); plt.plot(distance_map[:, 0], distance_map[:, 1], '.k')

    return distance_map


#%%
filename = ("D:/University of Curtin/EoR foreground mitigation with the CRAM/Data/Merged MWA tile coordinates_Phase2 - AW 2018-06-21.xlsx")


wb = xlrd.open_workbook(filename)
sheet = wb.sheet_by_name("Phase 2 Compact")
x = []; y = []
for i in range(93, sheet.nrows):
    x.append(sheet.cell_value(i, 4))
    y.append(sheet.cell_value(i, 5))

plt.figure(); 
plt.xlim([-40,60])
plt.ylim([80,170]) 
plt.xlabel("east (m) ", fontweight ='bold',fontsize=12)
plt.ylabel("north (m) ", fontweight ='bold',fontsize=12)
for i in range(36):
    plt.plot(x[i], y[i], '.k')
    plt.text(x[i], y[i], i, bbox=dict(facecolor='red', alpha=0.5), fontweight ='bold',fontsize=12 )

plt.plot(11.019,122.058, '*k')
plt.text(11.019, 122.058, "C", bbox=dict(facecolor='blue', alpha=0.5), fontweight ='bold',fontsize=12)

    
cram_distance_map = cram_distance_calculation(filename, 150*10**6)
d1 = np.sqrt((cram_distance_map[:, 0]**2) + (cram_distance_map[:, 1]**2))    
    

mwa_distance_map = mwa_distance_calculation(filename, 150*10**6)
d2 = np.sqrt((mwa_distance_map[:, 0]**2) + (mwa_distance_map[:, 1]**2))

plt.plot(d2, 'k')
plt.plot(d1, '.-r')
plt.xlabel('number of baselines')
plt.ylabel("baseline length (m)")
#%%

def visibility(N, S, F, distance_map, T):    
    nb = len(distance_map)
    l = np.sin(np.deg2rad(S[:, 0]))*np.cos(np.deg2rad(S[:, 1]))
    m = np.sin(np.deg2rad(S[:, 0]))*np.sin(np.deg2rad(S[:, 1]))
    V = np.zeros((nb, 1), dtype=complex)
    
    uv = np.zeros((nb,  2))
    for i in range(nb):
        [u, v, w ] = np.matmul(T, distance_map[i])
        uv[i, :] = [u, v]
        for k in range(N):
            V[i] = V[i] + F(np.deg2rad(S[k, 0]), np.deg2rad(S[k, 1]))[0]*S[k, 2]*np.exp(-2j*np.pi*(u*l[k] + v*m[k]))
            # print(F(np.deg2rad(S[k, 0]), np.deg2rad(S[k, 1]))[0])
            # V[i] = V[i] + P[int(S[k, 1]), int(S[k, 0])]*S[k, 2]*np.exp(-2j*np.pi*(u*l[k] + v*m[k]))
       
    # plt.figure(); plt.plot(uv[:,0], uv[:,1], '.k')
    return V, uv

# Sources kept between wedges. 
def source_position(pos1, pos2):
    theta = np.linspace(pos1, pos2, N_sources)
    for i in range(N_sources):     
        S[i] = [theta[i], 60, 1]
    return S


#%% Visibility comparison between MWA x MWA and MWA x CRAM. 

# General setup: 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
freq = (150*10**6) 

filename = ("D:/University of Curtin/EoR foreground mitigation with the CRAM/Data/Merged MWA tile coordinates_Phase2 - AW 2018-06-21.xlsx")


N_sources = 1
S = np.zeros((N_sources, 3)) #array defined to mention the details of the source under consideration. 
# the columns could be considered to define the values of (theta, phi, intensity)
S[0] = [90, 60, 1] # zenith 
string_name = 'zenith_src_'

# S[0] = [85, 60, 1] # horizon
# string_name = 'horizon_src_'

# S[0] = [26, 36, 1] # -ve cram response
# string_name = 'neg_cram_res_'

# two sources with -ve cram response. 
# S[0] = [36, 26, 1] 
# S[1] = [38, 27, 1] 
# string_name = '2neg_cram_res_'



# two sources, with one pos and one neg cram response. 
# S[0] = [10, 104, 1] 
# S[1] = [15, 104, 1] 
# string_name = 'pos_neg_cram_res_'

# sources spread across two given zenith angles.
# pos1 = 70
# pos2 = 90
# S = source_position(pos1, pos2)
# string_name = 'zenith_horizon_'
# string_name = '20deg_horizon_'
# single snapshot 
HA = np.deg2rad(0) 
DEC = np.deg2rad(90)
trans_matrix = np.array([[np.sin(HA), np.cos(HA), 0], [-np.sin(DEC)*np.cos(HA), np.sin(DEC)*np.sin(HA), np.cos(DEC)], [np.cos(DEC)*np.cos(HA), -np.cos(DEC)*np.sin(HA), np.sin(DEC)]])


#------------------------------MWA x MWA--------------------------------------#
n = 4
mwa_res, l, m = voltage_response(freq, n, theta, phi)
P_mwa =  np.real(mwa_res*np.conj(mwa_res))
P_mwa = P_mwa/np.max(P_mwa)

#interpolation
interp_function_mwa = interpolate.interp2d(theta, phi, P_mwa, kind = 'cubic')

#distance calculation
mwa_distance_map = mwa_distance_calculation(filename, freq)


#visibility calculation
V_mwa, uv_mwa = visibility(N_sources, S, interp_function_mwa, mwa_distance_map, trans_matrix)

# plotting the visibility for MWA x MWA. 
# plt.figure()
# plt.scatter(V_mwa.real, V_mwa.imag, marker =".", c ="black")
# V_conj = np.conj(V_mwa)
# plt.scatter(V_conj.real, V_conj.imag, marker =".", c ="red")
# plt.title("Visibility function: MWA x MWA, f = "+ str(freq/10**6) + " MHz" )
# plt.grid(True)
# plt.show()

# ------------------------- MWA x CRAM ---------------------------------------#
n = 4
mwa_res, l, m = voltage_response(freq, n, theta, phi)

n = 8
cram_res, l, m = voltage_response(freq, n, theta, phi)




P_cram = np.real(mwa_res*np.conj(cram_res))
P_cram = P_cram/np.max(P_cram)

#interpolation
interp_function_cram = interpolate.interp2d(theta, phi, P_cram, kind = 'cubic')

#distance calculation
cram_distance_map = cram_distance_calculation(filename, freq)

#visibility calculation
V_cram, uv_cram = visibility(N_sources, S, interp_function_cram, cram_distance_map, trans_matrix)

# #plotting the visibility for MWA x CRAM. 
# plt.figure()
# plt.scatter(V_cram.real, V_cram.imag, marker =".", c ="black")
# V_conj = np.conj(V_cram)
# plt.scatter(V_conj.real, V_conj.imag, marker =".", c ="red")
# plt.title("Visibility function: MWA x CRAM, f = "+ str(freq/10**6) + " MHz" )
# plt.grid(True)
# plt.show()

#------------------------Plotting for visibility function --------------------#

fig = plt.figure()
plt.scatter(V_mwa.real, V_mwa.imag, marker =".", c ="black" , label = 'MWA x MWA')
plt.scatter(np.conj(V_mwa).real, np.conj(V_mwa).imag, marker =".", c ="black" )
plt.scatter(V_cram.real, V_cram.imag, marker ="o", c ="red", label = 'MWA x CRAM' )
plt.scatter(np.conj(V_cram).real, np.conj(V_cram).imag, marker ="o", c ="red" )
plt.legend( bbox_to_anchor=(1.11,1.), loc = 'best')
plt.xlabel('real components')
plt.ylabel('imaginary components')
plt.grid(True)

# plt.title("Visibility function for {N} source at (\u03B8, \u03C6) = ({val1}, {val2}) degrees".format(N = N_sources, val1=S[0,0], val2=S[0,1]))
# plt.title("Visibility function for {N} source at \n (\u03B8, \u03C6) = ({val1}, {val2}) degrees, \n (\u03B8, \u03C6) = ({val3}, {val4}) degrees".format(N = N_sources, val1=S[0,0], val2=S[0,1], val3 = S[1,0], val4 = S[1,1] ))
# plt.title("Visibility function for {N} sources kept between angles of ({val1} & {val2}) degrees".format(N = N_sources, val1=pos1, val2=pos2))
# plt.title("Visibility function for {N} sources spread around horizon, kept between angles of ({val1} & {val2}) degrees".format(N = N_sources, val1=pos1, val2=pos2))


# fig.savefig(save_results_to + string_name +'visibility.png',  bbox_inches="tight", dpi = 300)

# fig.savefig(save_results_to + 'wedge'+'_'+str(pos1)+'_'+str(pos2)+'_'+'visibility.png',  bbox_inches="tight", dpi = 300)
# plt.show()
# plt.close()

fig = plt.figure()
plt.scatter(np.sqrt(uv_mwa[:, 0]**2 + uv_mwa[:, 1]**2), abs(V_mwa), marker =".", c ="black", label = 'MWA x MWA' )
plt.scatter(np.sqrt(uv_mwa[:, 0]**2 + uv_mwa[:, 1]**2), abs(np.conj(V_mwa)), marker =".", c ="black" )
plt.scatter(np.sqrt(uv_cram[:, 0]**2 + uv_cram[:, 1]**2), abs(V_cram), marker =".", c ="red", label = 'MWA x CRAM' )
plt.scatter(np.sqrt(uv_cram[:, 0]**2 + uv_cram[:, 1]**2), abs(np.conj(V_cram)), marker =".", c ="red" )
plt.xlabel(str("$\sqrt{u^2 + v^2}$"))
plt.ylabel(str("$|V|$"))
plt.legend( loc = 'best')


# plt.title("Amplitude of visibility function Vs uv-plane: {N} source at (\u03B8, \u03C6) = ({val1}, {val2}) degrees".format(N = N_sources, val1=S[0,0], val2=S[0,1]))
# plt.title("Amplitude of visibility function Vs uv-plane: {N} source at \n (\u03B8, \u03C6) = ({val1}, {val2}) degrees, \n (\u03B8, \u03C6) = ({val3}, {val4}) degrees".format(N = N_sources, val1=S[0,0], val2=S[0,1], val3 = S[1,0], val4 = S[1,1] ))
# plt.title("Amplitude of visibility function Vs uv-plane for {N} sources kept between angles of ({val1} & {val2}) degrees".format(N = N_sources, val1=pos1, val2=pos2))
# plt.title("Amplitude of visibility function Vs uv-plane for {N} sources spread around horizon, kept between angles of ({val1} & {val2}) degrees".format(N = N_sources, val1=pos1, val2=pos2))



# fig.savefig(save_results_to + string_name+'visibility_amp.png', bbox_inches="tight", dpi = 300)
# fig.savefig(save_results_to + 'wedge'+'_'+str(pos1)+'_'+str(pos2)+'_'+'visibility_amp.png', bbox_inches="tight", dpi = 300)

# plt.show()

# plt.close()


#%% Amplitude of Visibility wrt the change of theta from horizon to zenith. 

'''
Wedges considered between 0 and 90 degrees and then plotting the visibility 
function wrt the change of theta. 

'''
# General setup: 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
freq = (150*10**6) 

filename = ("D:/University of Curtin/EoR foreground mitigation with the CRAM/Data/Merged MWA tile coordinates_Phase2 - AW 2018-06-21.xlsx")


N_sources = 1
S = np.zeros((N_sources, 3)) #array defined to mention the details of the source under consideration. 
delta_theta = 2
N_itr =len(np.arange(0, 90, delta_theta))

Vmwa = np.zeros((630, N_itr))
Vcram = np.zeros((36, N_itr))

# single snapshot 
HA = np.deg2rad(0)
DEC = np.deg2rad(90)
trans_matrix = np.array([[np.sin(HA), np.cos(HA), 0], [-np.sin(DEC)*np.cos(HA), np.sin(DEC)*np.sin(HA), np.cos(DEC)], [np.cos(DEC)*np.cos(HA), -np.cos(DEC)*np.sin(HA), np.sin(DEC)]])

for itr in range(N_itr):
    # sources spread across two given zenith angles.
    pos1 = itr*delta_theta
    pos2 = pos1 + delta_theta
    S = source_position(pos1, pos2)
    

    
    #------------------------------MWA x MWA--------------------------------------#
    n = 4
    mwa_res, l, m = voltage_response(freq, n, theta, phi)
    P_mwa =  np.real(mwa_res*np.conj(mwa_res))
    P_mwa = P_mwa/np.max(P_mwa)
    
    #interpolation
    interp_function_mwa = interpolate.interp2d(theta, phi, P_mwa, kind = 'cubic')
    
    #distance calculation
    mwa_distance_map = mwa_distance_calculation(filename, freq)
    
    
    #visibility calculation
    V_mwa, uv_mwa = visibility(N_sources, S, interp_function_mwa, mwa_distance_map, trans_matrix)
    
    # ------------------------- MWA x CRAM ---------------------------------------#
    n = 4
    mwa_res, l, m = voltage_response(freq, n, theta, phi)
    
    n = 8
    cram_res, l, m = voltage_response(freq, n, theta, phi)
    
    
    P_cram = np.real(mwa_res*np.conj(cram_res))
    P_cram = P_cram/np.max(P_cram)
    
    #interpolation
    interp_function_cram = interpolate.interp2d(theta, phi, P_cram, kind = 'cubic')
    
    #distance calculation
    cram_distance_map = cram_distance_calculation(filename, freq)
    
    #visibility calculation
    V_cram, uv_cram = visibility(N_sources, S, interp_function_cram, cram_distance_map, trans_matrix)

    Vmwa[:, itr] = abs(V_mwa[:, 0])
    Vcram[:, itr] = abs(V_cram[:,0])
    print("Completed iteration: ", itr)
    
#%%
plt.figure()
vmwa = np.mean(Vmwa, axis = 0)
vmwa = np.flip(vmwa)
vcram = np.mean(Vcram, axis = 0)   
vcram = np.flip(vcram)
xaxis = (np.arange(0, 90, delta_theta))
plt.plot(xaxis, vmwa, 'ok')
plt.plot(xaxis, vcram, 'or')
plt.yscale("log")
plt.ylim([10**-7, 10**1])
plt.ylabel("flux density (Jy)", fontweight = 'bold', fontsize=12)
plt.xlabel("elevation angle (deg.)", fontweight ='bold',fontsize=12)
plt.legend(['MWA-MWA', 'MWA-CRAM'])
#%%
plt.figure()
plt.xlim([0, 90])
count = N_itr - 1
y1 = []
y2 = []
for itr in range(N_itr): 
    pos = itr*delta_theta
    x1 = np.linspace(0, pos+delta_theta, 630*(itr+1))
    x2 = np.linspace(0, pos+delta_theta, 36*(itr+1))
    y1 = np.append(y1, Vmwa[:, count])
    y2 = np.append(y2, Vcram[:, count])
    plt.scatter(x1, y1, c = 'black')
    plt.scatter(x2, y2, c = 'red')
    # print(itr, pos, count)
    count = count-1
    
plt.yscale("log")
# plt.xscale("log")
plt.ylim([10**-7, 10**1])
# plt.xlim([10**0, 10**2])

plt.ylabel("visibility measurements", fontweight = 'bold', fontsize=12)
plt.xlabel("elevation angle", fontweight ='bold',fontsize=12)
plt.legend(['MWA-MWA', 'MWA-CRAM'])
# plt.title('Amplitude of visibility vs \u0394\u03B8')

# plt.savefig(save_results_to + 'theta_change_visibility_amp.png',  bbox_inches="tight", dpi = 300)


#%% Integrating LoBES catelogue. 

from PyAstronomy import pyasl # type: ignore
import matplotlib.pylab as plt
import datetime
import numpy as np

# Convert calendar date to JD
# use the datetime package
jd = datetime.datetime(2022, 3, 8)
jd = pyasl.jdcnv(jd)
# Specific RA and DEC
ra = 0.13175539914686585
dec =-2.861792715380515
print()
print("Get horizontal coordinates (alt, az, ha) from JD, RA,")

print(pyasl.eq2hor(jd, ra, dec, lon=116.67083333, lat=-26.70331941, alt=377.827))

print()
print("From a list of Julian dates ...")
jds = np.arange(jd, jd+1, .2)
ras = np.zeros(jds.size) + ra
decs = np.zeros(jds.size) + dec
alt, az, ha = pyasl.eq2hor(jd, ra, dec, lon=116.67083333, lat=-26.70331941, alt=377.827)

for i in range(alt.size):
    print("JD = %g : alt = % g,  az = % g,  ha = % g" % (jds[i], alt[i], az[i], ha[i]))


print()
print("For one object and different times at the VLT...")
jds = np.arange(jd-.25, jd+.25, .01)
ras = np.zeros(jds.size) + 130.
decs = np.zeros(jds.size) - 30.
res = pyasl.eq2hor(jds, ras, decs, lon=-70.4042, lat=-24.6272, alt=2635.)

plt.plot(jds, res[0])
plt.xlabel("Julian date")
plt.ylabel("Altitude [deg]")
plt.show()



#%% Calculating the FWHM using the 1D beam. 

delta = 0 # phase difference of one source with respect to another. 
theta = np.deg2rad(np.arange(0, 90+1))
phi = np.deg2rad(np.arange(0, 360+1))
[theta, phi] = np.meshgrid(theta, phi)
f = (150*10**6) 
lamda = ((3*10**8)/f)
beta = (2*np.pi/ lamda)
# d = lamda/2
l = np.sin(theta)*np.cos(phi)
m = np.sin(theta)*np.sin(phi)

# Considering the element pattern: 
h = 0.3/lamda
ground_plane_effect = np.sin(2*h*np.pi*np.cos(theta))


# voltage response for MWA. 
n = 4
afx = 0; afy = 0
drange = ((n-1)*(1.1))/2.
d = np.linspace(-drange, drange, n)
# d = np.linspace(0, (n-1)*lamda/2, n)
mwa_af = 0
for p in range(n):
    for q in range(n):
        mwa_af = mwa_af + np.exp(1j*beta*(l*d[q]+m*d[p]))

mwa_res = mwa_af*ground_plane_effect
mwa_res = mwa_res/np.max(mwa_res)

# voltage response for CRAM. 
n = 8
afx = 0; afy = 0
drange = ((n-1)*(1.1))/2.
d = np.linspace(-drange, drange, n)
# d = np.linspace(0, (n-1)*lamda/2, n)
cram_af = 0
for p in range(n):
    for q in range(n):
        cram_af = cram_af + np.exp(1j*beta*(l*d[q]+m*d[p]))



cram_res = cram_af*ground_plane_effect
cram_res = cram_res/np.max(cram_res)

# MWA x MWA
power_1 = np.real(mwa_res*np.conj(mwa_res))

# MWA x CRAM
power_2 = np.real(mwa_res*np.conj(cram_res))


plt.figure()
cs = plt.pcolormesh(l, m, power_1, norm=colors.SymLogNorm(
            linthresh= 0.3, 
            linscale= 0.3,
            vmin=np.nanmin(power_1), 
            vmax= np.nanmax(power_1), 
            base=10), 
            cmap='jet_r', 
            shading='auto')
bar = plt.colorbar(cs, orientation='vertical', extend='both')
plt.title('Power Response of MWA x MWA')
plt.xlabel('l')
plt.ylabel('m')


plt.figure()
cs = plt.pcolormesh(l, m, power_2, norm=colors.SymLogNorm(
            linthresh= 0.3, 
            linscale= 0.3, 
            vmin=np.nanmin(power_2), 
            vmax= np.nanmax(power_2), base=10), cmap='jet_r', 
            shading='auto')
bar = plt.colorbar(cs, orientation='vertical', extend='both')
plt.title('Power Response of MWA x CRAM')
plt.xlabel('l')
plt.ylabel('m')
    

#%%
import numpy as np

# Generate some mock data
np.random.seed(0)
x = np.linspace(0, 10, 100)
y = 2 * x + 1 + np.random.normal(0, 1, 100)

# Perform linear regression
X = np.vstack([x, np.ones(len(x))]).T
beta_hat = np.linalg.inv(X.T @ X) @ X.T @ y
residuals = y - X @ beta_hat
variance = np.var(residuals)

# Compute the Fisher matrix
F = (X.T @ X) / variance

# Compute the Cramer-Rao bounds
CRLB = np.linalg.inv(F)

print("Estimated coefficients:", beta_hat)
print("Fisher matrix:")
print(F)
print("Cramer-Rao lower bounds:")
print(CRLB)


#%%
import numpy as np

# Define the model parameters
m = 2.0  # slope
c = 1.0  # intercept

# Generate some mock data
x = np.linspace(0, 10, 100)  # x values
y = m * x + c #+ np.random.normal(0, 1, 100)  # y values with noise

# Compute the Fisher matrix
F = np.zeros((2, 2))  # Fisher matrix
for i in range(len(x)):
    # Compute the derivatives of the observables
    dL_dm = x[i]
    dL_dc = 1.0

    # Update the Fisher matrix
    F[0, 0] += dL_dm ** 2
    F[0, 1] += dL_dm * dL_dc
    F[1, 0] += dL_dc * dL_dm
    F[1, 1] += dL_dc ** 2

# Invert the Fisher matrix
F_inv = np.linalg.inv(F)

print("Fisher matrix:")
print(F)
print("Inverse Fisher matrix:")
print(F_inv)

#%%
import numpy as np

# Generate some example data
x = np.linspace(0, 10, 100)
y = 2 * x + 1
sigma = 0.5  # Gaussian uncertainty of the data points

# Define the model function
def model(x, m, c):
    return m * x + c

# Calculate the Fisher matrix elements
N = len(x)
fisher_11 = np.sum((x**2) / sigma**2) / N
fisher_22 = np.sum(1 / sigma**2) / N
fisher_12 = np.sum(x / sigma**2) / N

# Assemble the Fisher matrix
fisher_matrix = np.array([[fisher_11, fisher_12], [fisher_12, fisher_22]])

# Invert the Fisher matrix to obtain the covariance matrix
covariance_matrix = np.linalg.inv(fisher_matrix)

print("Covariance matrix:")
print(covariance_matrix)
