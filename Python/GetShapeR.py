#! /usr/bin/python3
# -*- coding: utf-8 -*-

""" 
GetShapeR.py provides an estimate of the R-parameter of the cosine scattering law 
R(C+1)(cos(th))^(2C), and computes the radar cross sections and radar albedos over one
rotation in N user-given orientations for a given shape model. Plotting is also an option.

(C) Anne Virkki, Jan 2023.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import linalg

##### INPUT VALUES #####
# Shape model (OBJ) file:
modelname = '98KY26/model01.obj'  
# Size: a0 is the longest axis dimension. Use 0 to not scale the model size.
a0 = 14.0 
# If you want to keep the shape model's axis ratios, use b0, c0 = 0, 0, 
# or scale all axes by providing the intermediate and short axis dimensions
b0, c0 = 0, 0
# Give either the subradar latitude or the set of observer and polar ecliptic coordinates 
# (make delta False for the latter! Note: 0 = False so use e.g. 0.01 for an equatorial view!)
delta = 0					# subradar latitude ([-90,90])
ObsEcLat = 35.0 #-20.3           # Asteroid's observer-centered ecliptic latitude
ObsEcLon = 324.0          		# Asteroid's observer-centered ecliptic longitude
PoleEcLat = -44.3 #-88.0			# Pole ecliptic latitude in degrees
PoleEcLon = 35.6 #45				# Pole ecliptic longitude in degrees
x0 = 10 						# Map dimension (number of pixels per side)
law = 'Cos'						# Scattering law: 'Cos'/'Hag' for Cosine/Hagfors
C = 0.77 						# scattering law parameter (roughness)
rcs_obs = 25					# The expected radar albedo
step_ang = 45					# Angle interval (e.g., 90 for 4 orientations)
P = 321.1                         # Rotation period in seconds
wl = 0.035                      # Observation wavelength
Shadowing = False				# Shadowing (recommended only for significantly concave shapes)
PrintAll = False				# Print all information for every orientation
MakePlots = False				# Plot the cross sections and albedos for every orientation
##### END INPUT #####


def readObj(modelname):
    f = open('%s' % modelname, "r")
    raw = f.readlines()
    f.close()
    vertl = []
    facel = []
    for lines in raw: 
        line = lines.strip().split()
        if line[0] == 'v':
            vertl.append(line[1:])
        elif line[0] == 'f':
            facel.append(line[1:])
        
    verts = np.array(vertl, dtype=float)
    faces = np.array(facel, dtype=int)
	
    return verts, faces;
###
    	
def rotellipsoid(t, axis):
    # Computes the rotation matrix based on the given angle t (in radians)
    cosphi = np.cos(t)
    sinphi = np.sin(t)
    if axis == 'x' or axis == 0:      # roll
        R = np.array([[1.0,  0.0,  0.0],
        			  [0.0,cosphi,-sinphi],
                      [0.0,sinphi, cosphi]])
            
    elif axis == 'y' or axis == 1:    # pitch
        R = np.array([[cosphi,0.0,sinphi],
                      [0.0,  1.0,  0.0],
                      [-sinphi,0.0,cosphi]]) 
             
    elif axis == 'z' or axis == 2:    # yaw
        R = np.array([[cosphi,-sinphi, 0],
                      [sinphi, cosphi, 0],
                      [0.0,  0.0,  1.0]])
    return R;
###	
	
def mesh_vectors(verts,faces):
    # Combine the indices of the triangle corners to the correct coordinates 

    msh = np.zeros((np.shape(faces)[0],3,3))

    for i, face in enumerate(faces):
        for j in range(3):
            msh[i][j] = verts[face[j],:]

    return msh
###
    
def getNormals(verts, faces):
	# Number of orientations and faces
	Nf = np.shape(faces)[0]

	# Find the normals and center points of each facet
	# Combine the indices of the triangle corners to the correct coordinates
	msh = np.zeros((Nf, 3, 3))
	normal = np.zeros((Nf, 3))
	centroid = np.zeros((Nf, 3))
	for i, face in enumerate(faces):
		for k in range(3):
			msh[i][k] = verts[face[k],:]
		
		# Center point of each facet
		centroid[i] = np.mean(msh[i,:,:], axis=0) #np.column_stack((v1, v2, v3)),axis = 1)
		
		# Unit vectors for the normals of the triangles
		cv1 = msh[i,1,:] - msh[i,0,:]
		cv2 = msh[i,2,:] - msh[i,0,:]
		normal[i] = np.cross(cv1, cv2)
		
		# Check if the normal vector is pointing inward and turn it if it is
		if modelname == 'None':
			if np.dot(normal[i], centroid[i]) < 0:
				normal[i] = np.cross(cv2,cv1)
				if savename != 'None':
					# Fix the order of the faces (should be counterclockwise)
					tmp = face[2]
					faces[i,2] = face[1]
					faces[i,1] = tmp

	return normal, centroid, faces
###

def RayTracer(msh, X0, K):
    EPSILON = 1E-10
    edge1 = msh[1,:] - msh[0,:]
    edge2 = msh[2,:] - msh[0,:]
    h = np.cross(K,edge2)
    a = np.dot(edge1,h)
    if (a > -EPSILON) and (a < EPSILON):
        return 0   # This ray is parallel to this triangle.
    f = 1.0/a
    s = X0 - msh[0,:]
    u = f * np.dot(s,h)
    if (u < 0.0) or (u > 1.0):
        return 0
    q = np.cross(s,edge1)
    v = f * np.dot(K,q)
    if (v < 0.0) or (u + v > 1.0):
        return 0
    # At this stage we can compute t to find out where the intersection point is on the line.
    t = f * np.dot(edge2,q)
    if (t > EPSILON):   # ray intersection
        return t
    else: 
    	#This means that there is a line intersection but not a ray intersection.
        return 0
###

def plot_mesh(verts, faces):
    # plot the ellipsoid (or sphere)
    fig = plt.figure(figsize=(8,8))
    ax = mplot3d.Axes3D(fig)
    meshvectors = mesh_vectors(verts, faces)
    
    ax.add_collection3d(mplot3d.art3d.Poly3DCollection(meshvectors, 
    facecolor = [0.5,0.5,0.5], lw = 0.5, edgecolor = [0,0,0], alpha = .8, antialiaseds = True))

    scale = verts.flatten('F')
    ax.auto_scale_xyz(scale, scale, scale)
###
    
def rangeDoppler(a, C, normal, centroid, msh, e_obs):
    # Number of orientations and facets
	Nf = np.shape(centroid)[0]

    # Range dimension scalers
	RanSc = x0/a
	
	# Initiate for the R-D map and sigma0, projected area and number of visible faces
	xsec = 0
	AreaSum = 0
	if PrintAll == True:
		A_r = 0
		Nvis = 0
	
	# Integrate the brightness over all triangles when observing in a phase angle alpha
	for i, centr in enumerate(centroid):
		nu = linalg.norm(normal[i])
		nHat = normal[i] / nu

		# Define the cosines of the incident and emergence angles using dot products
		u_0 = np.dot(nHat, e_obs)

		# For triangles that are illuminated and visible, find the reflectivity
		# and area, and sum them over the illuminated surface
		if (u_0 > 0):
			# Distance between the reference plane and the centroid
			dcei = np.dot(centr, e_obs)
			t0 = 1.1*a - dcei

			# Shadowing (not working properly)
			if Shadowing == True:
				isect = False
				for j in range(Nf):
					# Go through all other triangles to check if 
					# another triangle is shadowing Triangle i
					dcej = np.dot(centroid[j,:], e_obs)
					if j != i and dcei < dcej and dcej < 0.1:
						tj = RayTracer(msh[j], centroid[j,:], e_obs)
						if tj > 0:
							isect = True
							print('%d blocked %d, t1,t2: %.3f, %.3f' % (j,i,tj,t0))
							break
				if isect == True:
					continue
					
			# Continue from here if t0 is the shortest distance to the reference plane
			# Range: Reference range x0 - dot product of a triangle center and observation dir
			# np.dot(centr,e_obs) = centr[0] when e_obs = (1,0,0)
			# Leading edge distance is 0.1*x0 pixels
			ra = int(round(RanSc * t0)) 
			if ra > 1.5*x0:
			    continue  
			
			# Incidence angle
			th_i = np.arccos(u_0)
			# Area
			dA = 0.5 * nu			
			
			# Backscatter cross section of an element
			if law == 'Cos':
				BSC = (C+1) * u_0 ** (2*C)    						 # Cosine law
			elif law == 'Hag':
				BSC = 0.5 * C * (u_0 ** 4 + C - C * u_0**2)**(-1.5) # Hagfors          	
			
			# Reflectivity and area sums
			xsec += dA * BSC
			AreaSum += dA * u_0
			if PrintAll == True:
				A_r += dA			
				Nvis += 1 
	
	if PrintAll == True:	   
		print('%s elements of %s are visible' % (Nvis, Nf))
		print('Effective diameter is %.3f' % (2*(AreaSum/np.pi)*0.5))
		print('Total area of illuminated hemisphere is %.3f' % (A_r))
		#print('Average vertex side length is %.3f' % (4.0*(A_r/Nvis/3)**0.5))
		print('Total projected area is %.3f' % (AreaSum))
		print('Geometric gain is %.3f' % (xsec/AreaSum))
		#print('Range resolution is %.4f' % (a/x0))
    
	return xsec, AreaSum
###

def ShepardDensity(ralb):
    if ralb > 0.1:
        return (ralb+0.156)/0.144
    else:
        return 6.4 * np.arctanh(np.sqrt(ralb/1.2))
###

def main(delta):
	
    verts, faces = readObj(modelname)
    faces -= 1
#       verts /= np.max(verts)
    xd = np.max(verts[:,0])-np.min(verts[:,0])
    yd = np.max(verts[:,1])-np.min(verts[:,1])
    zd = np.max(verts[:,2])-np.min(verts[:,2])
    if a0 == 0 and b0 == 0 and c0 == 0:
        a = xd/2.0
        b = yd/2.0
        c = zd/2.0  
    elif a0 > 0 and b0 == 0 and c0 == 0:
        scal = a0/xd
        a = xd*scal/2.0
        b = yd*scal/2.0
        c = zd*scal/2.0
        verts *= scal
    else:
        a = a0/2.0
        b = b0/2.0
        c = c0/2.0
        verts[:,0] *= a0/xd
        verts[:,1] *= b0/yd
        verts[:,2] *= c0/zd            
    q, q2 = a/c, a/b	
    
    # Find the normals and center points of each facet
    normal, centroid, faces = getNormals(verts, faces)
		    
    # Plot the object
#     plot_mesh(verts, faces)

    # Print the input
    print(modelname)
    print('Dimensions: %.4f, %.4f, %.4f' % (2*a,2*b,2*c))
    print('Axis ratios: %.2f, %.2f' % (q, q2))
    print('%d vertices, %d faces' % (len(verts), len(faces)) )

    # Subradar point (See Shapeintro.pdf, 2.3. Spin properties)
    # Subradarpoint components (!!!difference in lat, lon w.r.t. the pole!!!)
    deg2rad = np.pi/180.0
    if delta == False:
        dB = (PoleEcLat + ObsEcLat) * deg2rad
        dL = PoleEcLon - ObsEcLon
        if dL >= 180.0:
            dL = (dL - 180.0) * deg2rad
        else:
            dL = (dL + 180.0) * deg2rad
        cosdelta = np.cos(dB) * np.cos(dL)							# [-1,1]
        delta = np.arccos(cosdelta)	
        delta = 90.0 - delta/deg2rad   # [0°,180°]
        # Print w.r.t. the equator
        print('Subradar latitude (delta): %.2f' % (delta))  
     
    else: 
        cosdelta = np.cos((90.0 - delta) * deg2rad)
        delta = np.arccos(cosdelta)/deg2rad
        print('Subradar latitude (delta): %.2f' % (delta))

    BperD = 8.0 * 3.1416 * np.cos(delta*deg2rad) / (wl*P)
    print('Doppler bandwidth: %.2f-%.2f' 
    % (BperD * b, BperD * a ))
  
    # Observer is in an angle delta from the target's spin pole (sin, 0, cos)
    sindelta = (1.0 - cosdelta*cosdelta) ** 0.5						# [0,1]
    e_obs = np.array([sindelta, 0.0, cosdelta])   
         
    # Find the normals and center points of each facet
    normal, centroid, faces = getNormals(verts, faces)
            
    # Range-Doppler map
    msh = mesh_vectors(verts, faces)
    xsec1, AreaSum1 = rangeDoppler(a, C, normal, centroid, msh, e_obs)
    
    RF = rcs_obs / xsec1
    
    angles = np.arange(0.0, 360.0, step_ang)
    N_ang = len(angles)
    rcsarray = np.zeros((N_ang,3))
    
    # Orientation-averaging   
    ang = angles[1]-angles[0]   
    rcsarray[0,0] = RF*xsec1
    rcsarray[0,1] = RF*rcsarray[0,0]/AreaSum1
    rcsarray[0,2] = AreaSum1
    print('Phase 0°, Projected area %.3f \n' % AreaSum1)
    
    for t in range(N_ang-1):
        if t == 0:
            Rz = rotellipsoid(ang * deg2rad, 2)
        for i,coord in enumerate(verts):
            verts[i,:] = np.dot(Rz,verts[i,:])
            
        normal, centroid, faces = getNormals(verts, faces)
        msh = mesh_vectors(verts, faces)
        xsec2, AreaSum2 = rangeDoppler(a, C, normal, centroid, msh, e_obs)

        # Per scan
        rcsarray[t+1,0] = RF*xsec2
        rcsarray[t+1,1] = RF*xsec2/AreaSum2
        rcsarray[t+1,2] = AreaSum2
        
        print('Phase %d°, Projected area %.3f \n' % (angles[t+1], AreaSum2))

    fixR = rcs_obs / np.mean(rcsarray[:,0])
    RF *= fixR
    rcsarray *= fixR
    eps = ((1.0+RF**0.5)/(1.0-RF**0.5))**2
    
    # Mean and other stats
    print()
    print('Statistics for %d orientations:' % N_ang)
    print('Mean radar cross section (OC, km2): %.3f' % np.mean(rcsarray[:,0]))
    print('Full range: %.3f -> %.3f' % (np.min(rcsarray[:,0]), np.max(rcsarray[:,0])))
    print('Mean radar albedo (OC): %.3f' % np.mean(rcsarray[:,1]))
    print('Full range: %.3f -> %.3f' % (np.min(rcsarray[:,1]), np.max(rcsarray[:,1])))
    print('R parameter is: %.3f' % RF)
    print('C parameter is: %.3f' % C)
    print('Electric permittivity based on R parameter: %.3f' % eps)
    print('Near-surface bulk density based on permittivity (Hickson et al. 2018): %.2f g cm^-3' 
    % ((eps**0.333-1.0)/0.307))
    print('Near-surface bulk density based on radar albedos (Shepard et al. 2010): %.2f-%.2f g cm^-3'
    % (ShepardDensity(np.min(rcsarray[:,1])), ShepardDensity(np.max(rcsarray[:,1]))))
    
    if MakePlots == True:
        # Radar cross sections and albedos as a function of orientation.
        figP, axP = plt.subplots(3, sharex=True)
        figP.subplots_adjust(left=0.12, bottom=0.11, top=0.95, right=0.95, hspace=0.08)
    
        axP[0].set_ylabel(r"$\sigma$ (km$^2$)", fontsize=22)
        axP[0].tick_params(axis='both',labelsize=22)
        axP[1].set_ylabel(r"$A_{proj}$ (km$^2$)", fontsize=22)
        axP[1].tick_params(axis='both',labelsize=22)
        axP[2].set_ylabel(r"$\hat{\sigma}$", fontsize=22)
        axP[2].set_xlabel("Orientation (°)", fontsize=22)
        axP[2].tick_params(axis='both',labelsize=22)  
    
        axP[0].plot(angles, rcsarray[:,0], 'k-', lw=2)
        axP[1].plot(angles, rcsarray[:,2], 'k-', lw=2)
        axP[2].plot(angles, rcsarray[:,1], 'k-', lw=2)  
        
        plt.show() 
###

if __name__=="__main__":
    main(delta)