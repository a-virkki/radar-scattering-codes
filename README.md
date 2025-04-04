# radar-scattering-codes
This repository includes a Python code for simulating the radar scattering for specific shape models or a generic ellipsoid. It can provide the "R" value for the scattering law better than SHAPE, interprets it in terms of the near-surface density, and shows the range of variation in the radar cross section and radar albedo over a rotation. You can provide either the subradar latitude or the set of observer and polar ecliptic coordinates (make delta False for the latter! Note: 0 = False so use e.g. 0.01 for an equatorial view!)
The required input includes:

modelname = 'shapemodel.obj'  # OBJ-format models can be found at https://sbn.psi.edu/pds/shape-models/

a0 = 14.0               # Size: a0 is the longest axis dimension. Use 0 to not scale the model size.

b0, c0 = 0, 0           # Intermediate and short axis dimensions. Use 0 to not scale the model size.

delta = 0					      # subradar latitude ([-90,90])

ObsEcLat = 35.0         # Asteroid's observer-centered ecliptic latitude

ObsEcLon = 324.0        # Asteroid's observer-centered ecliptic longitude

PoleEcLat = -44.3 			# Pole ecliptic latitude in degrees

PoleEcLon = 35.6 				# Pole ecliptic longitude in degrees

x0 = 10 						    # Map dimension (number of pixels per side)

law = 'Cos'						  # Scattering law: 'Cos'/'Hag' for Cosine/Hagfors

C = 0.77 						    # scattering law parameter (roughness)

rcs_obs = 25					  # The observed radar cross section

step_ang = 45					  # Angle interval (e.g., 90 for 4 orientations, 45 for 8 orientations)

P = 321.1               # Rotation period in seconds

wl = 0.126              # Observation wavelength

Shadowing = False				# Shadowing (recommended only for significantly concave shapes)

PrintAll = False				# Print all information for every orientation

MakePlots = False				# Plot the cross sections and albedos for every orientation
