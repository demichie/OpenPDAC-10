path = '../'
DEM_file = 'str_5m.asc'

resample = 2  # subsampling factor of the original DEM
dist0 = 25.0  # semi-width of the fissure/radius of the cylinder
enne = 3.0  # shape parameter (1 linear, 2 spherical)
depth = 0.0  # depth of the fissure
nlevels = 4  # levels of refinements of the subsampled grid

xb = 0.0  # horizontal x-translation of the base of the fissure/conduit
yb = 0.0  # horizontal y-translation of the base of the fissure/conduit

xc = 518307  # x of first point of fissure/x-center of cylinder (UTM)
yc = 4293852  # y of first point of fissure/y-center of cylinder (UTM)

# FOR CYLINDRICAL FISSURE
points = [(xc, yc), (xc, yc)]

# FOR LINEAR FISSURE
# points = [(xc, yc), (xc+20, yc+10), (xc+40, yc+15)]

saveDicts_flag = True
z_atm = 2000
offset_mesh = 50.0
delta_mesh = 200.0
domain_size_x = 1500.0
domain_size_y = 1500.0


