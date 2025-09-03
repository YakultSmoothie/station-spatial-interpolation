dset ^./output/nc/interpolated_PP01_%y4%m2%d2%h2.nc
dtype netcdf
options template

title station_spatial_interpolation.py out
undef NaN

xdef 72 linear 120 0.03
ydef 118 linear 21.9 0.03
zdef 1 linear 1000 1
tdef 999 linear 01z30JUL2025 60mn

vars 1
   PP01=>r1  1  y,x  [mm/h]
endvars
