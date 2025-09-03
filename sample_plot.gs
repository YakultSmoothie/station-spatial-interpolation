'reinit'
rc = gsfallow('on')

'sdfopen etopo1_ice_g_f4-slrain.nc'
'define hgt1 = hgt.1'
'close 1'


'open interpolated_PP01.ctl'
'define hgtlt = lterp(hgt1, lat)'
'define landmask = maskout(1, hgtlt - 0)'

'set xlint 1'
'set ylint 1'
'set mpdset hires'
'd sum(r1, t=1, t=3) * landmask'

'draw title Accumulated \ Precipitation (mm)' 

