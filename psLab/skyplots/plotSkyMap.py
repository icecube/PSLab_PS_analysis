from SkyMap import *

s = SkyMap()

#s.AddHistogramFile("/net/user/jdumm/testlab/macro_llh/ic40_fix/AllSkyBasic_FineSkyMap.root")
s.AddHistogramFile("/net/user/mfbaker/psLab/trunk/macro_llh/ic59/maps/flare/mapIC4059_unblind.root")

s.AddGalPlane("resources/GalPlane3600.coords")

s.set_latitude_grid(4)
s.Plot("hammer")

savefig("hammer_map.png")

show()



