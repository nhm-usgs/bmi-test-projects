
from pymt.models import PRMSSurface, PRMSSoil
from pathlib import Path

run_dir = '../prms/pipestem'
config_surf = 'control_surface.simple1'
config_soil = 'control_soil.simple1'
print(Path(run_dir).exists())
print((Path(run_dir) / config_surf).exists())
print((Path(run_dir) / config_soil).exists())

msurf = PRMSSurface()
msoil = PRMSSoil()

print(msoil.name)

msurf.initialize(config_surf, run_dir)
msoil.initialize(config_soil, run_dir)

msurf.update()
msoil.update()

msoil.finalize()
msurf.finalize()

tmp = 0

