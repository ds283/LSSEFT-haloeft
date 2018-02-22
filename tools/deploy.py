import os
import version as release

# DEPLOYMENT PARAMETERS

release_name = "LSSEFT-haloeft_{tag}".format(tag=release.version)
deploy_root = os.path.join("LSSEFT-haloeft")
local_site = None
cosmosis_executable = "./bin/cosmosis"
sampler = "maxlike"
MPI_processes = 8
email = "D.Seery@sussex.ac.uk"
