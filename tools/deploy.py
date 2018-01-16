import os
import version as release

# DEPLOYMENT PARAMETERS

release_name = "LSSEFT-haloeft_{tag}".format(tag=release.version)
deploy_root = os.path.join("/", "home", "d", "ds", "ds283", "LSSEFT-haloeft", release_name)
cosmosis_executable = "cosmosis"
MPI_processes = 64
email = "D.Seery@sussex.ac.uk"
