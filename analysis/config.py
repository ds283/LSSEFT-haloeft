import os


def make_config_block(realization, use_Halofit, fit_kmin=0.02, fit_kmax=0.30, ren_kmin=0.02, ren_kmax=0.30):

    # use a dictionary to mimic the CosmoSIS datablock API
    my_config = {}

    root = os.path.join("..", "assets")

    my_config["HaloEFT", "h1_means"] = os.path.join(root, "data/means/pkpole_wizcola_1hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h3_means"] = os.path.join(root, "data/means/pkpole_wizcola_3hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h9_means"] = os.path.join(root, "data/means/pkpole_wizcola_9hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h11_means"] = os.path.join(root, "data/means/pkpole_wizcola_11hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h15_means"] = os.path.join(root, "data/means/pkpole_wizcola_15hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h22_means"] = os.path.join(root, "data/means/pkpole_wizcola_22hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h1_matrix"] = os.path.join(root, "data/covariances/covar_1hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h3_matrix"] = os.path.join(root, "data/covariances/covar_3hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h9_matrix"] = os.path.join(root, "data/covariances/covar_9hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h11_matrix"] = os.path.join(root, "data/covariances/covar_11hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h15_matrix"] = os.path.join(root, "data/covariances/covar_15hr_z0pt2_0pt6.dat")
    my_config["HaloEFT", "h22_matrix"] = os.path.join(root, "data/covariances/covar_22hr_z0pt2_0pt6.dat")

    my_config["HaloEFT", "realization"] = realization

    if use_Halofit:
        my_config["HaloEFT", "theory_db"] = os.path.join(root, "theory/WizCOLA_HALOFIT_halo@z=0_kWiggleZ.sqlite")
    else:
        my_config["HaloEFT", "theory_db"] = os.path.join(root, "theory/WizCOLA_CAMB_halo_full2@z=0_kWiggleZ.sqlite")

    my_config["HaloEFT", "model"] = 0
    my_config["HaloEFT", "growth_params"] = 0
    my_config["HaloEFT", "loop_params"] = 0
    my_config["HaloEFT", "XY_params"] = 0
    my_config["HaloEFT", "zid"] = 0
    my_config["HaloEFT", "init_Pk"] = 0
    my_config["HaloEFT", "final_Pk"] = 1
    my_config["HaloEFT", "IR_cutoff"] = 0
    my_config["HaloEFT", "UV_cutoff"] = 0
    my_config["HaloEFT", "IR_resum"] = 0
    my_config["HaloEFT", "fit_kmin"] = fit_kmin
    my_config["HaloEFT", "fit_kmax"] = fit_kmax
    my_config["HaloEFT", "renormalize_kmin"] = ren_kmin
    my_config["HaloEFT", "renormalize_kmax"] = ren_kmax

    return my_config
