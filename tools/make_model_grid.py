import os
import stat

bias_models = [ "coevolution", "full", "linear", "MR" ]
RSD_models = [ "EFT", "KaiserHalofit", "KaiserTree", "ZhengSong" ]

bias_folders = { "coevolution": "coevolution", "full": "full", "linear": "linear", "MR": "MR" }
RSD_folders = { "EFT": "EFT", "KaiserHalofit": "KaiserHalofit", "KaiserTree": "KaiserTree", "ZhengSong": "ZS" }
RSD_config_files = { "EFT": "EFT_common.ini", "KaiserHalofit": "KaiserHalofit_common.ini", "KaiserTree": "KaiserTree_common.ini", "ZhengSong": "ZhengSong_common.ini"}

regions = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10']

numbers = {'r01': 1, 'r02': 2, 'r03': 3, 'r04': 4, 'r05': 5, 'r06': 6, 'r07': 7, 'r08': 8, 'r09': 9, 'r10': 10}

inis = {'r01': 'r01.ini', 'r02': 'r02.ini', 'r03': 'r03.ini',
        'r04': 'r04.ini', 'r05': 'r05.ini', 'r06': 'r06.ini',
        'r07': 'r07.ini', 'r08': 'r08.ini', 'r09': 'r09.ini',
        'r10': 'r10.ini'}

outputs = {'r01': 'r01.txt', 'r02': 'r02.txt', 'r03': 'r03.txt',
          'r04': 'r04.txt', 'r05': 'r05.txt', 'r06': 'r06.txt',
          'r07': 'r07.txt', 'r08': 'r08.txt', 'r09': 'r09.txt',
          'r10': 'r10.txt'}


# DEPLOYMENT PARAMETERS

deploy_root = os.path.join("LSSEFT-haloeft")
cosmosis_executable = "./bin/cosmosis"
MPI_processes = 8


# build global paths

local_models_root = os.path.join("..", "models")
deploy_models_root = os.path.join(deploy_root, "models")

local_scripts_root = os.path.join("..", "scripts")


def populate_region_ini_files():

    for b in bias_models:

        for r in RSD_models:

            local_folder = os.path.join(local_models_root, bias_folders[b], RSD_folders[r])

            if not os.path.exists(local_folder):
                print 'error -- did not find expected model local_folder {p}'.format(p=local_folder)
                raise RuntimeError

            output_folder = os.path.join(local_folder, "output")

            if not os.path.exists(output_folder):

                try:
                    os.makedirs(output_folder)
                except OSError, e:
                    if e.errno != os.errno.EEXIST:
                        raise

            for reg in regions:

                ini_file = os.path.join(local_folder, inis[reg])
                common_file = os.path.join(deploy_root, "haloeft_common.ini")
                config_file = os.path.join(deploy_models_root, bias_folders[b], RSD_folders[r], "config.ini")
                output_file = os.path.join(deploy_models_root, bias_folders[b], RSD_folders[r], "output", outputs[reg])

                with open(ini_file, "w") as f:

                    f.write("[runtime]\n")
                    f.write("sampler = emcee\n")
                    f.write("\n")
                    f.write("%include {p}\n".format(p=common_file))
                    f.write("%include {p}\n".format(p=config_file))
                    f.write("\n")
                    f.write("[output]\n")
                    f.write("filename = {p}\n".format(p=output_file))
                    f.write("\n")
                    f.write("[HaloEFT]\n")
                    f.write("realization = {n}\n".format(n=numbers[reg]))


def populate_RSD_config_files():

    for b in bias_models:

        for r in RSD_models:

            bias_folder = os.path.join(deploy_models_root, bias_folders[b])
            RSD_folder = os.path.join(deploy_models_root, bias_folders[b], RSD_folders[r])

            local_folder = os.path.join(local_models_root, bias_folders[b], RSD_folders[r])

            if not os.path.exists(local_folder):
                print 'error -- did not find expected model local_folder {p}'.format(p=local_folder)
                raise RuntimeError

            config_file = os.path.join(local_folder, "config.ini")

            with open(config_file, "w") as f:

                common_config = os.path.join(deploy_root, RSD_config_files[r])
                bias_values_file = os.path.join(bias_folder, "parameters.ini")
                pipeline_file = os.path.join(RSD_folder, "pipeline.py")

                f.write("%include {common}\n".format(common=common_config))
                f.write("\n")
                f.write("[pipeline]\n")
                f.write("values = {file}\n".format(file=bias_values_file))
                f.write("\n")
                f.write("[HaloEFT]\n")
                f.write("file = {file}\n".format(file=pipeline_file))


def make_scripts_directory():

    if not os.path.exists(local_scripts_root):

        try:
            os.makedirs(local_scripts_root)
        except OSError, e:
            if e.errno != os.errno.EEXIST:
                raise



def populate_batch_scripts():

    # now generate shell scripts to do batch jobs -- organized by bias model:
    for b in bias_models:

        script = os.path.join(local_scripts_root, "run-{b}".format(b=b))

        with open(script, "w") as f:

            for r in RSD_models:

                for reg in regions:

                    ini_file = os.path.join(deploy_models_root, bias_folders[b], RSD_folders[r], inis[reg])
                    f.write("mpiexec -n {np} {cosmosis} --mpi {config}\n".format(np=MPI_processes, cosmosis=cosmosis_executable, config=ini_file))

        st = os.stat(script)
        os.chmod(script, st.st_mode | stat.S_IEXEC)


    # now generate shell scripts to do batch jobs -- organized by RSD model:
    for r in RSD_models:

        script = os.path.join(local_scripts_root, "run-{r}".format(r=r))

        with open(script, "w") as f:

            for b in bias_models:

                for reg in regions:

                    ini_file = os.path.join(deploy_models_root, bias_folders[b], RSD_folders[r], inis[reg])
                    f.write("mpiexec -n {np} {cosmosis} --mpi {config}\n".format(np=MPI_processes, cosmosis=cosmosis_executable, config=ini_file))

        st = os.stat(script)
        os.chmod(script, st.st_mode | stat.S_IEXEC)


    # individual grid cell
    for b in bias_models:

        for r in RSD_models:

            script = os.path.join(local_scripts_root, "run-{b}-{r}".format(b=b, r=r))

            with open(script, "w") as f:

                for reg in regions:

                    ini_file = os.path.join(deploy_models_root, bias_folders[b], RSD_folders[r], inis[reg])
                    f.write("mpiexec -n {np} {cosmosis} --mpi {config}\n".format(np=MPI_processes, cosmosis=cosmosis_executable, config=ini_file))

            st = os.stat(script)
            os.chmod(script, st.st_mode | stat.S_IEXEC)


    # now generate shell scripts to do batch jobs -- single script to run entire model grid
    script = os.path.join(local_scripts_root, "run-grid")

    with open(script, "w") as f:

        for b in bias_models:

            for r in RSD_models:

                for reg in regions:

                    ini_file = os.path.join(deploy_models_root, bias_folders[b], RSD_folders[r], inis[reg])
                    f.write("mpiexec -n {np} {cosmosis} --mpi {config}\n".format(np=MPI_processes, cosmosis=cosmosis_executable, config=ini_file))

    st = os.stat(script)
    os.chmod(script, st.st_mode | stat.S_IEXEC)


populate_region_ini_files()

populate_RSD_config_files()

make_scripts_directory()
populate_batch_scripts()
