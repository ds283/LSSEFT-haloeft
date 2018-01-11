import os
import stat

bias_models = [ "coevolution", "full", "linear", "MR" ]
RSD_models = [ "EFT", "KaiserHalofit", "KaiserTree", "ZhengSong" ]

bias_folders = { "coevolution": "coevolution", "full": "full", "linear": "linear", "MR": "MR" }
RSD_folders = { "EFT": "EFT", "KaiserHalofit": "KaiserHalofit", "KaiserTree": "KaiserTree", "ZhengSong": "ZS" }

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

local_models_root = os.path.join("..", "models")
deploy_common_root = os.path.join("/", "home", "d", "ds", "ds283", "LSSEFT-haloeft")
deploy_models_root = os.path.join(deploy_common_root, "models")

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
            common_file = os.path.join(deploy_common_root, "haloeft_common.ini")
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


local_scripts_root = os.path.join("..", "scripts")
deploy_model_root = os.path.join("/", "home", "d", "ds", "ds283", "LSSEFT-haloeft", "models")

if not os.path.exists(local_scripts_root):

    try:
        os.makedirs(local_scripts_root)
    except OSError, e:
        if e.errno != os.errno.EEXIST:
            raise


# now generate shell scripts to do batch jobs -- bias model:
for b in bias_models:

    script = os.path.join(local_scripts_root, "run-{b}".format(b=b))

    with open(script, "w") as f:

        for r in RSD_models:

            for reg in regions:

                ini_file = os.path.join(deploy_model_root, bias_folders[b], RSD_folders[r], inis[reg])
                f.write("mpiexec -n 8 ./bin/cosmosis --mpi {config}\n".format(config=ini_file))

    st = os.stat(script)
    os.chmod(script, st.st_mode | stat.S_IEXEC)


# now generate shell scripts to do batch jobs -- RSD model:
for r in RSD_models:

    script = os.path.join(local_scripts_root, "run-{r}".format(r=r))

    with open(script, "w") as f:

        for b in bias_models:

            for reg in regions:

                ini_file = os.path.join(deploy_model_root, bias_folders[b], RSD_folders[r], inis[reg])
                f.write("mpiexec -n 8 ./bin/cosmosis --mpi {config}\n".format(config=ini_file))

    st = os.stat(script)
    os.chmod(script, st.st_mode | stat.S_IEXEC)


# now generate shell scripts to do batch jobs -- entire model grid
script = os.path.join(local_scripts_root, "run-grid")

with open(script, "w") as f:

    for b in bias_models:

        for r in RSD_models:

            for reg in regions:

                ini_file = os.path.join(deploy_model_root, bias_folders[b], RSD_folders[r], inis[reg])
                f.write("mpiexec -n 8 ./bin/cosmosis --mpi {config}\n".format(config=ini_file))

st = os.stat(script)
os.chmod(script, st.st_mode | stat.S_IEXEC)


# individual grid cell
for b in bias_models:

    for r in RSD_models:

        script = os.path.join(local_scripts_root, "run-{b}-{r}".format(b=b, r=r))

        with open(script, "w") as f:

            for reg in regions:

                ini_file = os.path.join(deploy_model_root, bias_folders[b], RSD_folders[r], inis[reg])
                f.write("mpiexec -n 8 ./bin/cosmosis --mpi {config}\n".format(config=ini_file))

        st = os.stat(script)
        os.chmod(script, st.st_mode | stat.S_IEXEC)
