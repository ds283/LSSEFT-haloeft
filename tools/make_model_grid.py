import os

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

root_path = os.path.join("..", "models")

for b in bias_models:

    for r in RSD_models:

        folder = os.path.join(root_path, bias_folders[b], RSD_folders[r])

        if not os.path.exists(folder):
            print 'error -- did not find expected model folder {p}'.format(p=folder)
            raise RuntimeError

        output_folder = os.path.join(folder, "output")

        if not os.path.exists(folder):

            try:
                os.makedirs(folder)
            except OSError, e:
                if e.errno != os.errno.EEXIST:
                    raise

        for reg in regions:

            ini_file = os.path.join(folder, inis[r])
            common_file = os.path.join("LSSEFT-haloeft", "haloeft_common.ini")
            config_file = os.path.join("LSSEFT-haloeft", bias_folders[b], RSD_folders[r], "config.ini")
            output_file = os.path.join("LSSEFT-haloeft", bias_folders[b], RSD_folders[r], "output", outputs[r])

            with open(ini_file, "w") as f:

                f.write("[runtime]")
                f.write("sampler = emcee")
                f.write("")
                f.write("%include {p}".format(p=common_file))
                f.write("%include {p}".format(p=config_file))
                f.write("")
                f.write("[output]")
                f.write("filename = {p}".format(p=config_file))
                f.write("")
                f.write("[HaloEFT]")
                f.write("realization = {n}".format(n=numbers[r]))
