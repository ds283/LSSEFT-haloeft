__author__ = 'ds283'

import os
import tarfile
import subprocess


def list_files(path):
    # returns a list of names (with extension, without full path) of all files
    # in folder path
    files = []
    for name in os.listdir(path):

        if os.path.isfile(os.path.join(path, name)):

            name_lower = name.lower()

            if name_lower[0] != '.' \
                and not name_lower.endswith('.pyc') \
                and name_lower != 'r01.ini' \
                and name_lower != 'r02.ini' \
                and name_lower != 'r03.ini' \
                and name_lower != 'r04.ini' \
                and name_lower != 'r05.ini' \
                and name_lower != 'r06.ini' \
                and name_lower != 'r07.ini' \
                and name_lower != 'r08.ini' \
                and name_lower != 'r09.ini' \
                and name_lower != 'r10.ini' \
                and name_lower != 'package_tool.py':

                    files.append(name)

    return files


def list_directories(path):
    dirs = []
    for name in os.listdir(path):
        if os.path.isdir(os.path.join(path, name)):
            dirs.append(name)
    return dirs


def add_folder(tree_path, archive_path, archive):
    # get list of files at this level
    files = list_files(tree_path)

    # get list of directories at this level
    dirs = list_directories(tree_path)

    for file in files:

        file_lower = file.lower()

        if file_lower != '.DS_Store':

            tree_name = os.path.join(tree_path, file)
            archive_name = os.path.join(archive_path, file)

            archive.add(name=tree_name, arcname=archive_name, recursive=False)

    for dir in dirs:

        dir_lower = dir.lower()

        if dir[0] != '.' \
            and dir_lower != '.idea' \
            and dir_lower != 'data' \
            and dir_lower != 'theory'\
            and dir_lower != 'scripts' \
            and dir_lower != 'packages' :

            tree_name = os.path.join(tree_path, dir)
            archive_name = os.path.join(archive_path, dir)

            # add this directory to the archive
            archive.add(name=tree_name, arcname=archive_name, recursive=False)

            # and recursively add its contents
            add_folder(tree_name, archive_name, archive)


def package_LSSEFT_haloeft(package_dir, archive_file, version_string):

    # tree_path points to the position of some object within the source tree
    tree_path = os.path.join("..")

    # archive path is the corresponding position of the object within the tar archive
    archive_path = "LSSEFT-haloeft" + "_" + version_string

    # write out all variables
    tarname = os.path.join(package_dir, archive_file)
    with tarfile.open(name=tarname, mode='w:gz', format=tarfile.PAX_FORMAT) as archive:

        # add top-level directory to archive
        archive.add(name=tree_path, arcname=archive_path, recursive=False)

        add_folder(tree_path, archive_path, archive)

        archive.close()


version = "2018_01_Apollo"

cwd = os.getcwd()
package_dir = os.path.join(cwd, "..", "packages")
version_dir = os.path.join(package_dir, version)

# ensure packaging directory exists
if not os.path.exists(version_dir):
    os.makedirs(version_dir)

package_LSSEFT_haloeft(version_dir, "LSSEFT-haloeft" + "_" + version + ".tar.gz", version)
