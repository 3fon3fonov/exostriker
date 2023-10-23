import desktop_file as df
from pathlib import Path
import sys
import os


def get_exo_dir():
    def anaconda_error():
        sys.exit("Error: anaconda environment detected but path not found.")

    pre = sys.prefix
    if "anaconda" in pre:
        print("Anaconda detected.")

        pre = Path(pre).parts

        ind_anaconda = -1
        for i in range(len(pre)):
            if "anaconda" in pre[i]:
                ind_anaconda = i
                break

        if ind_anaconda == -1:
            anaconda_error()

        env_name, anaconda_path = "", ""
        try:
            env_name = pre[ind_anaconda + 2]
            anaconda_path = pre[:ind_anaconda + 1]
        except IndexError:
            anaconda_error()

        anaconda_path = Path(*anaconda_path).joinpath("bin", "conda")

        conda_activate = "{} run --no-capture-output --name {} ".format(anaconda_path, env_name)
    else:
        conda_activate = ""

    path = Path(__file__).parts[:-1]
    path = Path(path[0]).joinpath(*path[1:])

    return conda_activate, str(path)


def create():
    def module_create(path):
        print("Creating shortcut at: ", path)
        shortcut = df.Shortcut(path, "Exostriker", "{}python{} -m exostriker".format(conda, vers))
        shortcut.setTitle("Exostriker")
        shortcut.setWorkingDirectory(dir)
        shortcut.setComment("Transit and Radial velocity Interactive Fitting tool for Orbital analysis and N-body "
                            "simulations")
        shortcut.setIcon(icon_path)
        shortcut.setCategories("Science;X-Astrophysics;X-Astronomy;X-Education;")
        shortcut.save()

    conda, dir = get_exo_dir()
    icon_path = str(Path(dir).joinpath("source", "png", "33_striker.ico"))

    # Get the actual Python version
    vers = str(sys.version_info.major) + "." + str(sys.version_info.minor)

    if "linux" in sys.platform:
        module_create(df.getMenuPath())
    elif "win" in sys.platform:
        import win32com.client

        fls = os.listdir(dir)
        for fl in fls:
            if "Exostriker" in fl and ".bat" in fl:
                os.remove(str(Path(dir).joinpath(fl)))

        module_create(dir)
        os.remove(str(Path(dir).joinpath("Exostriker.lnk")))

        file = open(str(Path(dir).joinpath("Exostriker.vbs")), "w")
        execFile = "Set WshShell = CreateObject(\"WScript.Shell\")\n" \
                   "WshShell.Run chr(34) & \"{}\" & Chr(34), 0\n" \
                   "Set WshShell = Nothing".format(str(Path(dir).joinpath("Exostriker.bat")))
        file.write(execFile)
        file.close()

        shell = win32com.client.Dispatch("WScript.Shell")
        shortcut = shell.CreateShortcut(str(Path(dir).joinpath("Exostriker.lnk")))
        shortcut.TargetPath = str(Path(dir).joinpath("Exostriker.vbs"))
        shortcut.IconLocation = icon_path
        shortcut.Save()

        src = str(Path(dir).joinpath("Exostriker.lnk"))
        print("Creating shortcut at: ", df.getDesktopPath())
        os.system("copy {} \"{}\"".format(src, df.getDesktopPath()))
        print("Creating shortcut at: ", df.getMenuPath())
        os.system("copy {} \"{}\"".format(src, df.getMenuPath()))

    print("No errors reported. Logout and login again for changes to take effect.")


def remove():
    def erasing(path_fl, fl):
        path_fl = str(Path(path_fl).joinpath(fl))
        print("Removing file: ", path_fl)

        try:
            os.remove(path_fl)
            print("Done.")
        except FileNotFoundError:
            print("File {} not found.".format(path_fl))

    if "linux" in sys.platform:
        erasing(df.getMenuPath(), "Exostriker.desktop")
    elif "win" in sys.platform:
        path_to_fl = Path(__file__).parts[:-1]
        path_to_fl = Path(path_to_fl[0]).joinpath(*path_to_fl[1:])

        erasing(df.getDesktopPath(), "Exostriker.lnk")
        erasing(df.getMenuPath(), "Exostriker.lnk")
        erasing(path_to_fl, "Exostriker.vbs")
        erasing(path_to_fl, "Exostriker.lnk")
        erasing(path_to_fl, "Exostriker.bat")
