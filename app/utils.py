import os


def get_R_script_base_path():
    return f"{os.getcwd()}/../app" if '/tests' in os.getcwd() else f"{os.getcwd()}/app"


def get_directory_files(directory):
    return [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
