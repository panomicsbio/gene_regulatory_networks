import os


def get_R_script_base_path():
    return f"{os.getcwd()}/../app" if '/tests' in os.getcwd() else f"{os.getcwd()}/app"
