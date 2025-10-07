import os, re, time, inspect
from functools import wraps
from importlib.metadata import version


global LOG_DIR
global VERBOSE
global session_init

LOG_DIR = None
VERBOSE = True
session_init = False


def init_log_dir():
    global LOG_DIR

    if LOG_DIR == None:
        LOG_DIR = os.path.join(os.path.expanduser("~"), ".CALINS")
    if not os.path.exists(LOG_DIR):
        os.mkdir(LOG_DIR)

    check_folder_space()


def check_folder_space():
    max_files_nb = 9
    log_files_path = []

    for root, dirs, files in os.walk(LOG_DIR):
        for file in files:
            if file.endswith(".log"):
                log_files_path.append(os.path.join(root, file))

    len_log_files = len(log_files_path)
    if len_log_files > max_files_nb:
        warn(f"Log directory {LOG_DIR} contains too many files ({len_log_files} > max log file nb ({max_files_nb})), removing oldest log files.")
    while len_log_files > max_files_nb:
        oldest_file = min(log_files_path, key=os.path.getctime)
        log_files_path.remove(oldest_file)
        os.remove(oldest_file)
        len_log_files -= 1


def init_big_header(pkg_name):
    len_header = 90

    try:
        __version__ = version(pkg_name)
    except:
        __version__ = f"no-version-found"

    pkg_fullname = f"{pkg_name} v. {__version__}"
    first_line = f"{len_header*'_'}" + "\n" + f"|{(len_header-2)*' '}|"
    end_line = f"|{(len_header-2)*'_'}|"
    info_str = f"(package : {pkg_fullname}) - python session"
    big_header = "\n" + first_line + "\n" + "CALINS".center(len_header) + "\n" + info_str.center(len_header) + "\n" + end_line + "\n"

    return big_header


def write_in_log(msg, type=None):
    with open(os.path.join(LOG_DIR, f'CALINS_{time.strftime("%Y-%m-%d")}.log'), "a") as f:
        f.write(f"{time.strftime('%H:%M:%S')} - {type+' : ' if type is not None else ''}{msg}\n")


def print_in_terminal(msg, type=None):
    print(f"[CALINS] {type+' : ' if type is not None else ''}{msg}")


def write_and_print(msg, type=None, bypass_verbose=False):
    write_in_log(msg, type)

    if VERBOSE or bypass_verbose:
        print_in_terminal(msg, type)


def warn(msg):
    write_and_print(msg, "WARNING", bypass_verbose=True)


def init_session(pkg_name=None):
    init_log_dir()
    big_header = init_big_header(pkg_name=pkg_name)
    write_and_print(big_header, bypass_verbose=True)


def log_exec():
    def log_exec_deco(func):

        @wraps(func)
        def log_exec_wrapper(*args, **kwargs):
            global session_init

            if not session_init:
                pkg_name = func.__module__.split(".")[0]
                init_session(pkg_name=pkg_name)
                session_init = True

            func_args = inspect.signature(func).bind(*args, **kwargs).arguments

            if "__init__" in func.__qualname__ and ".classes" in func.__module__:
                tab_str = ""
                parent_str = f"{tab_str}{func.__qualname__.split('.')[0]} - "
            elif ".classes" in func.__module__:
                tab_str = 2 * "    "
                parent_str = f"{tab_str}{func.__qualname__.split('.')[0]} - "
            else:
                tab_str = ""
                parent_str = ""

            arg_str = ""
            for key, val in func_args.items():
                if re.search("path|case", key) and type(val) == str:
                    arg_str += f"\n{11*' '}{tab_str}    {key}={val}"
                elif "case" in key and val.__class__.__name__ == "Case":
                    arg_str += f"\n{11*' '}{tab_str}    {key}=(Case){val.sdf_path}"
                elif re.search("iso|reac|threshold|normalize|integrate|energy", key):
                    arg_str += f"\n{11*' '}{tab_str}    {key}={val}"

            if not "calcul_bias" in inspect.stack()[1][3]:
                write_and_print(f"{parent_str}{func.__name__} : {arg_str}")

            ret = func(*args, **kwargs)

            if "calcul" in func.__name__ and not "calcul" in inspect.stack()[1][3]:
                write_and_print(f"{tab_str}    -> {ret.value if re.search('Uncertainty|Bias',ret.__class__.__name__) else ret}")

            return ret

        return log_exec_wrapper

    return log_exec_deco
