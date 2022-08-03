from appdata import AppDataPaths; _app_paths = AppDataPaths("mhd")
from pathlib import Path
mhd_settings_filepath = Path(_app_paths.app_data_path) / "mhd_settings.ini"


import configparser; _conf = configparser.ConfigParser(); _conf.read( mhd_settings_filepath )

import sys
MDSplus_python_path = _conf.get("mdsplus", "python_path") 
sys.path.append(MDSplus_python_path) # for Windows users, MDSplus is normally installed at "C:\\Program Files\\MDSplus\\python"

MDSplus_server_IP = _conf.get("mdsplus", "server_ip")