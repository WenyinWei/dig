from appdata import AppDataPaths; _app_paths = AppDataPaths("mhd")
from pathlib import Path
mhd_settings_filepath = Path(_app_paths) / "mhd_settings.ini"

import configparser; _conf = configparser.ConfigParser(); _conf.read( mhd_settings_filepath )
MDSplus_server_IP = _conf.get("mdsplus", "server_ip")