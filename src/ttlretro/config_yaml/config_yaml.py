from __future__ import division, unicode_literals

import yaml

class Config(object):
    def __init__(self, dpth):
        with open(dpth, "r", encoding="utf-8") as fp:
            conf_dict = yaml.load(fp, Loader=yaml.FullLoader)
            for key, val in conf_dict.items():
                setattr(self, key, val)




