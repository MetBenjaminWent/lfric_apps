import re
import sys

from metomi.rose.upgrade import MacroUpgrade

from .version20_21 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro
class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>
    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"
    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""


class vn21_t657(MacroUpgrade):
    """Upgrade macro for ticket #657 by Ian Boutle."""

    BEFORE_TAG = "vn2.1"
    AFTER_TAG = "vn2.1_t657"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        nml = "namelist:cloud"
        fsd_nonconv = self.get_setting_value(config, [nml, "fsd_nonconv_const"])
        self.remove_setting(config, [nml, "fsd_nonconv_const"])
        self.add_setting(config, [nml, "fsd_nonconv_ice_const"], fsd_nonconv)
        self.add_setting(config, [nml, "fsd_nonconv_liq_const"], fsd_nonconv)

        return config, self.reports
