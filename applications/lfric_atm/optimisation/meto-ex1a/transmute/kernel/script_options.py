# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
This file lifts optional overrides, where possible, into a localised,
single location.
The goal is to allow the creation of a simpler local.py which controls how OMP
is added around loops.
SCRIPT_OPTIONS_DICT is pulled into the local.py script, and overrides some
local functionality.
Each key represents a different kernel file which uses the local script.
Each file key in turn contains a dictionary of overrides.
'''

# Needs to be lifted and likely set by the build system longer term.
# The filename passed to PSyclone, this is the pre-processed FTN source.
FILE_EXTEN = ".xu90"

# Basic initialisation
SCRIPT_OPTIONS_DICT = {}

# File keys
# Disabled due build slow downs CCE, psyclone Issue to fix this is #3418
SCRIPT_OPTIONS_DICT["aerosol_ukca_kernel_mod"+str(FILE_EXTEN)] = {
        "ignore_dependencies_for": [
            "surf_wetness", "l_tile_active", "o3p", "o3",
            "n", "no", "no3", "lumped_n", "n2o5", "ho2no2",
            "hono2", "h2o2", "ch4", "co", "hcho", "meoo",
            "meooh", "h", "oh", "ho2", "cl", "cl2o2", "clo",
            "oclo", "br", "lumped_br", "brcl", "brono2",
            "n2o", "lumped_cl", "hocl", "hbr", "hobr", "hobr",
            "clono2", "cfcl3", "cf2cl2", "mebr", "hono", "c2h6",
            "etoo", "etooh", "mecho", "meco3", "pan", "c3h8",
            "n_proo", "i_proo", "i_proo", "i_prooh", "etcho",
            "etco3", "me2co", "mecoch2oo", "mecoch2ooh", "ppan",
            "meono2", "c5h8", "iso2", "isooh", "ison", "macr",
            "macrooh", "macro2", "mpan", "hacet", "mgly",
            "nald", "hcooh", "meco3h", "meco2h", "h2", "meoh",
            "n_prooh", "msa", "nh3", "cs2", "csul", "h2s", "so3",
            "passive_o3", "age_of_air", "dms", "so2", "h2so4",
            "dmso", "monoterpene", "secondary_organic", "n_nuc_sol",
            "nuc_sol_su", "nuc_sol_om", "n_ait_sol", "ait_sol_su",
            "nuc_sol_om", "n_ait_sol", "ait_sol_su", "ait_sol_om",
            "n_acc_sol", "acc_sol_su", "acc_sol_bc", "acc_sol_om",
            "acc_sol_ss", "n_cor_sol", "cor_sol_su", "cor_sol_bc",
            "cor_sol_om", "cor_sol_ss", "n_ait_ins", "ait_ins_bc",
            "ait_ins_om", "n_acc_ins", "acc_ins_du", "n_cor_ins",
            "cor_ins_du", "cloud_drop_no_conc", "drydp_ait_sol",
            "drydp_acc_sol", "drydp_cor_sol", "drydp_ait_ins",
            "drydp_acc_ins", "drydp_cor_ins", "wetdp_ait_sol",
            "wetdp_acc_sol", "wetdp_cor_sol", "rhopar_ait_sol",
            "rhopar_acc_sol", "rhopar_cor_sol", "rhopar_ait_ins",
            "rhopar_acc_ins", "rhopar_cor_ins", "pvol_su_ait_sol",
            "pvol_bc_ait_sol", "pvol_om_ait_sol", "pvol_wat_ait_sol",
            "pvol_su_acc_sol", "pvol_bc_acc_sol", "pvol_om_acc_sol",
            "pvol_ss_acc_sol", "pvol_wat_acc_sol", "pvol_su_cor_sol",
            "pvol_bc_cor_sol", "pvol_om_cor_sol", "pvol_ss_cor_sol",
            "pvol_wat_cor_sol", "pvol_bc_ait_ins", "pvol_om_ait_ins",
            "pvol_du_acc_ins", "pvol_du_cor_ins", "no2", "bro", "hcl",
            "o1d", "ait_sol_bc",
    ]
}

SCRIPT_OPTIONS_DICT["bl_exp_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "tnuc_nlcl", "dw_bl", "surf_interp", "bq_bl", "bt_bl",
        "dtrdz_tq_bl", "level_ent", "level_ent_dsc", "ent_we_lim",
        "ent_t_frac", "ent_zrzi", "ent_we_lim_dsc", "ent_t_frac_dsc",
        "ent_zrzi_dsc", "fd_taux", "fd_tauy", "gradrinr", "rhokm_bl",
        "rhokh_bl", "moist_flux_bl", "heat_flux_bl", "gradrinr", "lmix_bl",
        "ngstress_bl", "tke_bl", "bl_weight_1dbl", "zh_nonloc", "z_lcl",
        "inv_depth", "qcl_at_inv_top", "shallow_flag", "uw0_flux", "vw0_flux",
        "lcl_height", "parcel_top", "level_parcel_top", "wstar_2d", "thv_flux",
        "parcel_buoyancy", "qsat_at_lcl", "bl_type_ind", "visc_m_blend",
        "visc_h_blend", "zh_2d", "zhsc_2d", "ntml_2d", "cumulus_2d",
        "rh_crit", "mix_len_bm", "dsldzm", "wvar", "zht", "oblen",
        ]
}

SCRIPT_OPTIONS_DICT["bl_imp_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "dqw_wth", "dtl_wth", "dqw_nt_wth", "dtl_nt_wth",
        "ct_ctq_wth", "qw_wth", "tl_wth", "dqw1_2d",
        "dtl1_2d", "ct_ctq1_2d",
        ]
}

SCRIPT_OPTIONS_DICT["bl_imp2_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "fric_heating_blyr", "fric_heating_incv", "nblyr", "cca",
        "ccw", "z_lcl", "cf_area", "cf_bulk", "cf_ice", "cf_liq",
        "dtheta_bl", "m_v", "m_cl", "m_s", "fqw_star_w3", "ftl_star_w3",
        "heat_flux_bl", "moist_flux_bl", "rhokh_bl", "dqw_wth",
        "dtl_wth", "qw_wth", "tl_wth"
        ]
}

# Disabled due build slow downs CCE, psyclone Issue to fix this is #3418
SCRIPT_OPTIONS_DICT["conv_comorph_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "precfrac", "cca_2d", "cape_diluted", "lowest_cca_2d",
        "dt_conv", "dmv_conv", "dmcl_conv", "dms_conv", "m_ci",
        "m_r", "m_g", "massflux_up_half", "massflux_up", "massflux_down",
        "pressure_inc_env", "entrain_up", "entrain_down", "detrain_up",
        "du_conv", "dv_conv", "conv_prog_dtheta", "dmv_conv", "dbcf_conv",
        "dcff_conv", "dcfl_conv", "dt_conv", "dmv_conv", "dmcl_conv",
        "dms_conv", "dd_mf_cb", "cca", "ccw", "cv_top", "cv_base",
        "lowest_cv_top", "lowest_cv_base", "pres_cv_top", "pres_cv_base",
        "pres_lowest_cv_top", "pres_lowest_cv_base", "tke_bl", "o3p", "o1d",
        "o3", "nit", "no", "no3", "lumped_n", "n2o5", "ho2no2", "hono2",
        "h2o2", "ch4", "co", "hcho", "meoo", "meooh", "h", "oh", "ho2",
        "cl", "cl2o2", "clo", "oclo", "br", "lumped_br", "brcl", "brono2",
        "n2o", "lumped_cl", "hocl", "hbr", "hobr", "clono2", "cfcl3", "cf2cl2",
        "mebr", "hono", "c2h6", "etoo", "etooh", "mecho", "meco3", "pan",
        "c3h8", "n_proo", "i_proo", "n_prooh", "i_prooh", "etcho", "etcho",
        "me2co", "mecoch2oo", "mecoch2ooh", "ppan", "meono2", "c5h8", "iso2",
        "isooh", "ison", "macr", "macro2", "macrooh", "mpan", "hacet", "mgly",
        "nald", "hcooh", "hcooh", "meco2h", "h2", "meoh", "msa", "nh3", "cs2",
        "csul", "h2s", "so3", "passive_o3", "age_of_air", "dms", "so2",
        "h2so4", "dmso", "monoterpene", "secondary_organic", "n_nuc_sol",
        "nuc_sol_su", "nuc_sol_om", "n_ait_sol", "ait_sol_su", "ait_sol_bc",
        "ait_sol_om", "n_acc_sol", "acc_sol_su", "acc_sol_bc", "acc_sol_om",
        "acc_sol_ss", "n_cor_sol", "cor_sol_su", "cor_sol_bc", "cor_sol_bc",
        "cor_sol_om", "cor_sol_ss", "n_ait_ins", "ait_ins_bc", "ait_ins_om",
        "n_acc_ins", "acc_ins_du", "n_cor_ins", "cor_ins_du", "detrain_down",
        "conv_prog_dmv", "etco3", "etco3", "meco3h",
        ]
}

SCRIPT_OPTIONS_DICT["mphys_kernel_mod"+str(FILE_EXTEN)] = {
    "node_type_check": False,
    "ignore_dependencies_for": [
        "dml_wth", "dmv_wth", "dml_wth", "dtheta", "refl_1km", "dms_wth",
        "dmr_wth", "dmg_wth", "murk", "dbcf_wth", "dcff_wth", "dcfl_wth",
        "ls_rain_2d", "ls_snow_2d", "ls_graup_2d", "lsca_2d", "ls_rain_3d",
        "ls_snow_3d", "precfrac", "autoconv", "accretion", "rim_cry",
        "rim_agg", "refl_tot", "superc_liq_wth", "superc_rain_wth",
        "sfwater", "sfrain", "sfsnow",
    ]
}
