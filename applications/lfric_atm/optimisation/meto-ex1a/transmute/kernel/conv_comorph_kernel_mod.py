# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
A copy of the override for aerosol_ukca_kernel_mod, used with
conv_comorph_kernel_mod. Adding OMP parallel do to all
loops causes the CCE compiler to fail to build the file within a reasonable
timeframe, or at all.
We can still add them to standard i loops, but for the loops for the tracers,
it requires something a little more bespoke.
'''

import logging
from psyclone.psyir.transformations import (
    ArrayAssignment2LoopsTrans,
    OMPLoopTrans,
    OMPMinimiseSyncTrans,
    TransformationError,
    MaximalOMPParallelRegionTrans,
)
from psyclone.psyir.nodes import (
    Loop, Call,
    Routine,
    Schedule,
    Literal,
    Reference,
    IfBlock,
    Fparser2CodeBlock,
    Assignment,
    OMPParallelDoDirective,
    OMPParallelDirective,
    OMPDoDirective,)
from psyclone.psyir.symbols import DataSymbol, ScalarType
from transmute_psytrans.transmute_functions import (
    get_ancestors,
    get_children,
    get_all_children,
    get_outer_loops,
    OMP_PARALLEL_REGION_TRANS,
    OMP_DO_LOOP_TRANS_STATIC,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
)

ignore_dependencies_for = [
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

def trans(psyir):
    '''
    PSyclone function call, run through psyir object.
    * Insert manual parallel regions around specific nodes.
    * Re-organise some if blocks to allow spanning parallel regions.
    * Add OMP do inside these regions.
    * Insert parallel do around remaining loop nodes.

    :param psyir: the PSyIR of the provided file.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    '''

    ### Over the tracers, where safe with current PSyclone ###

    # Declare 'case_default_used'. To get the correct subroutine,
    # Just jump to the first loop, grab it's schedule and check
    # the symbol table, which itself points at the routine.
    for loop in psyir.walk(Loop):
        for schedule in loop.walk(Schedule):
            symtab = schedule.symbol_table
            case_default_used = symtab.find_or_create(
            "case_default_used",
                symbol_type=DataSymbol,
                datatype=ScalarType.boolean_type())
            break
        break

    # # Identify outer loops
    # outer_loops = [loop for loop in get_outer_loops(psyir)
    #             if not loop.ancestor(Loop)]

    # # Work through the outer loops, those that are over n,
    # # readjust the tracer loop to move the logging
    # # Then span a parallel section over the ifblock
    # for index, loop in enumerate(outer_loops):
    #     if (get_all_children(loop, node_type=Loop) and 
    #         loop.variable.name == 'n'):
    #         # This node needs the logging moved
    #         if index in [8]:
    #             move_default_case_contents(loop)
    #         nodes_potential = get_children(loop)
    #         # As this node didn't get moved, use different range
    #         if index in [56]:
    #             nodes_span = [0, 1]
    #         else:
    #             nodes_span = [1,-1]
    #         try:
    #             OMP_PARALLEL_REGION_TRANS.apply(nodes_potential[
    #                 nodes_span[0]:nodes_span[1]],
    #                 force_private=["case_default_used"])
    #         except (TransformationError, IndexError) as err:
    #             logging.warning(
    #                 f"{index}: Could not transform because:\n {err}")

    ### End over tracers ###

    # To reduce some parallel sections for CCE, we can group these nodes
    top_nodes = []
    for index, node in enumerate(psyir.walk(Routine)[0].children):
        top_nodes.append(node)
    try:
        OMP_PARALLEL_REGION_TRANS.apply(top_nodes[0:2])
    except (TransformationError, IndexError) as err:
        logging.warning(
            f"Could not transform because:\n {err}")
    # try:
    #     OMP_PARALLEL_REGION_TRANS.apply(top_nodes[7:13])
    # except (TransformationError, IndexError) as err:
    #     logging.warning(
    #         f"Could not transform because:\n {err}")
    # try:
    #     OMP_PARALLEL_REGION_TRANS.apply(top_nodes[17:20])
    # except (TransformationError, IndexError) as err:
    #     logging.warning(
    #         f"Could not transform because:\n {err}")
    # try:
    #     OMP_PARALLEL_REGION_TRANS.apply(top_nodes[23:25])
    # except (TransformationError, IndexError) as err:
    #     logging.warning(
    #         f"Could not transform because:\n {err}")
    # try:
    #     OMP_PARALLEL_REGION_TRANS.apply(top_nodes[51:56])
    # except (TransformationError, IndexError) as err:
    #     logging.warning(
    #         f"Could not transform because:\n {err}")
    # try:
    #     OMP_PARALLEL_REGION_TRANS.apply(top_nodes[69:75])
    # except (TransformationError, IndexError) as err:
    #     logging.warning(
    #         f"Could not transform because:\n {err}")
    # try:
    #     OMP_PARALLEL_REGION_TRANS.apply(top_nodes[78:80])
    # except (TransformationError, IndexError) as err:
    #     logging.warning(
    #         f"Could not transform because:\n {err}")
    # try:
    #     OMP_PARALLEL_REGION_TRANS.apply(top_nodes[81:89])
    # except (TransformationError, IndexError) as err:
    #     logging.warning(
    #         f"Could not transform because:\n {err}")

    # To add do inside any spanned parallel sections
    for loop in psyir.walk(Loop):
        # For each loop which is inside a OMPParallelDirective, and not a OMPDoDirective
        # parallelise
        if loop.ancestor(OMPParallelDirective) and not loop.ancestor(OMPDoDirective):
            if loop.variable.name in ['i', 'j']:
                # To add do inside any spanned parallel sections
                try:
                    OMP_DO_LOOP_TRANS_STATIC.apply(
                        loop,
                        ignore_dependencies_for=ignore_dependencies_for)
                except (TransformationError, IndexError) as err:
                    logging.warning(
                        f"{index}: Could not transform because:\n {err}")
                    print(err)

    # # To add parallel do (sparingly) around any remaining 'i' or 'k' loops
    # for loop in psyir.walk(Loop):
    #     if (
    #         loop.ancestor(OMPParallelDoDirective) is not None
    #         or loop.ancestor(OMPDoDirective) is not None
    #         or loop.ancestor(OMPParallelDirective) is not None
    #     ):
    #         continue

    #     # We don't want to add parallel do's around loops inside 'n' loops
    #     # Where possible, this has already been done.
    #     loop_ancestors = get_ancestors(loop, node_type=Loop)
    #     found_outer_ancestor = False
    #     for loop_ancestor in loop_ancestors:
    #         if str(loop_ancestor.variable.name) in ['n']:
    #             found_outer_ancestor = True
    #     if found_outer_ancestor:
    #         continue

    #     # To add parallel do (sparingly) around any remaining 'i' loops
    #     if loop.variable.name == 'i':
    #         try:
    #             OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(
    #                 loop,
    #                 ignore_dependencies_for=ignore_dependencies_for)
    #         except (TransformationError, IndexError) as err:
    #             logging.warning(
    #                 f"{index}: Could not transform because:\n {err}")
    #             print(err)


def move_default_case_contents(loop):
    """
    This kernel requires the moving of a logging call out of an generated
    ifblock (which originally an case) for the else (default) clause, so that
    a parallel section can be spanned over the ifblock. This reduces the number
    parallel sections present in the file.
    So that the moved logging statement is still called as needed, if the
    else clause is used (by a tracer not being present), it flips a boolean,
    which is a flag to control whether the following if clause is called,
    which contains the moved logging statement.
    :arg loop: the Loop node for which we will move the logging statement.
    :type loop: :py:class:`Loop`
    """
    issue_nodes = get_all_children(loop, node_type=Fparser2CodeBlock)
    for issue_node in issue_nodes:
        lhs_false = Reference(DataSymbol(
                "case_default_used",
                ScalarType.boolean_type()))
        lhs_true = Reference(DataSymbol(
                "case_default_used",
                ScalarType.boolean_type()))
        condition = Reference(DataSymbol(
                "case_default_used",
                ScalarType.boolean_type()))
        rhs_false = Literal("false", ScalarType.boolean_type())
        assign_false = Assignment.create(lhs_false, rhs_false)
        rhs_true = Literal("true", ScalarType.boolean_type())
        assign_true = Assignment.create(lhs_true, rhs_true)

        ## case_default_used will alow us to span a parallel section ##
        ## Add it at the start, outside the region. It does not matter ##
        ## which thread sets the value ##
        # loop is the parents parent
            # add the assignment at 0
        loop.loop_body.addchild(assign_false, 0)

        ## Move the case default (else) contents to a new ifblock, ##
        ## controlled by case_default_used ##
        # loop is the parents parent
            # need to add an if block around the new assignment
                # and the Fparser2CodeBlock node here
        node_parent_children = issue_node.parent.children
        # as we detach nodes, node_parent_children shrinks
        # so position 1 is now actually 0
        # if_body = [node.detach() for node in node_parent_children]
        # didn't work, it only got node 1
        if_body = []
        for index2 in range(len(node_parent_children)):
            if_body.append(node_parent_children[0].detach())
        ifblock = IfBlock.create(condition, if_body)
        loop.loop_body.addchild(ifblock, 2)
        ## ##

        ## Find the else node and add ## 
        # loop is the parents parent case_default_used is true
            # actually needs to be added to else body of if node
        if_nodes = get_children(loop, node_type=IfBlock)
        cursor = if_nodes[0].else_body
        # cursor is the schedule of the node, cursor[0] is the node
        while len(cursor.children) == 1 and isinstance(cursor[0], IfBlock):
            if cursor[0].else_body is None:
                cursor = None # The nested if don't end with an else block
                break
            cursor = cursor[0].else_body
        # One it breaks at the case default (else),
        # we have the right part of the if block
        cursor.addchild(assign_true, 0)
        ## ##
