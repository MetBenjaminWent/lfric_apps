# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
A local.py script for all kernels, where instead of adding OMP across the
outermost loop, it is placed around the i loop, or across the l loop.
This script imports a SCRIPT_OPTIONS_DICT which can be used to override
small aspects of this script per file it is applied to.
Overrides currently include:
* ignore_dependencies_for
* node_type_check
* safe_pure_calls
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
    get_children,
    get_all_children,
    get_outer_loops,
    OMP_PARALLEL_REGION_TRANS,
    OMP_DO_LOOP_TRANS_STATIC,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
)

ignore_dependencies_for = [
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

def trans(psyir):
    '''
    PSyclone function call, run through psyir object,
    each schedule (or subroutine) and apply paralleldo transformations
    to each loop.
    '''

    # declare 'case_default_used'. To get the correct subroutine,
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

    # Work through each loop in the file and OMP PARALLEL DO
    #for loop in psyir.walk(Loop):

    # Identify outer loops
    outer_loops = [loop for loop in get_outer_loops(psyir)
                if not loop.ancestor(Loop)]

    for index, loop in enumerate(outer_loops):
        #if loop.variable.name == 'm' and index in [0]:
        if loop.variable.name == 'm':
            # Are there any loops?
            if get_all_children(loop, node_type=Loop):
                move_default_case_contents(loop)
                print(index)
                nodes_potential = get_children(loop)
                print(nodes_potential)
                try:
                    OMP_PARALLEL_REGION_TRANS.apply(nodes_potential[1:-1])
                except (TransformationError, IndexError) as err:
                    logging.warning(
                        "%s: Could not transform because:\n %s", index, err)

    print(outer_loops)
    print(len(outer_loops))

    # # Loop m 1:
    # nodes_potential = get_children(outer_loops[0])
    # try:
    #     OMP_PARALLEL_REGION_TRANS.apply(nodes_potential[1:-1])
    # except (TransformationError, IndexError) as err:
    #     logging.warning(
    #         "Could not transform because:\n %s", err)

    # # Loop m 2:
    # nodes_potential = get_children(outer_loops[1])
    # try:
    #     OMP_PARALLEL_REGION_TRANS.apply(nodes_potential[1:-1])
    # except (TransformationError, IndexError) as err:
    #     logging.warning(
    #         "Could not transform because:\n %s", err)

    # Loops 3,4:
    try:
        OMP_PARALLEL_REGION_TRANS.apply(outer_loops[2:4])
    except (TransformationError, IndexError) as err:
        logging.warning(
            "Could not transform because:\n %s", err)


    # Work through all loops
    for loop in psyir.walk(Loop):
        # For each loop which is inside a OMPParallelDirective, and not a OMPDoDirective
        # parallelise
        if loop.ancestor(OMPParallelDirective) and not loop.ancestor(OMPDoDirective):
            try:
                OMP_DO_LOOP_TRANS_STATIC.apply(
                    loop)
            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)


        # if loop.variable.name == 'm':
        #     # print(loop.variable.name)
        #     # children = get_all_children(loop, node_type=Loop)
        #     generic_children = get_children(loop)
        #     loop_grandchildren = get_children(generic_children[0], node_type=Loop)
        #     if loop_grandchildren:
        #         print(loop_grandchildren)
        #         if isinstance(loop_grandchildren[0], Loop):
        #             print("Success")


        # # If there is an OMP ancestor skip.
        # if (
        #     loop.ancestor(OMPParallelDoDirective) is not None
        #     or loop.ancestor(OMPDoDirective) is not None
        #     or loop.ancestor(OMPParallelDirective) is not None
        # ):
        #     continue
        # # Allow loops over 'i' and 'l' indexes to be parallelised.
        # if loop.variable.name in ['i', 'l']:
        #     try:
        #         OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(
        #             loop, 
        #             ignore_dependencies_for=ignore_dependencies_for,
        #             node_type_check=node_type_check)
        #     except (TransformationError, IndexError) as err:
        #         logging.warning(
        #             "Could not transform because:\n %s", err)

def move_default_case_contents(loop):
    """
    Seems to be safe atm
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
