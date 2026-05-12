# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''

'''

import logging
from psyclone.psyir.transformations import (
        ArrayAssignment2LoopsTrans,
        OMPLoopTrans,
        OMPMinimiseSyncTrans,
        TransformationError,
        MaximalOMPParallelRegionTrans
)
from psyclone.psyir.nodes import (
        Assignment,
        Directive,
        Loop,
        Routine,
)
#from transmute_psytrans.transmute_functions import (
from transmute_psytrans.transmute_functions import (
    loop_replacement_of,
    set_pure_subroutines,
    OMP_DO_LOOP_TRANS_STATIC)

Loop.set_loop_type_inference_rules({
    "i_wt": {"variable": "i_wt"},  # For ex_flux_tq.F90
    "ient": {"variable": "ient"},  # For tr_mix.F90
    "ii": {"variable": "ii"},  # For bdy_impl3.F90
    "k": {"variable": "k"},  # For all files which use script
    "j": {"variable": "j"},  # For all files which use script
    "i": {"variable": "i"}})  # For bdy_impl3.F90

# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
# pylint: disable=too-many-branches
def trans(psyir):
    ''' Adds OpenMP Loop directives with nowait to Nemo loops over levels.
    This is followed by applying OpenMP parallel directives as required
    with the OMPMaximalParallelRegionTrans, before removing barriers where
    possible.

    :param psyir: the PSyIR of the provided file.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    #loop_trans = OMPLoopTrans()
    minsync_trans = OMPMinimiseSyncTrans()

    # Remove any j loops and add an init for j
    remove_loop_type = ["j"]

    # Remove any loops relating to specified loop type
    for node in psyir.walk(Routine):
        for removal_type in remove_loop_type:
            loop_replacement_of(node, removal_type)

    safe_pure_calls = ["oneover_v"]

    ignore_dependencies_for =[
            "ct_ctq",
            "dqw",
            "dtl",
            "temp",
            "temp_out",
            "ctctq1",
            "dqw1",
            "dtl1",
            "l"
            ]

    force_private = [
            "temp",
            "temp_out",
            "l"
            ]

    # Set the pure calls if needed
    if safe_pure_calls:
        set_pure_subroutines(psyir, safe_pure_calls)            

    # First convert assignments to loops whenever possible
    for assignment in psyir.walk(Assignment):
        try:
            ArrayAssignment2LoopsTrans().apply(assignment)
        except TransformationError:
            pass

    # Apply OMP_DO_LOOP_TRANS_STATIC to all the loops possible.
    for loop in psyir.walk(Loop):
        if not loop.ancestor(Directive):
            try:
                OMP_DO_LOOP_TRANS_STATIC.apply(loop, ignore_dependencies_for=ignore_dependencies_for, nowait=True)
            except TransformationError as err:
                # Not all of the loops in the example can be parallelised.
                print(f"Add OMP loop: {err}")
                pass

    # Apply the largest possible parallel regions and remove any barriers that
    # can be removed.
    for routine in psyir.walk(Routine):
        #try:
        MaximalOMPParallelRegionTrans().apply(routine, force_private=force_private)
        # except TransformationError as err:
        #     # Not all of the loops in the example can be parallelised.
        #     print(f"Span Parallel: {TransformationError}")
        #     pass
        minsync_trans.apply(routine)
