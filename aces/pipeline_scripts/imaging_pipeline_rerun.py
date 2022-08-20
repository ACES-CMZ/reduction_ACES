import glob
import shutil
import os
from casarecipes.almahelpers import fixsyscaltimes  # SACM/JAO - Fixes
context = h_init()
context.set_state('ProjectSummary', 'proposal_code', '2021.1.00172.L')
context.set_state('ProjectSummary', 'proposal_title', 'unknown')
context.set_state('ProjectSummary', 'piname', 'unknown')

cwd = os.getcwd()
mous = cwd.split("/")[-3].split(".")[1].split("_")
ous_entity_id = f"{mous[0]}://{mous[3]}/{mous[4]}/{mous[5]}"

context.set_state('ProjectStructure', 'ous_entity_id', ous_entity_id)
context.set_state('ProjectStructure', 'ous_part_id', 'Undefined')
context.set_state('ProjectStructure', 'ous_title', 'Undefined')
context.set_state('ProjectStructure', 'ousstatus_entity_id', ous_entity_id)
context.set_state('ProjectStructure', 'recipe_name', 'hifa_calimage_renorm')

# don't re-run this on things that already inclue 'target'
# (mstransform will try again and again...)
mses = [x for x in glob.glob("*.ms") if 'target' not in x]
if len(mses) == 0:
    raise ValueError("No MSes found")

print(f"Starting imaging pipeline rerun for mses {mses} and for mouse {mous} in {cwd}")

try:
    hifa_importdata(vis=mses, dbservice=False)

    # if _target.ms exists, we don't need to transform
    if not all(os.path.exists(x.replace(".ms", "_target.ms")) for x in mses):
        hif_mstransform(pipelinemode="automatic")

    # per Tafoya, flagtargets is probably not needed (and it was already done once)
    # hifa_flagtargets(pipelinemode="automatic")

    hifa_imageprecheck(pipelinemode="automatic")

    if os.path.exists("../../calibration/cont.dat") and not os.path.exists('cont.dat'):
        shutil.copy("../../calibration/cont.dat", '.')

    hif_findcont(pipelinemode="automatic")
    hif_uvcontfit(pipelinemode="automatic", fitorder=0)
    hif_uvcontsub(pipelinemode="automatic")

    hif_makeimlist(specmode='mfs')
    hif_makeimages(pipelinemode="automatic")

    hif_makeimlist(specmode='cont')
    hif_makeimages(pipelinemode="automatic")

    hif_makeimlist(specmode='cube')
    hif_makeimages(pipelinemode="automatic")

    hifa_exportdata(pipelinemode="automatic")
finally:
    h_save()
