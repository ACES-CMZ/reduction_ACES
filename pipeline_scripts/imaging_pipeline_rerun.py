from casarecipes.almahelpers import fixsyscaltimes # SACM/JAO - Fixes
from casatasks import fixplanets
context = h_init()
context.set_state('ProjectSummary', 'proposal_code', '2021.1.00172.L')
context.set_state('ProjectSummary', 'proposal_title', 'unknown')
context.set_state('ProjectSummary', 'piname', 'unknown')

import os
import sys
cwd = os.getcwd()
mous = cwd.split("/")[-3].split(".")[1].split("_")
ous_entity_id = f"{mous[0]}://{mous[3]}/{mous[4]}/{mous[5]}"

context.set_state('ProjectStructure', 'ous_entity_id', ous_entity_id)
context.set_state('ProjectStructure', 'ous_part_id', 'Undefined')
context.set_state('ProjectStructure', 'ous_title', 'Undefined')
context.set_state('ProjectStructure', 'ousstatus_entity_id', ous_entity_id)
context.set_state('ProjectStructure', 'recipe_name', 'hifa_calimage_renorm')

import glob
mses = glob.glob("*.ms")
if len(mses) == 0:
    raise ValueError("No MSes found")

try:
    hifa_importdata(vis=mses, dbservice=False)
    hif_mstransform(pipelinemode="automatic")
    #hifa_flagtargets(pipelinemode="automatic")
    #hifa_imageprecheck(pipelinemode="automatic")    
    #hif_makeimlist(specmode='mfs')
    #hif_findcont(pipelinemode="automatic")
    #hif_uvcontfit(pipelinemode="automatic", fitorder=0)
    #hif_uvcontsub(pipelinemode="automatic")
    #hif_makeimages(pipelinemode="automatic")
    #hif_makeimlist(specmode='cont')
    #hif_makeimages(pipelinemode="automatic")
    hif_makeimlist(specmode='cube')
    hif_makeimages(pipelinemode="automatic")
finally:
    h_save()
