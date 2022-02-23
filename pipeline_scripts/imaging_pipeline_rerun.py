from casarecipes.almahelpers import fixsyscaltimes # SACM/JAO - Fixes
from casatasks import fixplanets
h_init()
try:
    hif_makeimlist(specmode='mfs')
    hif_findcont(pipelinemode="automatic")
    hif_uvcontfit(pipelinemode="automatic", fitorder=0)
    hif_uvcontsub(pipelinemode="automatic")
    hif_makeimages(pipelinemode="automatic")
    hif_makeimlist(specmode='cont')
    hif_makeimages(pipelinemode="automatic")
    hif_makeimlist(specmode='cube')
    hif_makeimages(pipelinemode="automatic")
finally:
    h_save()
