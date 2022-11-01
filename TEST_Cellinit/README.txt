TEST_Cellinit

testing initial values of cell model parameters to confirm that optimization is finding a minimum.
Grid91
BGC parameters set to output of Grid91 run with GM15 C:P. 
loadfile = “/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK91/testNPP_CTL_He_PCfixb_DOC0.25_DOP0_xhat.mat”

output saves to folder:
/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK91/testCellinit/

filename notation
a = alpha
q= Q10Photo
f = fStorage
k = kST0
r = PStor_rCutoff
g = gammaDNA
s = alphaS

variation
0 = not optimized. set at default
1 = initial value option 1
2 = initial value option 2 

