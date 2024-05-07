# Read Me: OCIM-BGC with Cell Model
This directory contains the model code for running OCIM-BGC using a cellular growth model to compute the organic matter stoichiometry coupling the elemental cycles

## File structure

## Development / Notes

### list of src files that need to be modified to run the Cell model
- `eqCcycle_v2.m`: Update Gradient and HEssian to include Cell model parameters
- `neglogpost.m`: 
   - add code to run the cell model. 
   - define C2P and O2C outside of eqCcycle. 
   - update gradient hessian to check is par.cellmodel is on and use the corresponding parameter count
   - save additional output. 
- `driver.m`: 
   - new filenames
   - new fields in par for cellmodel options and parameters
   - fix option to run without optimizing
- `eqOcycle_v2.m` : add option for respiration quotient (O2C)
- `SetPar.m` : add cell model parameters
- `PackPar.m` : add cell model parameters
- `PrintPara.m` : add cell model parameters
- `ResetPar.m` : add limits for cell model parameters

Optional changes:
- `Fsea2air.m` : add switch for future_case_options, to use projected environmental fields instead of climatological fields.

### New files unique to the cell model
- `CellCHNOP.m`: computes the cellular quotas of carbon, hydrogen, nitrogen, oxygen, and phosphorus for a phytoplankton cell with the optimal growth strategy in the given environment
- `eqC2Puptake.m`

## Ideas for updates
1.  save a config file with the data and options used to run the model. maybe something like... 
   - model_config_Full: saves everything in par
   - model_config_Small: saves just the paths to the input files and all the on/off switches in par

2. Change alphaPhoto to be a field saved in the par structure instead of a function in CellCHNOP.m. compute alphaPhoto from the light field in SetUp.m, by calling CellPhoto.m. 