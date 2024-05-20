# Read Me: OCIM-BGC with Cell Model
This directory contains the model code for running OCIM-BGC using a cellular growth model to compute the organic matter stoichiometry coupling the elemental cycles

## File structure

## Development / Notes

### issues
Error in Hessian for parameter pairs with Q10Photo and RO2C
- Q10Photo gradient looks fine, so error is in second derivative w.r.t. Cparam, Q10Photo
- hessian seems to work for: d, R_Si, and rR pairs with Q10Photo, but none of the other pairs with Q10Photo. Maybe the deriv of these does not depend on stoichiometry? 
- Cxx RHS equations are identical for all cell model parameters. So error is probably not in the RHS eqns.  
- possible cause:
   - mismatch between order of parameter pair columns in Cxx,Oxx and neglogpost
   - if rO2C is entered before cell model params somewhere.
	- this would make the first cell model parameter (Q10Photo) and rO2C wrong. the rest would work
	- Test: if RO2C works for the model with the cell model off, then this is probably the problem 
		- Conclusion: it does.
		- how to fix it?
			- [ ] check order of indices in loops. order of pairs at index kk needs to match across 
				- eqOcycle_v2
				- neglogpost
				- pindx
			- [ ] in eqOcycle, all pairs with O should come at the end after all the pairs computed in eqCcycle
         - [ ] update order of neglogpost loops to match

----
### list of src files that were modified to run the Cell model
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