# sequentialScanner
**Description**

- A customisable sequential gene scanner for optimising combinations of both gene deletions and upregulations using flux balance analysis
- Allows adapting the number of affected genes, the loop selection percentile, the upregulation factor & the growth rate threshold
- Created by: [@LucoDevro](https://github.com/LucoDevro), Lucas De Vrieze, Department of Chemical Engineering, KU Leuven, Heverlee (Belgium)
- Developed as part of a MSc thesis under supervision of prof. dr. ir. Kristel Bernaerts and prof. dr. ir. Steffen Waldherr
- Last updated: 14th April 2021

**Theoretical basis**

Sequential scanning algorithm based on pFBA, including MOMA extension to deal with genetically perturbed networks
- *Scanning algorithm:* Alper et al., Met. Eng., doi: https://doi.org/10.1016/j.ymben.2004.12.003
- *MOMA integration:* Wang et al., Biochem Eng J, doi: https://doi.org/10.1016/j.bej.2017.03.017

**Contents**

- The master function (`sequentialScanner.m`)
- The worker function (`sequentialScannerWorker.m`)
- Auxiliary function altering the model flux bounds for a gene upregulation (`upregulateModelGenes.m`)
- Auxiliary functions interconverting consensus gene names and systematic gene names (`GetGeneNameFromSysGeneName.m` & `GetSysGeneNameFromGeneName`)
- `MOMA.m` from the COBRA Toolbox, with added bypass for wild-type solutions
- An example overhead script calling the main function with appropiate parameters (`sequentialScannerOverhead.m`): optimise ammonium efflux
- The COBRA model necessary for this example script (`model100.mat`): *Bacillus subtilis* in protein hydrolysate medium (genome-scale metabolic model iYO844 with custom-made GECKO extension)

**Installation & software dependencies**
---
**Required software**
- MATLAB 2019b or higher
- [COBRA toolbox](https://github.com/opencobra/cobratoolbox) 3.0 for MATLAB
- IBM Cplex optimiser for MATLAB (or Gurobi Optimiser, or any other software capable of solving quadratic optimisation problems)
