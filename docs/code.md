# Usage

Code structure of Fleche is described.  
This document has been updated.

<br/>

## Structure

- commandline.py  
- settings.py  
- design.py

<br/>

## Flowchart

1. Preprocessing  
  1.1. Get a list of target molecules (target_molecules.inp) by modified APDFT.  
  1.2. Get a stable geometry by performing a geometry optimization of
       the reference molecule by electronic structure theory to be used in modified APDFT.  

2. Design  
  2.1. Get weighted potential energy (and weighted atomic forces) by performing
       a modified APDFT calculation.  
  2.2. Perturbed participation coefficients.  
  2.3.1. Get a stable geometry and potential energies at the end of a geometry optimization
         using a modified APDFT calculation.  
  2.3.2. Get weight potential energies and derivative of objective function with respect to
         participation coefficients.  
  2.3.3. Update participation coefficients.  
  2.3.4. Back to 2.3.1.  

<br/>

## How to add input parameters?

1. change settings.py
2. change commandline.py
3. change design.py