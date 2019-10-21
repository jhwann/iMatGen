#After generating materials do following heuristic post-processing procedure. 

1. Do inverse transform for basis and cell images using 'img2basis.py' and 'img2cell.py'.

2. Check that all the cells are well transformed, and remove invalid unit-cell structures ('nan' in cif file or too largely deformed shape).

3. Combine cell and basis using 'combine_basis_cell.py'.

4. Check overlapped position due to periodic boundary conditions using 'combined_pbc_chk.py'.
