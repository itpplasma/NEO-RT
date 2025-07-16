#!/usr/bin/env python3
"""
Convert tok_circ format (inp_swi = 8) to ASDEX format (inp_swi = 9)
Converts 4-column magnetic field data to 8-column ASDEX format
"""

import re

def convert_tok_circ_to_asdex(input_file, output_file):
    """
    Convert tok_circ format to ASDEX format
    
    tok_circ format (4 columns after m,n):
    m    n    r[m]    z[m]    (phib-phi)*nper/twopi    bmn[T]
    
    ASDEX format (8 columns after m,n):
    m    n    rmnc[m]    rmns[m]    zmnc[m]    zmns[m]    bmnc[T]    bmns[T]
    """
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    with open(output_file, 'w') as f:
        for line in lines:
            # Skip comment lines and empty lines
            if line.strip().startswith('CC') or line.strip() == '':
                f.write(line)
                continue
            
            # Check if this is a header line
            if 'r/[m]' in line and 'z/[m]' in line and 'bmn/[T]' in line:
                # Convert header for ASDEX format
                f.write("    m    n      rmnc/[m]        rmns/[m]        zmnc/[m]        zmns/[m]     bmnc/[T]     bmns/[T]\n")
                continue
            
            # Check if this is a data line (starts with whitespace followed by numbers)
            if re.match(r'\s*\d+\s+\d+\s+', line):
                # Parse the data line
                parts = line.split()
                if len(parts) >= 6:  # m, n, r, z, phi, bmn
                    m, n, r, z, phi, bmn = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
                    
                    # For axisymmetric tokamak, sine components are zero
                    # rmnc = r, rmns = 0, zmnc = z, zmns = 0, bmnc = bmn, bmns = 0
                    asdex_line = f"{m:>5} {n:>4} {r:>15} {float(0.0):>15.8e} {z:>15} {float(0.0):>15.8e} {bmn:>12} {float(0.0):>12.8e}\n"
                    f.write(asdex_line)
                else:
                    # Keep line as is if it doesn't match expected format
                    f.write(line)
            else:
                # Keep line as is (headers, comments, etc.)
                f.write(line)

if __name__ == "__main__":
    # Convert the current in_file to ASDEX format
    convert_tok_circ_to_asdex('/home/ert/code/NEO-RT/in_file', '/home/ert/code/NEO-RT/in_file_asdex')
    print("Converted tok_circ format to ASDEX format")
    print("Input:  /home/ert/code/NEO-RT/in_file (tok_circ format)")
    print("Output: /home/ert/code/NEO-RT/in_file_asdex (ASDEX format)")