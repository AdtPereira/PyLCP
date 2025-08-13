# ****** PROGRAM RIBBON.PY ******
# To Determine the NxN Generalized Capacitance Matrix
# of an N-Wire Ribbon Cable. The (N-1)x(N-1) Transmission-Line
# Capacitance and Inductance Matrices are also Computed from these
# Results. The Number of Fourier Coefficients for the Charge
# Distribution Around the Wire Peripheries is Denoted by NF.
#
# Written by: Clayton R. Paul
#
# Translated to Python by Gemini.
# DEFINITIVE VERSION + VALIDATION: This version adds an automatic
# validation step to compare its output (PUL.TXT) against the
# reference Fortran output (PUL.DAT).


import os
import re
import numpy as np
from pathlib import Path

def parse_pul_file(filename):
    """Helper function to parse a PUL file and extract L, C, CGEN and C0 values."""
    results = {'L': {}, 'C': {}, 'C0': {}, 'CGEN': {}}
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # CORREÇÃO DEFINITIVA APLICADA AQUI:
                # Trocado re.match() por re.search() para encontrar o padrão
                # em qualquer lugar da linha, ignorando espaços no início.
                match = re.search(r"(\d+)\s+(\d+)\s+([0-9.E+-]+)\s+=\s*([A-Z0-9]+)\(", line)
                if match:
                    i, j, value, key = match.groups()
                    i, j, value = int(i), int(j), float(value)
                    
                    # Ensure keys are sorted for consistent comparison, e.g. (1,2) not (2,1)
                    if i > j:
                        i, j = j, i

                    if key in results:
                        results[key][(i, j)] = value
    except FileNotFoundError:
        print(f"\nWarning: Could not find file {filename} for parsing.")
        return None
    return results

def validate_results(python_output_file, fortran_output_file, tolerance=1e-5):
    """
    Compares the results from the Python script's output file with the reference Fortran script's output file.
    """
    print("\n--- STARTING RESULTS VALIDATION ---")
    
    python_results = parse_pul_file(python_output_file)
    fortran_results = parse_pul_file(fortran_output_file)

    if python_results is None or fortran_results is None:
        print("Validation could not be completed: one or both output files were not found.")
        return

    discrepancies = []
    
    all_keys = sorted(list(set(python_results.keys()) | set(fortran_results.keys())))

    for key in all_keys:
        py_matrix = python_results.get(key, {})
        fo_matrix = fortran_results.get(key, {})
        all_indices = sorted(list(set(py_matrix.keys()) | set(fo_matrix.keys())))

        for index in all_indices:
            py_val = py_matrix.get(index)
            fo_val = fo_matrix.get(index)

            if py_val is None or fo_val is None:
                discrepancies.append(
                    f"Index {index} for matrix {key} exists in one file but not the other."
                )
                continue

            if not np.isclose(py_val, fo_val, rtol=tolerance):
                diff = 100 * abs(py_val - fo_val) / abs(fo_val) if fo_val != 0 else float('inf')
                discrepancies.append(
                    f"Discrepancy in Matrix {key}, Index {index}:\n"
                    f"  > Python (PUL.TXT): {py_val:.6E}\n"
                    f"  > Fortran (PUL.DAT): {fo_val:.6E}\n"
                    f"  > Relative Difference: {diff:.4f}%"
                )

    if not discrepancies:
        print("\n-------------------------------------------------")
        print("  VALIDATION PASSED  ")
        print(f"  Results in '{python_output_file}' match '{fortran_output_file}' (tolerance={tolerance}).")
        print("-------------------------------------------------")
    else:
        print("\n-------------------------------------------------")
        print("  VALIDATION FAILED  ")
        print("  The following discrepancies were found:")
        print("-------------------------------------------------")
        for d in discrepancies:
            print(f"- {d}\n")
    
    print("\n--- VALIDATION COMPLETE ---")

def ribbon_dot_for():
    """
    Main function to analyze a ribbon cable and calculate its
    per-unit-length capacitance and inductance matrices.
    """
    
    # Determine the path to the input file
    script_dir = Path(__file__).resolve().parent
    input_file_path = script_dir / 'RIBBON.IN'
    output_filename = script_dir / 'PUL.TXT'
    reference_output_file = script_dir / 'PUL.DAT'
    
    try:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            N = int(lines[0].split()[0])
            NF = int(lines[1].split()[0])
            IREF = int(lines[2].split()[0])
            RW = float(lines[3].split()[0])
            TD = float(lines[4].split()[0])
            ER = float(lines[5].split()[0])
            S = float(lines[6].split()[0])
    except FileNotFoundError:
        print("Error: Input file 'RIBBON.IN' not found.")
        return
    except (ValueError, IndexError):
        print("Error: Invalid 'RIBBON.IN' input file format.")
        return

    # Use float64 for all calculations for precision, similar to standard scientific computing.
    dtype = np.float64

    # --- Constants ---
    PI = np.pi
    V0 = 2.997925E8
    V02 = V0 * V0
    EPS = 1.0 / (V02 * 4.E-7 * PI)

    if N <= 1 or IREF > N or (2. * (RW + TD)) > S:
        with open('PUL.TXT', 'w') as f:
            f.write('INPUT DATA ERROR\n')
        return

    NNF = N * NF
    A = np.zeros((NNF, NNF), dtype=dtype)
    B = np.zeros((NNF, NNF), dtype=dtype)
    C = np.zeros((NNF, NNF), dtype=dtype)
    D = np.zeros((NNF, NNF), dtype=dtype)

    RD = RW + TD
    ANGLE = 2.0 * PI / NF
    DELTA = ANGLE / 4.0

    # --- Step 1: Calculate the self-interaction block (Block 1,1) ---
    A1 = np.zeros((NF, NF), dtype=dtype)
    B1 = np.zeros((NF, NF), dtype=dtype)
    C1 = np.zeros((NF, NF), dtype=dtype)
    D1 = np.zeros((NF, NF), dtype=dtype)
    for i in range(NF):
        ANG = i * ANGLE + DELTA
        for j in range(NF):
            j1 = j + 1
            if j == 0:
                A1[i, j] = -RW * np.log(RW) / EPS
                B1[i, j] = -RD * np.log(RD) / EPS
                C1[i, j] = (ER - 1.0) * RW / RD
                D1[i, j] = -1.0
            else:
                xm1 = dtype(j)
                cos_val = np.cos(xm1 * ANG)
                A1[i, j] = (RW**j1) * cos_val / (2. * EPS * xm1 * RW**j)
                C1[i, j] = (ER - 1.) * ((RW / RD)**j1) * cos_val / 2.0
                D1[i, j] = -(ER + 1.0) * cos_val / 2.0
                if j == 1:
                    B1[i, j] = (RW**j) * cos_val / (2. * EPS * xm1 * 1.0) # RD**(j-1) = RD**0 = 1
                else:
                    B1[i, j] = (RW**j) * cos_val / (2. * EPS * xm1 * (RD**(j - 1)))

    # --- Step 2: Calculate the first row/column of blocks (interaction with wire 1) ---
    for i_dist in range(1, N):
        sep = dtype(i_dist) * S

        # Calculate block (0, i_dist) and (i_dist, 0)
        A_block = np.zeros((NF, NF), dtype=dtype)
        B_block = np.zeros((NF, NF), dtype=dtype)
        C_block = np.zeros((NF, NF), dtype=dtype)
        D_block = np.zeros((NF, NF), dtype=dtype)
        
        A_block_sym = np.zeros((NF, NF), dtype=dtype)
        B_block_sym = np.zeros((NF, NF), dtype=dtype)
        C_block_sym = np.zeros((NF, NF), dtype=dtype)
        D_block_sym = np.zeros((NF, NF), dtype=dtype)

        for j_row in range(NF):
            ang = j_row * ANGLE + DELTA
            # Geometry for A, B blocks
            h_rw = sep - RW * np.cos(ang); v_rw = RW * np.sin(ang)
            th_rw = np.pi - np.arctan2(v_rw, h_rw); rp_rw = np.sqrt(h_rw**2 + v_rw**2)
            h_rw_s = sep + RW*np.cos(ang); v_rw_s = RW*np.sin(ang)
            th_rw_s = np.arctan2(v_rw_s, h_rw_s); rp_rw_s = np.sqrt(h_rw_s**2+v_rw_s**2)

            # Geometry for C, D blocks
            h_rd = sep - RD*np.cos(ang); v_rd = RD*np.sin(ang)
            th_rd = np.pi - np.arctan2(v_rd, h_rd); rp_rd = np.sqrt(h_rd**2+v_rd**2)
            rdn = np.cos(th_rd - ang); tdn = -np.sin(th_rd-ang)
            h_rd_s = sep+RD*np.cos(ang); v_rd_s = RD*np.sin(ang)
            th_rd_s = np.arctan2(v_rd_s,h_rd_s); rp_rd_s = np.sqrt(h_rd_s**2+v_rd_s**2)
            rdn_s = np.cos(th_rd_s - ang); tdn_s = -np.sin(th_rd_s - ang)

            for k_col in range(NF):
                k_fortran = k_col + 1
                if k_col == 0:
                    A_block[j_row,k_col] = -RW*np.log(rp_rw)/EPS
                    B_block[j_row,k_col] = -RD*np.log(rp_rw)/EPS
                    A_block_sym[j_row, k_col] = -RW*np.log(rp_rw_s)/EPS
                    B_block_sym[j_row, k_col] = -RD*np.log(rp_rw_s)/EPS
                    C_block[j_row,k_col] = (ER-1.)*(RW/rp_rd)*rdn
                    D_block[j_row,k_col] = (ER-1.)*(RD/rp_rd)*rdn
                    C_block_sym[j_row,k_col] = (ER-1.)*(RW/rp_rd_s)*rdn_s
                    D_block_sym[j_row,k_col] = (ER-1.)*(RD/rp_rd_s)*rdn_s
                else:
                    xkm1 = dtype(k_col)
                    A_block[j_row,k_col] = (RW**k_fortran)*np.cos(xkm1*th_rw)/(2.*EPS*xkm1*(rp_rw**(k_fortran-1)))
                    B_block[j_row,k_col] = (RD**k_fortran)*np.cos(xkm1*th_rw)/(2.*EPS*xkm1*(rp_rw**(k_fortran-1)))
                    A_block_sym[j_row,k_col] = (RW**k_fortran)*np.cos(xkm1*th_rw_s)/(2.*EPS*xkm1*(rp_rw_s**(k_fortran-1)))
                    B_block_sym[j_row,k_col] = (RD**k_fortran)*np.cos(xkm1*th_rw_s)/(2.*EPS*xkm1*(rp_rw_s**(k_fortran-1)))
                    
                    C_block[j_row,k_col] = (ER-1.)*(((RW/rp_rd)**k_fortran)/2.)*(np.cos(xkm1*th_rd)*rdn + np.sin(xkm1*th_rd)*tdn)
                    D_block[j_row,k_col] = (ER-1.)*(((RD/rp_rd)**k_fortran)/2.)*(np.cos(xkm1*th_rd)*rdn + np.sin(xkm1*th_rd)*tdn)
                    C_block_sym[j_row,k_col] = (ER-1.)*(((RW/rp_rd_s)**k_fortran)/2.)*(np.cos(xkm1*th_rd_s)*rdn_s + np.sin(xkm1*th_rd_s)*tdn_s)
                    D_block_sym[j_row,k_col] = (ER-1.)*(((RD/rp_rd_s)**k_fortran)/2.)*(np.cos(xkm1*th_rd_s)*rdn_s + np.sin(xkm1*th_rd_s)*tdn_s)
        
        # Place blocks in the first row/column
        A[0:NF, i_dist*NF:(i_dist+1)*NF] = A_block
        B[0:NF, i_dist*NF:(i_dist+1)*NF] = B_block
        C[0:NF, i_dist*NF:(i_dist+1)*NF] = C_block
        D[0:NF, i_dist*NF:(i_dist+1)*NF] = D_block
        A[i_dist*NF:(i_dist+1)*NF, 0:NF] = A_block_sym
        B[i_dist*NF:(i_dist+1)*NF, 0:NF] = B_block_sym
        C[i_dist*NF:(i_dist+1)*NF, 0:NF] = C_block_sym
        D[i_dist*NF:(i_dist+1)*NF, 0:NF] = D_block_sym


    # --- Step 3: Populate the entire matrix based on block-Toeplitz structure ---
    for i_blk in range(N):
        for j_blk in range(N):
            dist = abs(i_blk - j_blk)
            if dist == 0: # Diagonal block
                start = i_blk * NF
                end = (i_blk + 1) * NF
                A[start:end, start:end] = A1
                B[start:end, start:end] = B1
                C[start:end, start:end] = C1
                D[start:end, start:end] = D1
            else: # Off-diagonal block
                start_row, end_row = i_blk * NF, (i_blk + 1) * NF
                start_col, end_col = j_blk * NF, (j_blk + 1) * NF
                
                # Copy from first row or first column depending on relative position
                if i_blk < j_blk: # Upper triangle of blocks, copy from Block (0, dist)
                    src_start_col = dist * NF
                    src_end_col = (dist + 1) * NF
                    A[start_row:end_row, start_col:end_col] = A[0:NF, src_start_col:src_end_col]
                    B[start_row:end_row, start_col:end_col] = B[0:NF, src_start_col:src_end_col]
                    C[start_row:end_row, start_col:end_col] = C[0:NF, src_start_col:src_end_col]
                    D[start_row:end_row, start_col:end_col] = D[0:NF, src_start_col:src_end_col]
                else: # Lower triangle of blocks, copy from Block (dist, 0)
                    src_start_row = dist * NF
                    src_end_row = (dist + 1) * NF
                    A[start_row:end_row, start_col:end_col] = A[src_start_row:src_end_row, 0:NF]
                    B[start_row:end_row, start_col:end_col] = B[src_start_row:src_end_row, 0:NF]
                    C[start_row:end_row, start_col:end_col] = C[src_start_row:src_end_row, 0:NF]
                    D[start_row:end_row, start_col:end_col] = D[src_start_row:src_end_row, 0:NF]


    # --- Matrix calculations ---
    try:
        A_inv = np.linalg.inv(A)
        CGEN_0 = np.zeros((N, N), dtype=dtype)
        for i in range(N):
            for j in range(N):
                CGEN_0[i, j] = 2. * RW * PI * np.sum(A_inv[i*NF, j*NF:(j+1)*NF])

        SUMCAP_0 = np.sum(CGEN_0)
        NM1 = N - 1
        CAP0 = np.zeros((NM1, NM1), dtype=dtype)
        rows = [i for i in range(N) if i != IREF - 1]
        for i_idx, i in enumerate(rows):
            for j_idx, j in enumerate(rows):
                SUMR = np.sum(CGEN_0[i, :])
                SUMC = np.sum(CGEN_0[:, j])
                CAP0[i_idx, j_idx] = CGEN_0[i, j] - (SUMR * SUMC) / SUMCAP_0
        INDUCT = (1 / V02) * np.linalg.inv(CAP0)

        D_inv = np.linalg.inv(D)
        D_inv_C = D_inv @ C
        A_eff = A - B @ D_inv_C
        A_eff_inv = np.linalg.inv(A_eff)
        A_final = A_eff_inv
        B_final = -D_inv_C @ A_final
        CGEN = np.zeros((N, N), dtype=dtype)
        for i in range(N):
            for j in range(N):
                SUM1 = np.sum(A_final[i*NF, j*NF:(j+1)*NF])
                SUM2 = np.sum(B_final[i*NF, j*NF:(j+1)*NF])
                CGEN[i, j] = 2. * RW * PI * SUM1 + 2. * RD * PI * SUM2

        CAP = np.zeros((NM1, NM1), dtype=dtype)
        SUMCAP = np.sum(CGEN)
        for i_idx, i in enumerate(rows):
            for j_idx, j in enumerate(rows):
                SUMR = np.sum(CGEN[i, :])
                SUMC = np.sum(CGEN[:, j])
                CAP[i_idx, j_idx] = CGEN[i, j] - (SUMR * SUMC) / SUMCAP

    except np.linalg.LinAlgError as e:
        print(f"Error: Singular matrix encountered during linear algebra computation: {e}")
        return

    # --- Write the output output_filename 'PUL.TXT' ---
    with open(output_filename, 'w') as f:
        for i in range(NM1):
            for j in range(i, NM1):
                f.write(f"{i+1:3d}  {j+1:3d}  {INDUCT[i, j]:12.5E}        =L({i+1:3d},{j+1:3d})\n")
        # f.write("\n")
        for i in range(NM1):
            for j in range(i, NM1):
                f.write(f"{i+1:3d}  {j+1:3d}  {CAP[i, j]:12.5E}        =C({i+1:3d},{j+1:3d})\n")
        # f.write("\n")
        for i in range(NM1):
            for j in range(i, NM1):
                f.write(f"{i+1:3d}  {j+1:3d}  {CAP0[i, j]:12.5E}        =C0({i+1:3d},{j+1:3d})\n")
        # f.write("\n")
        for i in range(N):
            for j in range(i, N):
                f.write(f"{i+1:3d}  {j+1:3d}  {CGEN[i, j]:12.5E}        =CGEN({i+1:3d},{j+1:3d})\n")

        f.write("\n\n\n")
        f.write(f"NUMBER OF WIRES= {N:3d}\n")
        f.write(f"NUMBER OF FOURIER COEFFICIENTS= {NF:3d}\n")
        f.write(f"REFERENCE WIRE= {IREF:3d}\n")
        f.write(f"WIRE RADIUS (m)= {RW:10.3E}\n")
        f.write(f"DIELECTRIC INSULATION THICKNESS (m)= {TD:10.3E}\n")
        f.write(f"DIELECTRIC CONSTANT OF INSULATION= {ER:10.3E}\n")
        f.write(f"CENTER-TO-CENTER SEPARATION (m)= {S:10.3E}\n")

    print(f"Analysis complete. Results saved to {output_filename}.")

    # --- Final Validation Step ---
    # validate_results(output_filename, reference_output_file)

if __name__ == '__main__':
    os.system('cls' if os.name == 'nt' else 'clear')
    ribbon_dot_for()