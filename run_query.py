import re
import pandas as pd

log_path = r"Designed INVEST candidates/XSF-TDA/alda0_b3lyp_pbe0_cam.log"
xlsx_path = r"TADF_all.xlsx"

with open(log_path, encoding="utf-8", errors="ignore") as f:
    lines = f.readlines()

cur_label = None
cur_mol = None
cur_fun = None
label_re = re.compile(r'^(ALDA0|Multicollinear) (\d+)_optimized\.xyz ([\w-]+)$')
ex2_re = re.compile(r'^Excited state 2\s+(-?\d+\.\d+)\s+eV')

results = {}  # (fun,label) -> {mol: val}
for fun in ['b3lyp','pbe0','cam-b3lyp']:
    for label in ['ALDA0','Multicollinear']:
        results[(fun,label)] = {}

for line in lines:
    line = line.rstrip('\n')
    m = label_re.match(line)
    if m:
        cur_label = m.group(1)
        cur_mol = int(m.group(2))
        cur_fun = m.group(3)
        continue
    m = ex2_re.match(line.strip())
    if m and cur_label is not None and cur_fun in ('b3lyp','pbe0','cam-b3lyp'):
        val = float(m.group(1))
        d = results[(cur_fun,cur_label)]
        if cur_mol not in d:
            d[cur_mol] = val
        cur_label = None  # avoid re-capturing subsequent Excited state 2 lines before next label

xl = pd.ExcelFile(xlsx_path)
df = xl.parse('Peng', header=None)
mol_col = df.iloc[3:, 0]

col_map = {
    'G_b3lyp_ALDA0': (6, 'b3lyp', 'ALDA0'),
    'I_b3lyp_MC': (8, 'b3lyp', 'Multicollinear'),
    'O_pbe0_ALDA0': (14, 'pbe0', 'ALDA0'),
    'Q_pbe0_MC': (16, 'pbe0', 'Multicollinear'),
    'AE_cam_ALDA0': (30, 'cam-b3lyp', 'ALDA0'),
    'AG_cam_MC': (32, 'cam-b3lyp', 'Multicollinear'),
}

all_mismatches = {}
for name, (col_idx, fun, label) in col_map.items():
    print(f'=== {name} (excel col idx {col_idx}) vs log ({fun}, {label}) ===')
    excel_col = df.iloc[3:, col_idx]
    mismatches = []
    for idx, mol in mol_col.items():
        if pd.isna(mol):
            continue
        mol = int(mol)
        e_val = excel_col[idx]
        log_val = results[(fun,label)].get(mol)
        if log_val is None:
            print(f'{mol:>4} {e_val:>14}  MISSING')
            mismatches.append((mol,'MISSING'))
            continue
        diff = e_val - log_val
        flag = ' MISMATCH' if abs(diff) > 1e-4 else ''
        print(f'{mol:>4} {e_val:>14} {log_val:>14} {diff:>10.6f}{flag}')
        if flag:
            mismatches.append((mol,flag))
    all_mismatches[name] = mismatches
    print(f'-- {name} mismatch count: {len(mismatches)} {mismatches}')
    print()

print('SUMMARY:')
for name, mism in all_mismatches.items():
    print(name, '->', mism)

print()
print('Excel Peng molecule col A values:', list(mol_col.dropna().astype(int).values))
for fun in ['b3lyp','pbe0','cam-b3lyp']:
    for label in ['ALDA0','Multicollinear']:
        d = results[(fun,label)]
        print(f'{fun} {label} molecules found:', sorted(d.keys()), ' mol10:', d.get(10), 'mol11:', d.get(11))
