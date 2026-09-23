import pandas as pd
xlsx_path = r"TADF_all.xlsx"
xl = pd.ExcelFile(xlsx_path)
with open('sheet_dump.txt','w',encoding='utf-8') as out:
    for name in ['Peng','Fig5','Fig3']:
        df = xl.parse(name, header=None)
        out.write(f'=== SHEET: {name} shape: {df.shape}\n')
        for r in range(min(6, df.shape[0])):
            out.write(f'{r} {list(df.iloc[r, :].values)}\n')
        out.write(f'Molecule names col A: {list(df.iloc[3:, 0].values)}\n\n')
