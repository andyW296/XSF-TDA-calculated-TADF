"""Compare the candidate state-selection rules for the USF-TD-DFT column U."""
import mr_tadf_lib as L

TD, SF, SFSTATE, USTD, inv = L.collect()
df, names = L.load_excel()

ok = {"smaller": 0, "larger": 0, "state3": 0}
rows = []
for m in sorted(USTD.get("BHHLYP", {})):
    ex = df.iat[3 + m - 1, 20]
    states = USTD["BHHLYP"][m][0][0]
    vs, ns, ss = L.pick_state_2_3(states, keep="smaller")
    vl, nl, sl = L.pick_state_2_3(states, keep="larger")
    d = {n: (e, s) for n, e, s in states}
    s3 = d[3][0] if 3 in d else None
    if vs is not None and abs(vs - ex) <= 5e-4:
        ok["smaller"] += 1
    if vl is not None and abs(vl - ex) <= 5e-4:
        ok["larger"] += 1
    if s3 is not None and abs(s3 - ex) <= 5e-4:
        ok["state3"] += 1
    rows.append((m, names[m - 1], float(ex), vs, vl, s3, ns, nl))

print("column U (BHHLYP USF-TD-DFT), 35 molecules")
print("  rule 'smaller <S^2>' among states 2,3 :", ok["smaller"], "match")
print("  rule 'larger  <S^2>' among states 2,3 :", ok["larger"], "match")
print("  rule 'always the 3rd eigenvalue'      :", ok["state3"], "match")
print()
print(f"{'#':>2} {'molecule':<22} {'excel':>9} {'smaller':>9} {'larger':>9} {'state3':>9}"
      f"  {'rule picks':>10}")
for m, nm, ex, vs, vl, s3, ns, nl in rows:
    mark = "" if s3 is not None and abs(s3 - ex) <= 5e-4 else "  <-- excel != state3"
    pick = "state " + str(nl)
    s3s = f"{s3:.5f}" if s3 is not None else "   --   "
    print(f"{m:>2} {str(nm):<22} {ex:>9.5f} {vs:>9.5f} {vl:>9.5f} "
          f"{s3s:>9}  {pick:>10}{mark}")
