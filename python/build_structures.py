"""Step 1 (P0): build the Ti-7Al boxes and the B1 / B5a relaxation folders.

    python python/1_build/build_structures.py            # build everything that is missing
    python python/1_build/build_structures.py --force    # rebuild and overwrite

Each run folder gets: library.meam, AlTi.meam, box.data, in.relax
A summary of all boxes is written to results/structures.csv
"""

import argparse
import csv
import shutil
import sys

import numpy as np

from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))

import settings as S
from sro import hcp_positions, neighbor_list, random_types, swap_mc, warren_cowley


def write_data(path, pos, box, types, title):
    with open(path, "w") as f:
        f.write(f"{title}\n\n")
        f.write(f"{len(pos)} atoms\n")
        f.write("2 atom types\n\n")
        f.write(f"0.0 {box[0]:.6f} xlo xhi\n")
        f.write(f"0.0 {box[1]:.6f} ylo yhi\n")
        f.write(f"0.0 {box[2]:.6f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write(f"1 {S.MASS_TI}  # Ti\n")
        f.write(f"2 {S.MASS_AL}  # Al\n\n")
        f.write("Atoms  # atomic\n\n")
        for k, (t, (x, y, z)) in enumerate(zip(types, pos), start=1):
            f.write(f"{k} {t} {x:.6f} {y:.6f} {z:.6f}\n")


def write_input(template, path, subs):
    text = template
    for key, val in subs.items():
        text = text.replace(f"@{key}@", str(val))
    if "@" in text:
        raise ValueError(f"unfilled placeholder in {path}")
    path.write_text(text)


def build_set(case, cell, n_boxes, seed_base, batch, template, rows, force):
    pos, box = hcp_positions(cell, S.A_LAT, S.C_LAT)
    n_sites = len(pos)
    n_al = int(round(S.X_AL * n_sites))
    nbr1 = neighbor_list(cell, S.A_LAT, S.C_LAT, 0.1, S.RC1)
    nbr2 = neighbor_list(cell, S.A_LAT, S.C_LAT, S.RC1, S.RC2)
    target = S.ALPHA_TARGET[case]

    for k in range(1, n_boxes + 1):
        name = f"box{k:02d}"
        run_dir = S.SIM_DIR / case / batch / name
        if run_dir.exists() and not force:
            print(f"  {case}/{batch}/{name}: exists, skipped")
            continue
        run_dir.mkdir(parents=True, exist_ok=True)

        struct_seed = seed_base + k
        vel_seed = struct_seed + S.VEL_OFFSET
        rng = np.random.default_rng(struct_seed)

        types = random_types(n_sites, n_al, rng)
        types, attempts, accepted = swap_mc(types, nbr1, target, S.ALPHA_TOL, rng)

        a1 = warren_cowley(types, nbr1)
        a2 = warren_cowley(types, nbr2)

        title = f"Ti-7Al {case} {name} (seed {struct_seed}, alpha1 {a1:.4f})"
        write_data(run_dir / "box.data", pos, box, types, title)
        write_input(template, run_dir / "in.relax", {
            "CASE": case,
            "BOX": name,
            "LIB_ELEMENTS": S.LIB_ELEMENTS,
            "TYPE_ELEMENTS": S.TYPE_ELEMENTS,
            "VEL_SEED": vel_seed,
        })
        for fname in S.POT_FILES:
            shutil.copy(S.POT_DIR / fname, run_dir / fname)

        rows.append({
            "case": case, "batch": batch, "box": name,
            "nx": cell[0], "ny": cell[1], "nz": cell[2],
            "lx": f"{box[0]:.3f}", "ly": f"{box[1]:.3f}", "lz": f"{box[2]:.3f}",
            "natoms": n_sites, "n_al": n_al, "x_al": f"{n_al / n_sites:.4f}",
            "alpha1": f"{a1:.4f}", "alpha2": f"{a2:.4f}",
            "struct_seed": struct_seed, "vel_seed": vel_seed,
            "mc_attempts": attempts, "mc_accepted": accepted,
        })
        print(f"  {case}/{batch}/{name}: alpha1 = {a1:+.4f}  alpha2 = {a2:+.4f}")


def main():
    parser = argparse.ArgumentParser(description="Build Ti-7Al boxes and relaxation folders")
    parser.add_argument("--force", action="store_true", help="overwrite existing folders")
    args = parser.parse_args()

    missing = [f for f in S.POT_FILES if not (S.POT_DIR / f).is_file()]
    if missing:
        sys.exit(f"Missing potential files in {S.POT_DIR}: {', '.join(missing)}")

    template = (S.TEMPLATE_DIR / "in.relax.template").read_text()
    S.RESULTS_DIR.mkdir(exist_ok=True)
    rows = []

    for case in S.ALPHA_TARGET:
        print(f"{case}: main boxes {S.MAIN_CELL}")
        build_set(case, S.MAIN_CELL, S.N_BOXES, S.SEED_BASE[case],
                  "B1_relax", template, rows, args.force)
        if case in S.LARGE_CASES:
            print(f"{case}: large boxes {S.LARGE_CELL}")
            build_set(case, S.LARGE_CELL, S.N_BOXES_LARGE, S.SEED_BASE_LARGE[case],
                      "B5a_relax80", template, rows, args.force)

    if rows:
        out = S.RESULTS_DIR / "structures.csv"
        new_file = not out.exists() or args.force
        with open(out, "w" if new_file else "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            if new_file:
                w.writeheader()
            w.writerows(rows)
        print(f"\n{len(rows)} boxes written, summary in {out}")
    else:
        print("\nNothing to do.")


if __name__ == "__main__":
    main()
