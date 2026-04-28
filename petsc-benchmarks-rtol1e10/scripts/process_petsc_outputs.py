from pathlib import Path
import re
import csv
import statistics

OUT = Path.home() / "bc" / "petsc-benchmarks-rtol1e10"
LOGS = OUT / "logs"

CASES = {
    "AS": ("A-S", "2x2", 4, 2, 2, "S", 16, 2, 32, 32),
    "AM": ("A-M", "2x2", 4, 2, 2, "M", 32, 2, 64, 64),
    "AL": ("A-L", "2x2", 4, 2, 2, "L", 48, 2, 96, 96),
    "BS": ("B-S", "4x2", 8, 4, 2, "S", 16, 2, 64, 32),
    "BM": ("B-M", "4x2", 8, 4, 2, "M", 32, 2, 128, 64),
    "BL": ("B-L", "4x2", 8, 4, 2, "L", 48, 2, 192, 96),
    "CS": ("C-S", "4x4", 16, 4, 4, "S", 16, 2, 64, 64),
    "CM": ("C-M", "4x4", 16, 4, 4, "M", 32, 2, 128, 128),
    "CL": ("C-L", "4x4", 16, 4, 4, "L", 48, 2, 192, 192),
}

ORDER = ["AS", "AM", "AL", "BS", "BM", "BL", "CS", "CM", "CL"]


def parse_log(path):
    text = path.read_text(errors="replace")

    bddc_it = int(re.search(
        r"Linear physical_ solve converged due to CONVERGED_RTOL iterations\s+(\d+)",
        text
    ).group(1))

    feti_it = int(re.search(
        r"Linear fluxes_ solve converged due to CONVERGED_RTOL iterations\s+(\d+)",
        text
    ).group(1))

    dofs = [int(x) for x in re.findall(
        r"Number of degrees of freedom\s*:\s*(\d+)",
        text
    )]

    eigs = re.findall(
        r"Eigenvalues preconditioned operator\s*:\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)",
        text
    )

    bddc_lmin, bddc_lmax = map(float, eigs[0])
    feti_lmin, feti_lmax = map(float, eigs[1])

    time_match = re.search(r"^Time \(sec\):\s*([0-9.eE+-]+)", text, re.MULTILINE)
    total_time = float(time_match.group(1))

    return {
        "bddc_dofs": dofs[0],
        "feti_dofs": dofs[1],
        "bddc_iters": bddc_it,
        "feti_iters": feti_it,
        "bddc_lambda_min": bddc_lmin,
        "bddc_lambda_max": bddc_lmax,
        "bddc_kappa": bddc_lmax / bddc_lmin,
        "feti_lambda_min": feti_lmin,
        "feti_lambda_max": feti_lmax,
        "feti_kappa": feti_lmax / feti_lmin,
        "total_time_sec": total_time,
    }


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


raw_rows = []
summary_rows = []

for code in ORDER:
    label, layout, mpi_ranks, npx, npy, size_level, local_res, p, nex, ney = CASES[code]

    runs = []
    for run_id in [1, 2, 3]:
        path = LOGS / f"{code}_run{run_id}.txt"
        data = parse_log(path)
        runs.append(data)

        raw_rows.append({
            "case_id": label,
            "run_id": f"run{run_id}",
            "layout": layout,
            "mpi_ranks": mpi_ranks,
            "npx": npx,
            "npy": npy,
            "size_level": size_level,
            "local_res": local_res,
            "p": p,
            "nex": nex,
            "ney": ney,
            **data,
        })

    first = runs[0]

    summary_rows.append({
        "case_id": label,
        "layout": layout,
        "mpi_ranks": mpi_ranks,
        "npx": npx,
        "npy": npy,
        "size_level": size_level,
        "local_res": local_res,
        "p": p,
        "nex": nex,
        "ney": ney,
        "bddc_dofs": first["bddc_dofs"],
        "feti_dofs": first["feti_dofs"],
        "bddc_iters": first["bddc_iters"],
        "feti_iters": first["feti_iters"],
        "bddc_lambda_min": first["bddc_lambda_min"],
        "bddc_lambda_max": first["bddc_lambda_max"],
        "bddc_kappa": first["bddc_kappa"],
        "feti_lambda_min": first["feti_lambda_min"],
        "feti_lambda_max": first["feti_lambda_max"],
        "feti_kappa": first["feti_kappa"],
        "median_total_time_sec": statistics.median(r["total_time_sec"] for r in runs),
    })


raw_fields = list(raw_rows[0].keys())
summary_fields = list(summary_rows[0].keys())

write_csv(OUT / "tables" / "raw_results.csv", raw_rows, raw_fields)
write_csv(OUT / "tables" / "summary_results.csv", summary_rows, summary_fields)
write_csv(OUT / "final" / "tables" / "petsc_main_results.csv", summary_rows, summary_fields)


def extract_residual_csv(case):
    path = LOGS / f"{case}_residuals_rtol1e-10.txt"
    text = path.read_text(errors="replace")

    residual_lines = re.findall(
        r"^\s*(\d+)\s+KSP\s+(preconditioned|none)\s+resid norm\s+[0-9.eE+-]+\s+true resid norm\s+[0-9.eE+-]+\s+\|\|r\(i\)\|\|/\|\|b\|\|\s+([0-9.eE+-]+)",
        text,
        flags=re.MULTILINE
    )

    bddc = []
    feti = []

    for it, norm_type, relres in residual_lines:
        if norm_type == "preconditioned":
            bddc.append((int(it), float(relres)))
        elif norm_type == "none":
            feti.append((int(it), float(relres)))

    out_path = OUT / "figures" / f"{case}_residuals.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["iter_bddc", "bddc_relres", "iter_fetidp", "fetidp_relres"])

        n = max(len(bddc), len(feti))
        for k in range(n):
            row = []
            row += list(bddc[k]) if k < len(bddc) else ["", ""]
            row += list(feti[k]) if k < len(feti) else ["", ""]
            writer.writerow(row)


extract_residual_csv("CM")
extract_residual_csv("CL")

print("Done.")
print(f"Wrote tables to: {OUT / 'final' / 'tables'}")
print(f"Wrote residual CSV files to: {OUT / 'figures'}")
