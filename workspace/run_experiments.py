"""Batch experiment for MRVQ-copy.

Runs MRVQ-copy.exe on every test image for each (NB, LAMBDA) pair and
writes all results to results/results.csv (one row per run).

Usage (run from the workspace folder, MRCB_8192.raw must be here):
    gcc -O2 -o MRVQ-copy.exe MRVQ-copy.c -lm
    python run_experiments.py
"""
import csv
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
EXE = HERE / ("MRVQ-copy.exe" if sys.platform == "win32" else "MRVQ-copy")
IMAGE_DIR = HERE.parent / "Test IMAGE"
OUT_DIR = HERE / "results"

NB_LIST = [2, 4]
BM_LIST = [0, 1]
LAMBDA_LIST = [0, 10, 20, 30, 50, 75, 100, 150, 200, 300]

FIELDS = ["image", "NB", "BM_on", "LAMBDA",
          "MRVQ", "LEFT", "UP", "UL", "UR", "BM", "bpp", "PSNR"]


def run_one(image: pathlib.Path, nb: int, bm: int, lam: float) -> dict:
    out_img = OUT_DIR / f"{image.stem}_NB{nb}_BM{bm}_L{lam:g}.raw"
    proc = subprocess.run(
        [str(EXE), str(image), str(out_img), str(nb), str(lam), str(bm)],
        cwd=HERE, capture_output=True, text=True, check=True,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("RESULT,"):
            values = line.split(",")[1:]
            row = dict(zip(FIELDS[1:], values))
            row["image"] = image.stem
            return row
    raise RuntimeError(f"no RESULT line for {image.name}:\n{proc.stdout}")


def main() -> None:
    if not EXE.exists():
        sys.exit(f"{EXE.name} not found - compile MRVQ-copy.c first.")
    images = sorted(IMAGE_DIR.glob("*.raw")) + [HERE / "lenna512.raw"]
    OUT_DIR.mkdir(exist_ok=True)

    with open(OUT_DIR / "results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for image in images:
            for nb in NB_LIST:
                for bm in BM_LIST:
                    for lam in LAMBDA_LIST:
                        row = run_one(image, nb, bm, lam)
                        writer.writerow(row)
                        print(f"{image.stem:10s} NB={nb} BM={bm} LAMBDA={lam:<5g} "
                              f"bpp={row['bpp']}  PSNR={row['PSNR']}")
    print(f"\nSaved {OUT_DIR / 'results.csv'}")


if __name__ == "__main__":
    main()
