#!/usr/bin/env python3
"""Merge multiple filterStrainCube outputs into one BSE-compatible state set.

The filter code writes:
  - eval.dat: text, one state per line -> "<index> <eigenvalue> <sigma>"
  - psi.dat: raw binary doubles, no header, one contiguous block per state

The BSE code reads:
  - eval.par: same 3-column text format
  - psi.par: same raw binary-double layout

This script merges multiple filter runs by sorting all states by eigenvalue,
rewriting indices to 0..N-1, and copying the matching wavefunction blocks in
the same order so eigenvalue/state correspondence is preserved.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DOUBLE_SIZE = 8
COPY_CHUNK_SIZE = 1024 * 1024


class MergeError(RuntimeError):
    """Raised when the provided filter outputs cannot be merged safely."""


@dataclass(frozen=True)
class Dataset:
    root: Path
    eval_path: Path
    psi_path: Path
    input_path: Path | None
    conf_path: Path | None
    nstates: int


@dataclass(frozen=True)
class StateRecord:
    energy: float
    sigma: float
    source_index: int
    local_index: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine several filterStrainCube outputs into one eval.par/psi.par "
            "pair compatible with bse/."
        )
    )
    parser.add_argument(
        "datasets",
        nargs="+",
        help="Directories containing eval.dat and psi.dat from filterStrainCube runs.",
    )
    parser.add_argument(
        "--output-dir",
        default="merged_filter",
        help="Directory to write eval.par, psi.par, and a merge manifest into.",
    )
    parser.add_argument(
        "--grid-source",
        help=(
            "Path to an input.par file whose first three integers are nx ny nz. "
            "This can be either a filter input.par or a bse input.par."
        ),
    )
    parser.add_argument(
        "--ngrid",
        type=int,
        help="Number of real-space grid points per state. Use instead of --grid-source.",
    )
    parser.add_argument(
        "--eval-name",
        default="eval.dat",
        help="Filename for the filter eigenvalue file inside each dataset directory.",
    )
    parser.add_argument(
        "--psi-name",
        default="psi.dat",
        help="Filename for the filter wavefunction file inside each dataset directory.",
    )
    parser.add_argument(
        "--copy-conf",
        action="store_true",
        help="Copy conf.par from the first dataset into the output directory after verifying all copies match.",
    )
    args = parser.parse_args()

    if (args.grid_source is None) == (args.ngrid is None):
        parser.error("provide exactly one of --grid-source or --ngrid")

    if args.ngrid is not None and args.ngrid <= 0:
        parser.error("--ngrid must be positive")

    return args


def read_ngrid_from_input(path: Path) -> int:
    try:
        tokens = path.read_text().split()
    except OSError as exc:
        raise MergeError(f"failed to read grid source {path}: {exc}") from exc

    if len(tokens) < 3:
        raise MergeError(f"grid source {path} does not contain nx ny nz")

    try:
        nx = int(tokens[0])
        ny = int(tokens[1])
        nz = int(tokens[2])
    except ValueError as exc:
        raise MergeError(f"grid source {path} does not start with integer nx ny nz") from exc

    if nx <= 0 or ny <= 0 or nz <= 0:
        raise MergeError(f"grid source {path} contains a non-positive grid size")

    return nx * ny * nz


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(COPY_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_eval_file(path: Path) -> list[tuple[int, float, float]]:
    rows: list[tuple[int, float, float]] = []
    try:
        with path.open("r", encoding="utf-8") as handle:
            for lineno, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped:
                    continue
                parts = stripped.split()
                if len(parts) < 3:
                    raise MergeError(f"{path}:{lineno} does not have 3 columns")
                try:
                    idx = int(parts[0])
                    energy = float(parts[1])
                    sigma = float(parts[2])
                except ValueError as exc:
                    raise MergeError(f"{path}:{lineno} has a non-numeric field") from exc
                rows.append((idx, energy, sigma))
    except OSError as exc:
        raise MergeError(f"failed to read {path}: {exc}") from exc

    if not rows:
        raise MergeError(f"{path} is empty")

    return rows


def discover_dataset(root: Path, eval_name: str, psi_name: str, ngrid: int) -> Dataset:
    eval_path = root / eval_name
    psi_path = root / psi_name
    input_path = root / "input.par"
    conf_path = root / "conf.par"

    if not eval_path.is_file():
        raise MergeError(f"missing {eval_path}")
    if not psi_path.is_file():
        raise MergeError(f"missing {psi_path}")

    eval_rows = parse_eval_file(eval_path)
    nstates = len(eval_rows)

    if input_path.is_file():
        dataset_ngrid = read_ngrid_from_input(input_path)
        if dataset_ngrid != ngrid:
            raise MergeError(
                f"{input_path} implies ngrid={dataset_ngrid}, expected {ngrid}"
            )

    psi_size = psi_path.stat().st_size
    expected_size = nstates * ngrid * DOUBLE_SIZE
    if psi_size != expected_size:
        raise MergeError(
            f"{psi_path} size is {psi_size} bytes, expected {expected_size} "
            f"for {nstates} states and ngrid={ngrid}"
        )

    expected_indices = list(range(nstates))
    actual_indices = [row[0] for row in eval_rows]
    if actual_indices != expected_indices:
        raise MergeError(
            f"{eval_path} indices are not sequential 0..{nstates - 1}; "
            "the psi block mapping would be ambiguous"
        )

    return Dataset(
        root=root,
        eval_path=eval_path,
        psi_path=psi_path,
        input_path=input_path if input_path.is_file() else None,
        conf_path=conf_path if conf_path.is_file() else None,
        nstates=nstates,
    )


def iter_state_records(datasets: list[Dataset]) -> list[StateRecord]:
    records: list[StateRecord] = []
    for source_index, dataset in enumerate(datasets):
        for local_index, (_, energy, sigma) in enumerate(parse_eval_file(dataset.eval_path)):
            records.append(
                StateRecord(
                    energy=energy,
                    sigma=sigma,
                    source_index=source_index,
                    local_index=local_index,
                )
            )
    records.sort(
        key=lambda item: (item.energy, item.sigma, item.source_index, item.local_index)
    )
    return records


def ensure_matching_conf(datasets: list[Dataset]) -> None:
    available = [dataset for dataset in datasets if dataset.conf_path is not None]
    if len(available) <= 1:
        return

    reference = sha256(available[0].conf_path)
    for dataset in available[1:]:
        current = sha256(dataset.conf_path)
        if current != reference:
            raise MergeError(
                "conf.par files differ between datasets; these filter outputs do not "
                "describe the same nanostructure and cannot be merged for BSE"
            )


def copy_state_block(src_handle, dst_handle, block_size: int) -> None:
    remaining = block_size
    while remaining > 0:
        chunk = src_handle.read(min(COPY_CHUNK_SIZE, remaining))
        if not chunk:
            raise MergeError("unexpected end of psi.dat while copying a state block")
        dst_handle.write(chunk)
        remaining -= len(chunk)


def write_outputs(
    output_dir: Path,
    datasets: list[Dataset],
    records: list[StateRecord],
    ngrid: int,
    copy_conf: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    eval_out = output_dir / "eval.par"
    psi_out = output_dir / "psi.par"
    manifest_out = output_dir / "merge_manifest.txt"
    block_size = ngrid * DOUBLE_SIZE

    source_handles = [dataset.psi_path.open("rb") for dataset in datasets]
    try:
        with eval_out.open("w", encoding="utf-8") as eval_handle, psi_out.open("wb") as psi_handle:
            for merged_index, record in enumerate(records):
                eval_handle.write(
                    f"{merged_index} {record.energy:.16g} {record.sigma:.16g}\n"
                )
                src = source_handles[record.source_index]
                src.seek(record.local_index * block_size)
                copy_state_block(src, psi_handle, block_size)
    finally:
        for handle in source_handles:
            handle.close()

    with manifest_out.open("w", encoding="utf-8") as handle:
        handle.write(f"ngrid {ngrid}\n")
        handle.write(f"nstates {len(records)}\n")
        for source_index, dataset in enumerate(datasets):
            handle.write(
                f"source {source_index} {dataset.root} states={dataset.nstates}\n"
            )
        for merged_index, record in enumerate(records):
            handle.write(
                "merged_state "
                f"{merged_index} source={record.source_index} local_index={record.local_index} "
                f"energy={record.energy:.16g} sigma={record.sigma:.16g}\n"
            )

    if copy_conf:
        first_conf = datasets[0].conf_path
        if first_conf is None:
            raise MergeError("--copy-conf was requested, but the first dataset has no conf.par")
        shutil.copy2(first_conf, output_dir / "conf.par")


def summarize_sources(datasets: Iterable[Dataset]) -> str:
    return ", ".join(f"{dataset.root} ({dataset.nstates} states)" for dataset in datasets)


def main() -> int:
    args = parse_args()
    ngrid = args.ngrid if args.ngrid is not None else read_ngrid_from_input(Path(args.grid_source))

    try:
        datasets = [
            discover_dataset(Path(dataset).resolve(), args.eval_name, args.psi_name, ngrid)
            for dataset in args.datasets
        ]
        ensure_matching_conf(datasets)
        records = iter_state_records(datasets)
        output_dir = Path(args.output_dir).resolve()
        write_outputs(output_dir, datasets, records, ngrid, args.copy_conf)
    except MergeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"Merged {len(datasets)} datasets: {summarize_sources(datasets)}")
    print(f"Total merged states: {len(records)}")
    print(f"ngrid per state: {ngrid}")
    print(f"Wrote {output_dir / 'eval.par'}")
    print(f"Wrote {output_dir / 'psi.par'}")
    print(f"Wrote {output_dir / 'merge_manifest.txt'}")
    if args.copy_conf:
        print(f"Wrote {output_dir / 'conf.par'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
