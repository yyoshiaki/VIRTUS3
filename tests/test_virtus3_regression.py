from argparse import Namespace
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import virtus3  # noqa: E402


TG_MAP = REPO_ROOT / "data" / "NC_007605.1_CDS_EBER12.tgMap.tsv"
SCRATCH_ROOT = Path("/vast/palmer/scratch/hafler/yy693/PRJNA1301449/virtus3_runs")


def build_args(sample, fastqs_dir, output_dir):
    return Namespace(
        fastqs=str(fastqs_dir),
        chemistry_cr="ARC-v1",
        sample=sample,
        lib_alevin="-l ISR --chromiumV3",
        output=str(output_dir),
        index_human="/tmp/refdata",
        index_virus=str(REPO_ROOT / "data" / "NC_007605.1_CDS_EBER12_salmon_index"),
        tgMap=str(TG_MAP),
        cellranger="/tmp/cellranger",
        salmon="/tmp/salmon",
        cores=4,
        mem_per_core=1,
        skip_exist=True,
        use_filtered_bc=False,
        expect_cells=None,
    )


def make_dummy_fastqs(fastqs_dir, sample):
    fastqs_dir.mkdir(parents=True, exist_ok=True)
    for read in ["R1", "R2"]:
        (fastqs_dir / f"{sample}_S1_L001_{read}_001.fastq.gz").write_bytes(b"")


def copy_file(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def prepare_minimal_outs(sample, output_dir, include_stale_lane=False):
    source_root = SCRATCH_ROOT / sample / "cellranger_human" / "outs"
    target_root = output_dir / "cellranger_human" / "outs"
    target_root.mkdir(parents=True, exist_ok=True)

    # Files used only as existence checks can be empty placeholders.
    (target_root / "possorted_genome_bam.bam").write_bytes(b"")
    (target_root / "unmapped.bam").write_bytes(b"placeholder")
    (target_root / "raw_feature_bc_matrix").mkdir(exist_ok=True)
    (target_root / "raw_feature_bc_matrix" / "barcodes.tsv").write_text("AAACCCAAGAAACACT\n")
    unmapped_dir = target_root / "unmapped_fqs" / "cellranger_human_0_1_unknow_flowcell"
    unmapped_dir.mkdir(parents=True, exist_ok=True)
    for read in ["R1", "R2"]:
        (unmapped_dir / f"bamtofastq_S1_L000_{read}_001.fastq.gz").write_bytes(b"")

    copy_file(
        source_root / "alevin_virus" / "logs" / "salmon_quant.log",
        target_root / "alevin_virus" / "logs" / "salmon_quant.log",
    )
    copy_file(
        source_root / "alevin_virus" / "cmd_info.json",
        target_root / "alevin_virus" / "cmd_info.json",
    )

    if sample == "SRR35978549":
        for relative in [
            ("alevin_virus/alevin/quants_mat.mtx.gz", "alevin_virus/alevin/quants_mat.mtx.gz"),
            ("alevin_virus/alevin/quants_mat_rows.txt", "alevin_virus/alevin/quants_mat_rows.txt"),
            ("alevin_virus/alevin/quants_mat_cols.txt", "alevin_virus/alevin/quants_mat_cols.txt"),
        ]:
            copy_file(source_root / relative[0], target_root / relative[1])
    elif sample == "SRR35978576":
        for relative in [
            ("alevin_virus/alevin/alevin.log", "alevin_virus/alevin/alevin.log"),
            ("alevin_virus/alevin/quants_mat.mtx.gz", "alevin_virus/alevin/quants_mat.mtx.gz"),
            ("alevin_virus/alevin/quants_mat_cols.txt", "alevin_virus/alevin/quants_mat_cols.txt"),
        ]:
            copy_file(source_root / relative[0], target_root / relative[1])
        if include_stale_lane:
            copy_file(
                source_root / "alevin_virus_lane_001" / "logs" / "salmon_quant.log",
                target_root / "alevin_virus_lane_001" / "logs" / "salmon_quant.log",
            )
    else:
        raise ValueError(f"Unsupported sample for regression fixture: {sample}")

    return target_root


@pytest.fixture(autouse=True)
def restore_cwd():
    original = Path.cwd()
    try:
        yield
    finally:
        os.chdir(original)


@pytest.fixture
def stub_run_command(monkeypatch):
    monkeypatch.setattr(virtus3, "run_command", lambda command: "cellranger-9.0.1\n")


def test_pipeline_complete_with_reads(tmp_path, stub_run_command):
    sample = "SRR35978549"
    output_dir = tmp_path / "output"
    prepare_minimal_outs(sample, output_dir)
    fastqs_dir = tmp_path / "fastqs"
    make_dummy_fastqs(fastqs_dir, sample)

    args = build_args(sample, fastqs_dir, output_dir)
    log = virtus3.pipeline(args)

    csv_path = output_dir / "alevin_virus.csv"
    assert csv_path.exists()
    df = pd.read_csv(csv_path, index_col=0)
    assert df.shape[0] > 0
    assert df.to_numpy().sum() > 0
    assert "Total viral UMIs" in log


def test_pipeline_zero_read_partial_output_returns_empty_matrix(tmp_path, stub_run_command):
    sample = "SRR35978576"
    output_dir = tmp_path / "output"
    prepare_minimal_outs(sample, output_dir)
    fastqs_dir = tmp_path / "fastqs"
    make_dummy_fastqs(fastqs_dir, sample)

    args = build_args(sample, fastqs_dir, output_dir)
    log = virtus3.pipeline(args)

    csv_path = output_dir / "alevin_virus.csv"
    assert csv_path.exists()
    df = pd.read_csv(csv_path, index_col=0)
    expected_features = pd.read_csv(TG_MAP, sep="\t", header=None, usecols=[1]).drop_duplicates().shape[0]
    assert df.shape == (0, expected_features)
    assert "Total viral UMIs" in log


def test_stale_lane_output_is_ignored(tmp_path, stub_run_command):
    sample = "SRR35978576"
    output_dir = tmp_path / "output"
    outs_dir = prepare_minimal_outs(sample, output_dir, include_stale_lane=True)
    stale_log = outs_dir / "alevin_virus_lane_001" / "logs" / "salmon_quant.log"
    stale_log.write_text("BROKEN STALE LOG\n")
    fastqs_dir = tmp_path / "fastqs"
    make_dummy_fastqs(fastqs_dir, sample)

    args = build_args(sample, fastqs_dir, output_dir)
    virtus3.pipeline(args)

    csv_path = output_dir / "alevin_virus.csv"
    assert csv_path.exists()


def test_missing_cellranger_outs_raises_meaningful_error(tmp_path, monkeypatch, stub_run_command):
    sample = "SRR35978545"
    fastqs_dir = tmp_path / "fastqs"
    make_dummy_fastqs(fastqs_dir, sample)
    output_dir = tmp_path / "output"

    monkeypatch.setattr(
        virtus3,
        "run_command_result",
        lambda command: subprocess.CompletedProcess(command, 1, stdout="", stderr="cellranger failed"),
    )

    args = build_args(sample, fastqs_dir, output_dir)
    args.skip_exist = False

    with pytest.raises(RuntimeError, match="cellranger count failed and did not create the expected output directory"):
        virtus3.pipeline(args)
