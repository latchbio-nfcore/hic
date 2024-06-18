from dataclasses import dataclass
from enum import Enum
import os
import subprocess
import requests
import shutil
from pathlib import Path
import typing
import typing_extensions

from latch.resources.workflow import workflow
from latch.resources.tasks import nextflow_runtime_task, custom_task
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.ldata.path import LPath
from latch_cli.nextflow.workflow import get_flag
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.utils import urljoins
from latch.types import metadata
from flytekit.core.annotation import FlyteAnnotation

from latch_cli.services.register.utils import import_module_by_path

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata

@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_gib": 100,
        }
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]






@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(pvc_name: str, input: LatchFile, outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], email: typing.Optional[str], multiqc_title: typing.Optional[str], genome: typing.Optional[str], fasta: typing.Optional[LatchFile], bwt2_index: typing.Optional[str], digestion: typing.Optional[str], restriction_site: typing.Optional[str], ligation_site: typing.Optional[str], chromosome_size: typing.Optional[LatchFile], restriction_fragments: typing.Optional[LatchFile], save_reference: typing.Optional[bool], dnase: typing.Optional[bool], min_cis_dist: typing.Optional[int], split_fastq: typing.Optional[bool], save_aligned_intermediates: typing.Optional[bool], keep_dups: typing.Optional[bool], keep_multi: typing.Optional[bool], max_insert_size: typing.Optional[int], min_insert_size: typing.Optional[int], max_restriction_fragment_size: typing.Optional[int], min_restriction_fragment_size: typing.Optional[int], save_interaction_bam: typing.Optional[bool], save_pairs_intermediates: typing.Optional[bool], hicpro_maps: typing.Optional[bool], ice_filter_high_count_perc: typing.Optional[int], save_raw_maps: typing.Optional[bool], skip_maps: typing.Optional[bool], skip_dist_decay: typing.Optional[bool], skip_tads: typing.Optional[bool], skip_compartments: typing.Optional[bool], skip_balancing: typing.Optional[bool], skip_mcool: typing.Optional[bool], skip_multiqc: typing.Optional[bool], multiqc_methods_description: typing.Optional[str], fastq_chunks_size: typing.Optional[int], min_mapq: typing.Optional[int], bwt2_opts_end2end: typing.Optional[str], bwt2_opts_trimmed: typing.Optional[str], bin_size: typing.Optional[str], ice_filter_low_count_perc: typing.Optional[float], ice_eps: typing.Optional[float], ice_max_iter: typing.Optional[int], res_zoomify: typing.Optional[str], res_dist_decay: typing.Optional[str], tads_caller: typing.Optional[str], res_tads: typing.Optional[str], res_compartments: typing.Optional[str]) -> None:
    try:
        shared_dir = Path("/nf-workdir")



        ignore_list = [
            "latch",
            ".latch",
            "nextflow",
            ".nextflow",
            "work",
            "results",
            "miniconda",
            "anaconda3",
            "mambaforge",
        ]

        shutil.copytree(
            Path("/root"),
            shared_dir,
            ignore=lambda src, names: ignore_list,
            ignore_dangling_symlinks=True,
            dirs_exist_ok=True,
        )

        cmd = [
            "/root/nextflow",
            "run",
            str(shared_dir / "main.nf"),
            "-work-dir",
            str(shared_dir),
            "-profile",
            "docker",
            "-c",
            "latch.config",
                *get_flag('input', input),
                *get_flag('outdir', outdir),
                *get_flag('email', email),
                *get_flag('multiqc_title', multiqc_title),
                *get_flag('genome', genome),
                *get_flag('fasta', fasta),
                *get_flag('bwt2_index', bwt2_index),
                *get_flag('digestion', digestion),
                *get_flag('restriction_site', restriction_site),
                *get_flag('ligation_site', ligation_site),
                *get_flag('chromosome_size', chromosome_size),
                *get_flag('restriction_fragments', restriction_fragments),
                *get_flag('save_reference', save_reference),
                *get_flag('dnase', dnase),
                *get_flag('min_cis_dist', min_cis_dist),
                *get_flag('split_fastq', split_fastq),
                *get_flag('fastq_chunks_size', fastq_chunks_size),
                *get_flag('min_mapq', min_mapq),
                *get_flag('bwt2_opts_end2end', bwt2_opts_end2end),
                *get_flag('bwt2_opts_trimmed', bwt2_opts_trimmed),
                *get_flag('save_aligned_intermediates', save_aligned_intermediates),
                *get_flag('keep_dups', keep_dups),
                *get_flag('keep_multi', keep_multi),
                *get_flag('max_insert_size', max_insert_size),
                *get_flag('min_insert_size', min_insert_size),
                *get_flag('max_restriction_fragment_size', max_restriction_fragment_size),
                *get_flag('min_restriction_fragment_size', min_restriction_fragment_size),
                *get_flag('save_interaction_bam', save_interaction_bam),
                *get_flag('save_pairs_intermediates', save_pairs_intermediates),
                *get_flag('bin_size', bin_size),
                *get_flag('hicpro_maps', hicpro_maps),
                *get_flag('ice_filter_low_count_perc', ice_filter_low_count_perc),
                *get_flag('ice_filter_high_count_perc', ice_filter_high_count_perc),
                *get_flag('ice_eps', ice_eps),
                *get_flag('ice_max_iter', ice_max_iter),
                *get_flag('res_zoomify', res_zoomify),
                *get_flag('save_raw_maps', save_raw_maps),
                *get_flag('res_dist_decay', res_dist_decay),
                *get_flag('tads_caller', tads_caller),
                *get_flag('res_tads', res_tads),
                *get_flag('res_compartments', res_compartments),
                *get_flag('skip_maps', skip_maps),
                *get_flag('skip_dist_decay', skip_dist_decay),
                *get_flag('skip_tads', skip_tads),
                *get_flag('skip_compartments', skip_compartments),
                *get_flag('skip_balancing', skip_balancing),
                *get_flag('skip_mcool', skip_mcool),
                *get_flag('skip_multiqc', skip_multiqc),
                *get_flag('multiqc_methods_description', multiqc_methods_description)
        ]

        print("Launching Nextflow Runtime")
        print(' '.join(cmd))
        print(flush=True)

        env = {
            **os.environ,
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms2048M -Xmx8G -XX:ActiveProcessorCount=4",
            "K8S_STORAGE_CLAIM_NAME": pvc_name,
            "NXF_DISABLE_CHECK_LATEST": "true",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(urljoins("latch:///your_log_dir/nf_nf_core_hic", name, "nextflow.log"))
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)



@workflow(metadata._nextflow_metadata)
def nf_nf_core_hic(input: LatchFile, outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], email: typing.Optional[str], multiqc_title: typing.Optional[str], genome: typing.Optional[str], fasta: typing.Optional[LatchFile], bwt2_index: typing.Optional[str], digestion: typing.Optional[str], restriction_site: typing.Optional[str], ligation_site: typing.Optional[str], chromosome_size: typing.Optional[LatchFile], restriction_fragments: typing.Optional[LatchFile], save_reference: typing.Optional[bool], dnase: typing.Optional[bool], min_cis_dist: typing.Optional[int], split_fastq: typing.Optional[bool], save_aligned_intermediates: typing.Optional[bool], keep_dups: typing.Optional[bool], keep_multi: typing.Optional[bool], max_insert_size: typing.Optional[int], min_insert_size: typing.Optional[int], max_restriction_fragment_size: typing.Optional[int], min_restriction_fragment_size: typing.Optional[int], save_interaction_bam: typing.Optional[bool], save_pairs_intermediates: typing.Optional[bool], hicpro_maps: typing.Optional[bool], ice_filter_high_count_perc: typing.Optional[int], save_raw_maps: typing.Optional[bool], skip_maps: typing.Optional[bool], skip_dist_decay: typing.Optional[bool], skip_tads: typing.Optional[bool], skip_compartments: typing.Optional[bool], skip_balancing: typing.Optional[bool], skip_mcool: typing.Optional[bool], skip_multiqc: typing.Optional[bool], multiqc_methods_description: typing.Optional[str], fastq_chunks_size: typing.Optional[int] = 20000000, min_mapq: typing.Optional[int] = 10, bwt2_opts_end2end: typing.Optional[str] = "'--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder'", bwt2_opts_trimmed: typing.Optional[str] = "'--very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder'", bin_size: typing.Optional[str] = '1000000,500000', ice_filter_low_count_perc: typing.Optional[float] = 0.02, ice_eps: typing.Optional[float] = 0.1, ice_max_iter: typing.Optional[int] = 100, res_zoomify: typing.Optional[str] = '5000', res_dist_decay: typing.Optional[str] = '1000000', tads_caller: typing.Optional[str] = 'hicexplorer,insulation', res_tads: typing.Optional[str] = '40000,20000', res_compartments: typing.Optional[str] = '250000') -> None:
    """
    nf-core/hic

    Sample Description
    """

    pvc_name: str = initialize()
    nextflow_runtime(pvc_name=pvc_name, input=input, outdir=outdir, email=email, multiqc_title=multiqc_title, genome=genome, fasta=fasta, bwt2_index=bwt2_index, digestion=digestion, restriction_site=restriction_site, ligation_site=ligation_site, chromosome_size=chromosome_size, restriction_fragments=restriction_fragments, save_reference=save_reference, dnase=dnase, min_cis_dist=min_cis_dist, split_fastq=split_fastq, fastq_chunks_size=fastq_chunks_size, min_mapq=min_mapq, bwt2_opts_end2end=bwt2_opts_end2end, bwt2_opts_trimmed=bwt2_opts_trimmed, save_aligned_intermediates=save_aligned_intermediates, keep_dups=keep_dups, keep_multi=keep_multi, max_insert_size=max_insert_size, min_insert_size=min_insert_size, max_restriction_fragment_size=max_restriction_fragment_size, min_restriction_fragment_size=min_restriction_fragment_size, save_interaction_bam=save_interaction_bam, save_pairs_intermediates=save_pairs_intermediates, bin_size=bin_size, hicpro_maps=hicpro_maps, ice_filter_low_count_perc=ice_filter_low_count_perc, ice_filter_high_count_perc=ice_filter_high_count_perc, ice_eps=ice_eps, ice_max_iter=ice_max_iter, res_zoomify=res_zoomify, save_raw_maps=save_raw_maps, res_dist_decay=res_dist_decay, tads_caller=tads_caller, res_tads=res_tads, res_compartments=res_compartments, skip_maps=skip_maps, skip_dist_decay=skip_dist_decay, skip_tads=skip_tads, skip_compartments=skip_compartments, skip_balancing=skip_balancing, skip_mcool=skip_mcool, skip_multiqc=skip_multiqc, multiqc_methods_description=multiqc_methods_description)

