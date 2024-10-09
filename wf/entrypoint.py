import sys
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
from latch.executions import report_nextflow_used_storage
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
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage-ofs",
        headers=headers,
        json={
            "storage_expiration_hours": 0,
            "version": 2,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


class Digestion(Enum):
    hindiii = 'hindiii'
    mboi = 'mboi'
    dpnii = 'dpnii'
    arima = 'arima'


@dataclass
class Sample:
    sample: str
    fastq_1: LatchFile
    fastq_2: LatchFile

class Genome(Enum):
    mm10 = "mm10"
    GRCh37 = "GRCh37"


input_construct_samplesheet = metadata._nextflow_metadata.parameters['input'].samplesheet_constructor


@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(pvc_name: str, input: typing.List[Sample], outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], email: typing.Optional[str], multiqc_title: typing.Optional[str], genome: typing.Optional[Genome], fasta: typing.Optional[LatchFile], bwt2_index: typing.Optional[str], digestion: typing.Optional[Digestion], restriction_site: typing.Optional[str], ligation_site: typing.Optional[str], chromosome_size: typing.Optional[LatchFile], restriction_fragments: typing.Optional[LatchFile], min_cis_dist: typing.Optional[int], max_insert_size: typing.Optional[int], min_insert_size: typing.Optional[int], max_restriction_fragment_size: typing.Optional[int], min_restriction_fragment_size: typing.Optional[int], ice_filter_high_count_perc: typing.Optional[int], multiqc_methods_description: typing.Optional[str], save_reference: bool, dnase: bool, split_fastq: bool, fastq_chunks_size: typing.Optional[int], min_mapq: typing.Optional[int], bwt2_opts_end2end: typing.Optional[str], bwt2_opts_trimmed: typing.Optional[str], save_aligned_intermediates: bool, keep_dups: bool, keep_multi: bool, save_interaction_bam: bool, save_pairs_intermediates: bool, bin_size: typing.Optional[str], hicpro_maps: bool, ice_filter_low_count_perc: typing.Optional[float], ice_eps: typing.Optional[float], ice_max_iter: typing.Optional[int], save_raw_maps: bool, tads_caller: typing.Optional[str], res_tads: typing.Optional[str], skip_maps: bool, skip_dist_decay: bool, skip_tads: bool, skip_compartments: bool, skip_balancing: bool, skip_mcool: bool, skip_multiqc: bool) -> None:
    shared_dir = Path("/nf-workdir")

    exec_name = _get_execution_name()
    if exec_name is None:
        print("Failed to get execution name.")
        exec_name = "unknown"

    latch_log_dir = urljoins("latch:///nf_core_hic_logs/nf_nf_core_hic", exec_name)
    print(f"Log directory: {latch_log_dir}")


    input_samplesheet = input_construct_samplesheet(input)

    ignore_list = [
        "latch",
        ".latch",
        ".git",
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

    profile_list = []
    if False:
        profile_list.extend([p.value for p in execution_profiles])

    if len(profile_list) == 0:
        profile_list.append("standard")

    profiles = ','.join(profile_list)

    cmd = [
        "/root/nextflow",
        "run",
        str(shared_dir / "main.nf"),
        "-work-dir",
        str(shared_dir),
        "-profile",
        profiles,
        "-c",
        "latch.config",
        "-resume",
        *get_flag('input', input_samplesheet),
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
                *get_flag('save_raw_maps', save_raw_maps),
                *get_flag('tads_caller', tads_caller),
                *get_flag('res_tads', res_tads),
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

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
        }

        if False:
            env["LATCH_LOG_DIR"] = latch_log_dir

        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            remote = LPath(urljoins(latch_log_dir, "nextflow.log"))
            print(f"Uploading .nextflow.log to {remote.path}")
            remote.upload_from(nextflow_log)

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ['du', '-sb', str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print("Failed to compute storage size: Operation timed out after 5 minutes.")
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)


@workflow(metadata._nextflow_metadata)
def nf_nf_core_hic(input: typing.List[Sample], outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], email: typing.Optional[str], multiqc_title: typing.Optional[str], genome: typing.Optional[Genome], fasta: typing.Optional[LatchFile], bwt2_index: typing.Optional[str], digestion: typing.Optional[Digestion], restriction_site: typing.Optional[str], ligation_site: typing.Optional[str], chromosome_size: typing.Optional[LatchFile], restriction_fragments: typing.Optional[LatchFile], min_cis_dist: typing.Optional[int], max_insert_size: typing.Optional[int], min_insert_size: typing.Optional[int], max_restriction_fragment_size: typing.Optional[int], min_restriction_fragment_size: typing.Optional[int], ice_filter_high_count_perc: typing.Optional[int], multiqc_methods_description: typing.Optional[str], save_reference: bool = False, dnase: bool = False, split_fastq: bool = False, fastq_chunks_size: typing.Optional[int] = 20000000, min_mapq: typing.Optional[int] = 10, bwt2_opts_end2end: typing.Optional[str] = "'--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder'", bwt2_opts_trimmed: typing.Optional[str] = "'--very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder'", save_aligned_intermediates: bool = False, keep_dups: bool = False, keep_multi: bool = False, save_interaction_bam: bool = False, save_pairs_intermediates: bool = False, bin_size: typing.Optional[str] = '1000000,500000', hicpro_maps: bool = False, ice_filter_low_count_perc: typing.Optional[float] = 0.02, ice_eps: typing.Optional[float] = 0.1, ice_max_iter: typing.Optional[int] = 100, save_raw_maps: bool = False, tads_caller: typing.Optional[str] = 'hicexplorer,insulation', res_tads: typing.Optional[str] = '40000,20000', skip_maps: bool = False, skip_dist_decay: bool = False, skip_tads: bool = False, skip_compartments: bool = False, skip_balancing: bool = False, skip_mcool: bool = False, skip_multiqc: bool = False) -> None:
    """

    nf-core/hic is a bioinformatics best-practice analysis pipeline for Analysis of Chromosome Conformation Capture data (Hi-C).


    <html>
    <p align="center">
    <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>
    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>
    [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3476425-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.3476425)

    # nf-core/hic

    [![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/hic/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.2669512-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.2669512)

    [![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)


    ## Introduction

    **nf-core/hic** is a bioinformatics best-practice analysis pipeline for Analysis of Chromosome Conformation Capture data (Hi-C).

    The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

    On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/hic/results).

    ## Pipeline summary

    1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    2. Hi-C data processing
    1. [`HiC-Pro`](https://github.com/nservant/HiC-Pro)
        1. Mapping using a two steps strategy to rescue reads spanning the ligation
            sites ([`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
        2. Detection of valid interaction products
        3. Duplicates removal
        4. Generate raw and normalized contact maps ([`iced`](https://github.com/hiclib/iced))
    3. Create genome-wide contact maps at various resolutions ([`cooler`](https://github.com/open2c/cooler))
    4. Contact maps normalization using balancing algorithm ([`cooler`](https://github.com/open2c/cooler))
    5. Export to various contact maps formats ([`HiC-Pro`](https://github.com/nservant/HiC-Pro), [`cooler`](https://github.com/open2c/cooler))
    6. Quality controls ([`HiC-Pro`](https://github.com/nservant/HiC-Pro), [`HiCExplorer`](https://github.com/deeptools/HiCExplorer))
    7. Compartments calling ([`cooltools`](https://cooltools.readthedocs.io/en/latest/))
    8. TADs calling ([`HiCExplorer`](https://github.com/deeptools/HiCExplorer), [`cooltools`](https://cooltools.readthedocs.io/en/latest/))
    9. Quality control report ([`MultiQC`](https://multiqc.info/))

    ## Pipeline output

    To see the the results of a test run with a full size dataset refer to the [results](https://nf-co.re/hic/results) tab on the nf-core website pipeline page.
    For more details about the output files and reports, please refer to the
    [output documentation](https://nf-co.re/hic/output).

    ## Credits

    nf-core/hic was originally written by Nicolas Servant.

    ## Contributions and Support

    If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

    For further information or help, don't hesitate to get in touch on the [Slack `#hic` channel](https://nfcore.slack.com/channels/hic) (you can join with [this invite](https://nf-co.re/join/slack)).

    ## Citations

    If you use nf-core/hic for your analysis, please cite it using the following doi: doi: [10.5281/zenodo.2669512](https://doi.org/10.5281/zenodo.2669512)

    An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

    You can cite the `nf-core` publication as follows:

    > **The nf-core framework for community-curated bioinformatics pipelines.**
    >
    > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
    >
    > _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

    """

    pvc_name: str = initialize()
    nextflow_runtime(pvc_name=pvc_name, input=input, outdir=outdir, email=email, multiqc_title=multiqc_title, genome=genome, fasta=fasta, bwt2_index=bwt2_index, digestion=digestion, restriction_site=restriction_site, ligation_site=ligation_site, chromosome_size=chromosome_size, restriction_fragments=restriction_fragments, save_reference=save_reference, dnase=dnase, min_cis_dist=min_cis_dist, split_fastq=split_fastq, fastq_chunks_size=fastq_chunks_size, min_mapq=min_mapq, bwt2_opts_end2end=bwt2_opts_end2end, bwt2_opts_trimmed=bwt2_opts_trimmed, save_aligned_intermediates=save_aligned_intermediates, keep_dups=keep_dups, keep_multi=keep_multi, max_insert_size=max_insert_size, min_insert_size=min_insert_size, max_restriction_fragment_size=max_restriction_fragment_size, min_restriction_fragment_size=min_restriction_fragment_size, save_interaction_bam=save_interaction_bam, save_pairs_intermediates=save_pairs_intermediates, bin_size=bin_size, hicpro_maps=hicpro_maps, ice_filter_low_count_perc=ice_filter_low_count_perc, ice_filter_high_count_perc=ice_filter_high_count_perc, ice_eps=ice_eps, ice_max_iter=ice_max_iter, save_raw_maps=save_raw_maps, tads_caller=tads_caller, res_tads=res_tads, skip_maps=skip_maps, skip_dist_decay=skip_dist_decay, skip_tads=skip_tads, skip_compartments=skip_compartments, skip_balancing=skip_balancing, skip_mcool=skip_mcool, skip_multiqc=skip_multiqc, multiqc_methods_description=multiqc_methods_description)

