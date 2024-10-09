
from dataclasses import dataclass
import typing
import typing_extensions
from enum import Enum


from flytekit.core.annotation import FlyteAnnotation

from latch.types.metadata import NextflowParameter
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchRule,
    NextflowParameter,
    Params,
    Section,
    Spoiler,
    Text,
)
# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters

# sample,fastq_1,fastq_2

flow = [
    Section(
        "Inputs",
        Params("input"),
    ),
    Section(
        "Reference genome options",
        Params("genome", "fasta", "bwt2_index", "save_reference")
    ),
    Section(
        "Output Directory",
        Params("outdir")
    ),
    Spoiler(
        "Digestion Hi-C",
        Params("digestion",
               "restriction_site",
               "ligation_site",
               "chromosome_size",
               "restriction_fragments")
    ),
    Spoiler(
        "DNase Hi-C",
        Params("dnase", "min_cis_dist")
    ),
    Spoiler(
        "Alignments",
        Params("split_fastq",
               "fastq_chunks_size",
               "min_mapq",
               "bwt2_opts_end2end",
               "bwt2_opts_trimmed",
               "save_aligned_intermediates")
    ),
    Spoiler(
        "Valid Pairs Detection",
        Params("keep_dups",
               "keep_multi",
               "max_insert_size",
               "min_insert_size",
               "max_restriction_fragment_size",
               "min_restriction_fragment_size",
               "save_interaction_bam",
               "save_pairs_intermediates")
    ),
    Spoiler(
        "Contact Maps",
        Params("bin_size",
               "hicpro_maps",
               "ice_filter_low_count_perc",
               "ice_filter_high_count_perc",
               "ice_eps",
               "ice_max_iter",
               "save_raw_maps")
    ),
    Spoiler(
        "Downstream Analysis",
        Params("tads_caller",
               "res_tads")
    ),
    Spoiler(
        "Skip Options",
        Params("skip_maps",
               "skip_dist_decay",
               "skip_tads",
               "skip_compartments",
               "skip_balancing",
               "skip_mcool",
               "skip_multiqc")
    ),
    Spoiler(
        "Generic Options",
        Params("email",
               "multiqc_title",
               "multiqc_methods_description")
    )
]

@dataclass(frozen=True)
class Sample:
    sample: str
    fastq_1: LatchFile
    fastq_2: LatchFile

class Digestion(Enum):
    hindiii = "hindiii"
    mboi = "mboi"
    dpnii = "dpnii"
    arima = "arima"

class Genome(Enum):
    mm10 = "mm10"
    GRCh37 = "GRCh37"

generated_parameters = {
    'input': NextflowParameter(
        display_name="input",
        type=typing.List[Sample],
        samplesheet=True,
        samplesheet_type='csv',
        default=None,
        section_title='Input/output options',
        description='Path to comma-separated file containing information about the samples in the experiment.',
    ),
    'outdir': NextflowParameter(
        display_name="outdir",
        type=typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})],
        default=None,
        section_title=None,
        description='The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.',
    ),
    'email': NextflowParameter(
        display_name="email",
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Email address for completion summary.',
    ),
    'multiqc_title': NextflowParameter(
        display_name="multiqc_title",
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='MultiQC report title. Printed as page header, used for filename if not otherwise specified.',
    ),
    'genome': NextflowParameter(
        display_name="genome",
        type=typing.Optional[Genome],
        default=None,
        section_title='Reference genome options',
        description='Name of iGenomes reference.',
    ),
    'fasta': NextflowParameter(
        display_name="fasta",
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Path to FASTA genome file.',
    ),
    'bwt2_index': NextflowParameter(
        display_name="bwt2_index",
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Full path to directory containing Bowtie index including base name. i.e. `/path/to/index/base`.',
    ),
    'digestion': NextflowParameter(
        display_name="digestion",
        type=typing.Optional[Digestion],
        default=None,
        section_title='Digestion Hi-C',
        description='Name of restriction enzyme to automatically set the restriction_site and ligation_site options (hindiii, mboi, dpnii, arima)',
    ),
    'restriction_site': NextflowParameter(
        display_name="restriction_site",
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Restriction motifs used during digestion. Several motifs (comma separated) can be provided.',
    ),
    'ligation_site': NextflowParameter(
        display_name="ligation_site",
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Expected motif after DNA ligation.  Several motifs (comma separated) can be provided.',
    ),
    'chromosome_size': NextflowParameter(
        display_name="chromosome_size",
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Full path to file specifying chromosome sizes (tab separated with chromosome name and size)`.',
    ),
    'restriction_fragments': NextflowParameter(
        display_name="restriction_fragments",
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Full path to restriction fragment (bed) file.',
    ),
    'save_reference': NextflowParameter(
        display_name="save_reference",
        type=bool,
        default=False,
        section_title=None,
        description='If generated by the pipeline save the annotation and indexes in the results directory.',
    ),
    'dnase': NextflowParameter(
        display_name="dnase",
        type=bool,
        default=False,
        section_title='DNAse Hi-C',
        description='For Hi-C protocols which are not based on enzyme digestion such as DNase Hi-C',
    ),
    'min_cis_dist': NextflowParameter(
        display_name="min_cis_dist",
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='Minimum distance between loci to consider. Useful for --dnase mode to remove spurious ligation products. Only values > 0 are considered',
    ),
    'split_fastq': NextflowParameter(
        display_name="split_fastq",
        type=bool,
        default=False,
        section_title='Alignments',
        description='Split the reads into chunks before running the pipelne',
    ),
    'fastq_chunks_size': NextflowParameter(
        display_name="fastq_chunks_size",
        type=typing.Optional[int],
        default=20000000,
        section_title=None,
        description='Read number per chunks if split_fastq is used',
    ),
    'min_mapq': NextflowParameter(
        display_name="min_mapq",
        type=typing.Optional[int],
        default=10,
        section_title=None,
        description='Keep aligned reads with a minimum quality value',
    ),
    'bwt2_opts_end2end': NextflowParameter(
        display_name="bwt2_opts_end2end",
        type=typing.Optional[str],
        default="'--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder'",
        section_title=None,
        description='Option for HiC-Pro end-to-end bowtie mapping',
    ),
    'bwt2_opts_trimmed': NextflowParameter(
        display_name="bwt2_opts_trimmed",
        type=typing.Optional[str],
        default="'--very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder'",
        section_title=None,
        description='Option for HiC-Pro trimmed reads mapping',
    ),
    'save_aligned_intermediates': NextflowParameter(
        display_name="save_aligned_intermediates",
        type=bool,
        default=False,
        section_title=None,
        description='Save all BAM files during two-steps mapping',
    ),
    'keep_dups': NextflowParameter(
        display_name="keep_dups",
        type=bool,
        default=False,
        section_title='Valid Pairs Detection',
        description='Keep duplicated reads',
    ),
    'keep_multi': NextflowParameter(
        display_name="keep_multi",
        type=bool,
        default=False,
        section_title=None,
        description='Keep multi-aligned reads',
    ),
    'max_insert_size': NextflowParameter(
        display_name="max_insert_size",
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='Maximum fragment size to consider. Only values > 0 are considered',
    ),
    'min_insert_size': NextflowParameter(
        display_name="min_insert_size",
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='Minimum fragment size to consider. Only values > 0 are considered',
    ),
    'max_restriction_fragment_size': NextflowParameter(
        display_name="max_restriction_fragment_size",
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='Maximum restriction fragment size to consider. Only values > 0 are considered',
    ),
    'min_restriction_fragment_size': NextflowParameter(
        display_name="min_restriction_fragment_size",
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='Minimum restriction fragment size to consider. Only values > 0 are considered',
    ),
    'save_interaction_bam': NextflowParameter(
        display_name="save_interaction_bam",
        type=bool,
        default=False,
        section_title=None,
        description='Save a BAM file where all reads are flagged by their interaction classes',
    ),
    'save_pairs_intermediates': NextflowParameter(
        display_name="save_pairs_intermediates",
        type=bool,
        default=False,
        section_title=None,
        description='Save all types of non valid read pairs in distinct output files',
    ),
    'bin_size': NextflowParameter(
        display_name="bin_size",
        type=typing.Optional[str],
        default='1000000,500000',
        section_title='Contact maps',
        description='Resolution to build the maps (comma separated)',
    ),
    'hicpro_maps': NextflowParameter(
        display_name="hicpro_maps",
        type=bool,
        default=False,
        section_title=None,
        description='Generate raw and normalized contact maps with HiC-Pro',
    ),
    'ice_filter_low_count_perc': NextflowParameter(
        display_name="ice_filter_low_count_perc",
        type=typing.Optional[float],
        default=0.02,
        section_title=None,
        description='Filter low counts rows before HiC-Pro normalization',
    ),
    'ice_filter_high_count_perc': NextflowParameter(
        display_name="ice_filter_high_count_perc",
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='Filter high counts rows before HiC-Pro normalization',
    ),
    'ice_eps': NextflowParameter(
        display_name="ice_eps",
        type=typing.Optional[float],
        default=0.1,
        section_title=None,
        description='Threshold for HiC-Pro ICE convergence',
    ),
    'ice_max_iter': NextflowParameter(
        display_name="ice_max_iter",
        type=typing.Optional[int],
        default=100,
        section_title=None,
        description='Maximum number of iteraction for HiC-Pro ICE normalization',
    ),
    # 'res_zoomify': NextflowParameter(
    #     display_name="res_zoomify",
    #     type=typing.Optional[str],
    #     default='5000',
    #     section_title=None,
    #     description='Maximum resolution to build mcool file',
    # ),
    'save_raw_maps': NextflowParameter(
        display_name="save_raw_maps",
        type=bool,
        default=False,
        section_title=None,
        description='Save raw contact maps',
    ),
    # 'res_dist_decay': NextflowParameter(
    #     display_name="res_dist_decay",
    #     type=typing.Optional[str],
    #     default='1000000',
    #     section_title='Downstream Analysis',
    #     description='Resolution to build count/distance plot',
    # ),
    'tads_caller': NextflowParameter(
        display_name="tads_caller",
        type=typing.Optional[str],
        default='hicexplorer,insulation',
        section_title=None,
        description='Define methods for TADs calling',
    ),
    'res_tads': NextflowParameter(
        display_name="res_tads",
        type=typing.Optional[str],
        default='40000,20000',
        section_title=None,
        description='Resolution to run TADs callers (comma separated)',
    ),
    # 'res_compartments': NextflowParameter(
    #     display_name="res_compartments",
    #     type=typing.Optional[str],
    #     default='250000',
    #     section_title=None,
    #     description='Resolution for compartments calling',
    # ),
    'skip_maps': NextflowParameter(
        display_name="skip_maps",
        type=bool,
        default=False,
        section_title='Skip options',
        description='Do not build contact maps',
    ),
    'skip_dist_decay': NextflowParameter(
        display_name="skip_dist_decay",
        type=bool,
        default=False,
        section_title=None,
        description='Do not run distance/decay plot',
    ),
    'skip_tads': NextflowParameter(
        display_name="skip_tads",
        type=bool,
        default=False,
        section_title=None,
        description='Do not run TADs calling',
    ),
    'skip_compartments': NextflowParameter(
        display_name="skip_compartments",
        type=bool,
        default=False,
        section_title=None,
        description='Do not run compartments calling',
    ),
    'skip_balancing': NextflowParameter(
        display_name="skip_balancing",
        type=bool,
        default=False,
        section_title=None,
        description='Do not run cooler balancing normalization',
    ),
    'skip_mcool': NextflowParameter(
        display_name="skip_mcool",
        type=bool,
        default=False,
        section_title=None,
        description='Do not generate mcool file for Higlass visualization',
    ),
    'skip_multiqc': NextflowParameter(
        display_name="skip_multiqc",
        type=bool,
        default=False,
        section_title=None,
        description='Do not generate MultiQC report',
    ),
    'multiqc_methods_description': NextflowParameter(
        display_name="multiqc_methods_description",
        type=typing.Optional[str],
        default=None,
        section_title='Generic options',
        description='Custom MultiQC yaml file containing HTML including a methods description.',
    ),
}
