# CHANGELOG


## v1.17.0 (2026-02-11)

### Features

- Implement stubs for cli tests to streamline testing process
  ([`2bff724`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2bff724))

- Refactor alignment classes to streamline output handling and remove redundant fields
  ([`563c9b9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/563c9b9))

- Pixi lock
  ([`7a30ad0`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7a30ad0))

- Update alignment classes to use cram format and streamline commands
  ([`9bfdae3`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/9bfdae3))

- Rename bamfilter classes to samtoolsfilter and update output to cram format
  ([`d9db599`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d9db599))

- Add output directory option for tree command
  ([`5bcac09`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/5bcac09))

- Add samclip command and utility for filtering clipped reads from sam files
  ([`bd59abb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/bd59abb))

- Add aligner options for long and short read pipelines and adjust histogram height
  ([`75e815d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/75e815d))

- Pass bam through without copying
  ([`745e326`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/745e326))

- Allow global param override in order
  ([`8e8e622`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8e8e622))

- Better pipeline error handling
  ([`2252067`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2252067))

- Enhance aligner class and add test for zero mapped reads
  ([`fb4a1b1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fb4a1b1))

- Enhance shellcommand and basestage to support output files in commands and pipelines
  ([`a44a4af`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a44a4af))

- Add test for paf output file creation and validation
  ([`94f4744`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/94f4744))

- Update multi cli to support gathered config
  ([`cf8f615`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/cf8f615))

- Add snippy-ng gather command
  ([`5aaac45`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/5aaac45))

- Add runtime testing and test discovery to stages
  ([`fe15bf3`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fe15bf3))


### Bug Fixes

- Clarify error message for unmapped reads in bam file after filtering
  ([`10a85d4`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/10a85d4))

- Update samtools version constraint and drop samclip
  ([`073afef`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/073afef))

- Move mapped reads test to filtering stage
  ([`59242d6`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/59242d6))

- Improve error message formatting in stagetestfailure
  ([`bed8b74`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/bed8b74))

- Exclude samples.csv from version control
  ([`c33fecb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c33fecb))

- Update example commands to use genbank reference files
  ([`9ec0fba`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/9ec0fba))

- Update error messages for missing dependencies
  ([`994057a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/994057a))

- Remove unused 'phylip' output path from iqtreebuildtreeoutput
  ([`6be713d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6be713d))

- Add output validation and testing after each stage run
  ([`0068660`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0068660))

- Improve error message for missing dependencies
  ([`08cdaf6`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/08cdaf6))

- Add 'long' to .gitignore to exclude long output files
  ([`2eba508`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2eba508))

- Add caller and clair3_model to handle_ont entry
  ([`103121d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/103121d))

- Remove show_default from cpus-per-sample option in multi command
  ([`43febe6`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/43febe6))

- Read fconst from file if it is a path in tree command
  ([`2a61d1f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2a61d1f))

- Ensure terminal width is at least 40 columns for histogram display
  ([`cd8cc11`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/cd8cc11))


### Refactoring

- Update bam parameter type to path and streamline alignment logic
  ([`add399b`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/add399b))

- Streamline command construction and error handling across multiple stages
  ([`3745270`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3745270))

- Simplify exception handling in run_snippy_pipeline
  ([`4574d0e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4574d0e))

- Unify exception hierarchy under snippyerror base class
  ([`784fb90`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/784fb90))

- Multi-sample pipeline
  ([`a8c7021`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a8c7021))

- Update import paths for pipeline_runner
  ([`f540e0e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f540e0e))


### Chores

- Rename fastq reads
  ([`73aae93`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/73aae93))


## v1.16.0 (2026-02-02)

### Features

- Set clair3 as the default caller
  ([`d846de0`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d846de0))

- Add printvcfhistogram stage to visualize vcf data in terminal
  ([`4adb67b`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4adb67b))

- Add minimap_preset option for customizable minimap2 alignment
  ([`3dfc513`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3dfc513))


### Bug Fixes

- Clean up completion output in snippy class
  ([`b134bb7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b134bb7))

- Improve logging in dependency validation and welcome message
  ([`f9ed912`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f9ed912))

- Optimize debug mode check in logger class
  ([`2f3c78e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2f3c78e))

- Remove unused import of referencemetadata in base.py
  ([`c0c42e8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c0c42e8))

- Update seqkitcleanlongreads documentation and adjust min_length description
  ([`429fd98`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/429fd98))

- Add --remove-gaps option to seqkit command for long read cleaning
  ([`2f4a57f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2f4a57f))

- Update vcffilter commands to use binary
  ([`56ffab1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/56ffab1))

- Replace continue_last_run with create_missing
  ([`5d79604`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/5d79604))

- Add platform detection for clair3caller based on model type
  ([`5e60bb3`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/5e60bb3))

- Add fasta_file_matches_ref method for validating fasta file against reference metadata
  ([`d2a305c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d2a305c))

- Rename metadata to ref_metadata for consistency in rasusadownsamplereads tests
  ([`d6345ef`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d6345ef))

- Remove ref metadata from base stage
  ([`3034561`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3034561))

- Rename metadata to ref_metadata for consistency across pipeline stages
  ([`50e4f07`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/50e4f07))

- Update cleanup method in test stubs to accept directory parameter fix: update test for guess_format to raise error for non-existent files
  ([`a3ed236`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a3ed236))

- Replace valueerror with msavalidationerror for better error handling in combinefastafile
  ([`7c51d81`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7c51d81))

- Update multi command to improve cpu allocation and reference loading
  ([`1260239`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1260239))

- Refactor guess_format function for improved readability and structure
  ([`dac4cd2`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/dac4cd2))

- Update cleanup method to accept directory parameter for resetting working directory
  ([`c2fde82`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c2fde82))

- Pass current directory to snippy.cleanup for proper cleanup
  ([`b3e0f79`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b3e0f79))

- Update load_or_prepare_reference to accept output_directory parameter
  ([`f7e8ce4`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f7e8ce4))

- Rename variable for clarity in loadreferencefrommetadatafile output method
  ([`c4094fa`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c4094fa))

- Ensure command order in bugcatchinggroup
  ([`0b2a985`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0b2a985))

- Add usage example to multi command docstring for clarity
  ([`0852988`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0852988))


### Documentation

- Add usage examples for short and long commands in readme
  ([`413cc4e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/413cc4e))


## v1.15.0 (2026-01-27)

### Features

- Add core-snp-filter and iqtree packages to bioconda dependencies
  ([`fa602c3`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fa602c3))

- Add multi command for multi-sample processing and sample configuration yaml
  ([`b6f66f6`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b6f66f6))

- Add clair3_fast_mode option to long read pipeline configuration
  ([`a8c6ba4`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a8c6ba4))

- Add core-snp-filter and iqtree dependencies to bioconda
  ([`5cf573e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/5cf573e))

- Add tree command for phylogenetic tree construction and integrate iq-tree stages
  ([`0c232de`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0c232de))

- Enhance alignment processing with constant sites extraction and improve error handling
  ([`70ce77e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/70ce77e))

- Add aln command to cli for alignment processing
  ([`147ab4f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/147ab4f))

- Implement alignment creation and pipeline stages for multi-sample processing
  ([`fbe1dee`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fbe1dee))

- Add not implemented callback for new options and enhance global option definitions
  ([`6637af9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6637af9))

- Add reference preparation command to cli for genome processing
  ([`8f29435`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8f29435))

- Add option to exclude insertions in bcftools pseudo-alignment
  ([`8f4be5a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8f4be5a))

- Add command attribute to dependency class for customizable command execution
  ([`c721de5`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c721de5))

- Enhance vcffilterlong to manage vcf tags and annotations
  ([`2b36db6`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2b36db6))

- Add alias for --assembly option in asm command
  ([`50089c5`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/50089c5))


### Bug Fixes

- Reorder imports in multi_cli.py for better organization and readability
  ([`8b922e9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8b922e9))

- Remove bwa package entries from pixi.lock to clean up dependencies
  ([`43697bf`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/43697bf))

- Add 'multi' directory to .gitignore to prevent unnecessary tracking
  ([`2525f41`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2525f41))

- Update outdir default value to path(".") for consistency and improve error handling with stageexecutionerror
  ([`4f4c844`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4f4c844))

- Correct typo in development instructions and add example for multi command usage
  ([`740cd49`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/740cd49))

- Rename 'fbopt' argument to 'freebayes_opts' for clarity and consistency
  ([`af3a043`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/af3a043))

- Change outdir type from str to path for better path handling
  ([`07e014c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/07e014c))

- Remove bwa dependency from bioconda in pyproject.toml
  ([`6bae243`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6bae243))

- Rename argument 'aln' to 'alignment' for consistency in tree command
  ([`1f73f4e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1f73f4e))

- Improve error handling for incomplete outputs in snippy class
  ([`bcf6349`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/bcf6349))

- Remove unused import of path in downsample_reads.py
  ([`cb6713e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/cb6713e))

- Return a list of commands in copyfasta class for consistency
  ([`88d0ac7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/88d0ac7))

- Streamline command handling in copyfasta class for better readability
  ([`c335ffe`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c335ffe))

- Add is_eager option to global_defs for immediate flag evaluation
  ([`5e8010e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/5e8010e))


### Documentation

- Add example usage for run_clair3_docker.sh
  ([`83026af`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/83026af))

- Add example usage for snippy-ng command
  ([`1c6b728`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1c6b728))


### Refactoring

- Replace basemodel with baseoutput
  ([`7ef6bdb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7ef6bdb))

- Replace snippy_global_options with add_snippy_global_options for consistency across cli commands
  ([`481620d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/481620d))

- Update bug report template messages for clarity and consistency
  ([`b79768f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b79768f))

- Streamline reference preparation logic and improve error handling for reference format detection
  ([`1c1d233`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1c1d233))

- Update reference option help text to include prepared reference directory
  ([`b798a0c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b798a0c))

- Update load_or_prepare_reference to use loadreferencefrommetadatafile and improve error handling for metadata.json
  ([`a31cec4`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a31cec4))

- Enhance referencemetadata class with improved error handling and data validation
  ([`913ee06`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/913ee06))

- Add stageexecutionerror exception and update error handling in snippy pipeline
  ([`180d962`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/180d962))

- Update run_snippy_pipeline to accept individual parameters instead of a config dictionary
  ([`0bed04a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0bed04a))

- Standardize parameter name in metadata class constructor
  ([`8b611ef`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8b611ef))

- Update metadata handling to use path instead of metadata class
  ([`d3ace9d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d3ace9d))

- Replace at_run_time with metadata class
  ([`50e7f54`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/50e7f54))

- Simplify bcftoolsconsequencescaller class
  ([`a3fd5eb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a3fd5eb))


## v1.14.0 (2025-12-15)

### Features

- Pass clair3_fast_mode to clair3caller for enhanced variant calling
  ([`51e8b8d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/51e8b8d))

- Add --prefix option to globals
  ([`6c88b86`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6c88b86))

- Add --clair3-fast-mode option for quicker variant calling in clair3
  ([`89a8bd6`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/89a8bd6))

- Add clair3 dependency
  ([`c125ae5`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c125ae5))

- Add run_clair3_docker.sh script to run clair3 in docker
  ([`f60ad7a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f60ad7a))

- Add long command to snippy-ng cli
  ([`19fe72e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/19fe72e))

- Implement long read snp calling pipeline with configurable options
  ([`ca73c0b`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/ca73c0b))

- Add freebayescallerlong class for long-read variant calling with parallel processing
  ([`134f55c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/134f55c))

- Update depthmask and hetmask classes to use fasta instead of reference and improve min-depth masking logic
  ([`341180f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/341180f))

- Add vcffilterlong class for long-read variant filtering using bcftools
  ([`3242ab9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3242ab9))

- Add seqkit stage for cleaning long reads and update dependencies
  ([`fb3bc6e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fb3bc6e))

- Update basestage model configuration and add prefix field
  ([`dc9d01f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/dc9d01f))

- Add test modules for rasusa and seqkit read downsampling and statistics stages
  ([`8e624ef`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8e624ef))

- Refactor pipeline stage creation functions
  ([`bad38cf`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/bad38cf))

- Add additional sequences to wildtype.fasta for comprehensive testing
  ([`37622d8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/37622d8))

- Enhance pipeline stages with detailed configuration and improve error handling
  ([`55ba4d7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/55ba4d7))

- Implement create_short_stages function and refactor short cli command
  ([`e2f0680`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e2f0680))


### Bug Fixes

- Remove aligner_opts parameter from asm function
  ([`a1b9068`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a1b9068))

- Remove unused aligner_opts parameter from create_asm_pipeline_stages
  ([`5a7f969`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/5a7f969))

- Remove min_frac from vcffilterlong to simplify variant filtering
  ([`d5f1370`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d5f1370))

- Remove unused prefix option from various stages
  ([`74bbdc8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/74bbdc8))

- Remove unused header option and outdir from asm and short pipeline stages
  ([`d78ed57`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d78ed57))

- Remove unused header option from asm and short commands
  ([`f2512b0`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f2512b0))

- Improve stage logging format for better readability
  ([`0f58de0`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0f58de0))

- Remove unused prefix field from pseudoalignment class
  ([`77750eb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/77750eb))

- Adjust logging output format for shell command execution
  ([`79398f2`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/79398f2))


### Refactoring

- Move clean reads stage to maintain current reads flow
  ([`4fe1c3a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4fe1c3a))

- Run_snippy_pipeline function
  ([`dc58acb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/dc58acb))

- Rename shellcommandpipeline to shellcommandpipe for consistency
  ([`7ab673c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7ab673c))


### Chores

- Add unit tests for dependency and snippy classes
  ([`8115a4d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8115a4d))

- Add parameterized tests for asm, long, and short commands
  ([`d9d4092`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d9d4092))


## v1.13.0 (2025-11-18)

### Features

- Implement atruntime class for deferred execution and update related tests
  ([`583a9c1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/583a9c1))

- Global debug option
  ([`7bdb4b5`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7bdb4b5))


### Bug Fixes

- Import path from pathlib in runtime.py
  ([`784ea7a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/784ea7a))

- Update environment variable name for debug logging
  ([`40f9552`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/40f9552))

- Remove min-depth and min-qual options from asm cli and hardcode values in vcffilter
  ([`b88f046`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b88f046))

- Remove skipstageerror exception handling and related class
  ([`2a6806a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2a6806a))


### Refactoring

- Replace kwargs with config in asm function for clarity
  ([`a0bac67`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a0bac67))

- Update imports for at_run_time
  ([`efa4401`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/efa4401))

- Update parameter names for clarity in asm command and global options
  ([`424ade4`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/424ade4))

- Mark og command as deprecated and update help documentation for clarity
  ([`3cc10ed`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3cc10ed))

- Move atruntime class and at_run_time function to runtime.py for better organization
  ([`2d2a06f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2d2a06f))

- Move load_or_prepare_reference to common utilities and update imports in asm and short cli
  ([`e5bc9c2`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e5bc9c2))

- Replace pipeline with snippy class in cli modules and tests
  ([`75bc80d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/75bc80d))

- Streamline output directory handling in asm and short cli functions
  ([`dbe8317`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/dbe8317))

- Remove debug environment variable handling and unused exception
  ([`3f275c1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3f275c1))


## v1.12.2 (2025-10-17)

### Bug Fixes

- Update temporary fasta file suffix to .tmp in applymask class
  ([`c287288`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c287288))

- Update output file extension from .afa to .raw.fna in bcftoolspseudoalignment
  ([`eb88125`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/eb88125))


## v1.12.1 (2025-10-16)

### Bug Fixes

- Remove deprecated --continue option from short command
  ([`b3b0e40`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b3b0e40))


## v1.12.0 (2025-10-16)

### Features

- Import assemblyaligner and pafcaller in asm function for alignment and calling stages
  ([`e3407cb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e3407cb))

- Add asm_out to .gitignore to exclude assembly output files
  ([`5697bc4`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/5697bc4))

- Add keep_incomplete parameter to run method in test stub
  ([`b223db8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b223db8))

- Enhance reference handling by adding reference dictionary output and updating output structure
  ([`b89dc5f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b89dc5f))

- Refactor user masking stage to apply bed mask and manage temporary files
  ([`c8d78b7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c8d78b7))

- Update vcffilter parameters for improved variant calling criteria
  ([`63a3fd4`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/63a3fd4))

- Refactor freebayescaller output structure and add pafcaller for paf alignments
  ([`d40e4e5`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d40e4e5))

- Add assemblyaligner for aligning assemblies to references using minimap2
  ([`fff91e2`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fff91e2))

- Refactor short read snp calling pipeline to streamline reference loading and update output file naming
  ([`3b4733a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3b4733a))

- Implement reference loading and preparation utility in cli
  ([`98d623e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/98d623e))

- Add asm command for assembly-based snp calling pipeline
  ([`b4bdfb6`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b4bdfb6))

- Add version pattern for paftools.js dependency
  ([`aae7063`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/aae7063))

- Add options to continue from last run and keep incomplete outputs in pipeline
  ([`c4d147f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c4d147f))

- Update tmpdir field in basestage to use a default temporary directory
  ([`a9eea13`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a9eea13))

- Enhance stub_everything to support additional output parameters and modify test_run_cli argument
  ([`a30fb1a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a30fb1a))

- Reorder kwargs assignment and conditionally append seqkit read statistics stage
  ([`7492d96`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7492d96))

- Refactor genome length retrieval to use a closure for runtime metadata access
  ([`31c49c9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/31c49c9))

- Implement bam and vcf filtering stages with new bamfilter and vcffilter classes
  ([`67785b4`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/67785b4))

- Refactor command registration and optimize imports in cli
  ([`8249159`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8249159))

- Add paftools.js dependency
  ([`3ab9ccf`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3ab9ccf))

- Change default mask character from 'n' to 'x' in usermask class
  ([`38b0707`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/38b0707))

- Update minimapaligner to include samclip in dependencies and streamline alignment pipeline
  ([`cd4d732`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/cd4d732))


### Bug Fixes

- Update benchmark command to use 'asm' instead of 'og'
  ([`2c95c59`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2c95c59))

- Update min-qual threshold for variant calling to improve accuracy
  ([`eef6e7d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/eef6e7d))

- Improve error logging for failed commands in basestage
  ([`fadf387`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fadf387))


## v1.11.0 (2025-10-13)

### Features

- Enhance freebayes variant calling commands with improved filtering and normalization
  ([`e5a75f6`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e5a75f6))

- Add '--mark-del' option to bcftools consensus command for deletion marking
  ([`342df9f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/342df9f))


## v1.10.0 (2025-10-13)

### Features

- Add og_out to .gitignore to exclude output files from version control
  ([`1439068`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1439068))

- Enhance masking functionality by adding depth and heterozygous site masking options
  ([`ace21a6`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/ace21a6))

- Remove mask option from pseudoalignment and related command adjustments
  ([`59d79e0`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/59d79e0))

- Improve code formatting and enhance clarity in freebayescaller class
  ([`0c83a8a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0c83a8a))

- Refactor masking stages to enhance depth and user masking functionality
  ([`65cdca1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/65cdca1))

- Add check_outputs method to verify existence of expected output files
  ([`3e13814`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3e13814))

- Add continue option to skip completed stages in pipeline execution
  ([`0125763`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0125763))


## v1.9.0 (2025-10-13)

### Features

- Implement comprehensive masking and update output file handling in short read snp pipeline
  ([`0983687`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0983687))

- Implement comprehensive masking stage with depth-based mask generation and application
  ([`cd2bc65`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/cd2bc65))

- Enhance copyfile and copyfasta classes for improved file copying functionality
  ([`5d97d3c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/5d97d3c))


### Bug Fixes

- Update output file extension from .fasta to .afa in bcftoolspseudoalignment
  ([`ea560b9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/ea560b9))


## v1.8.0 (2025-10-10)

### Features

- Add more whitespace
  ([`e3b72cd`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e3b72cd))

- Add copyfile stage to rename final consensus output
  ([`71229dd`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/71229dd))

- Add copyfile stage to handle file copying
  ([`7cc901b`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7cc901b))

- Enhance short cli with masking options and integrate bcftools pseudo-alignment
  ([`8df0137`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8df0137))

- Add bedtools dependency with version constraints
  ([`d77a2a7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d77a2a7))

- Add bedtools dependency with citation and version constraints
  ([`63cb361`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/63cb361))

- Implement bcftools pseudo-alignment stage with command construction
  ([`e280c73`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e280c73))

- Add bgzipcompressor class for gzip file compression
  ([`e7882d5`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e7882d5))

- Implement depthmask and applymask stages for bed file creation and masking
  ([`006bd36`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/006bd36))

- Add option to exclude insertions from freebayes variant calls
  ([`65a17ca`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/65a17ca))


### Bug Fixes

- Correct command assertion in test_seqkit_read_stats
  ([`4d9e7dd`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4d9e7dd))

- Ensure shell command arguments are properly quoted in string representation
  ([`d863f8b`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d863f8b))


### Refactoring

- Remove unused imports in alignment_filtering.py
  ([`ad3c763`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/ad3c763))


## v1.7.0 (2025-10-08)

### Features

- Add mkdocs_autorefs dependency for improved documentation support
  ([`450c34f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/450c34f))


### Bug Fixes

- Update installation script url and improve development setup instructions
  ([`6ed1e69`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6ed1e69))


### Refactoring

- Remove outdated test scripts for at_run_time and fastpcleanreadsoutput
  ([`bd6d846`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/bd6d846))

- Import loadreference in short command for improved pipeline setup
  ([`84f6a81`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/84f6a81))

- Update help text for --quiet option to clarify functionality
  ([`53092f9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/53092f9))

- Rename shellpipeline to shellcommandpipeline for consistency
  ([`bbb1278`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/bbb1278))

- Improve logging in run method to display command descriptions and format
  ([`477d899`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/477d899))

- Change logging of stages from info to debug level for improved verbosity
  ([`1e0ac64`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1e0ac64))


## v1.6.0 (2025-10-08)

### Features

- Refactor command construction to use shellcommand for improved safety and consistency
  ([`35e4d66`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/35e4d66))

- Add reference_index to kwargs in short command for improved pipeline setup
  ([`224bfa2`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/224bfa2))

- Add invalidcommandtypeerror exception for better error handling
  ([`fe7a2ec`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fe7a2ec))

- Add tstrings-backport dependency to project
  ([`2a2d2fa`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2a2d2fa))


### Bug Fixes

- Convert command outputs to string for consistent handling in tests
  ([`806cf1d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/806cf1d))

- Remove tstrings-backport dependency from project requirements
  ([`f641b24`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f641b24))

- Improve output validation in pipeline execution
  ([`970b004`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/970b004))


### Refactoring

- Enhance command construction and pipeline integration across all stages
  ([`7f0d864`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7f0d864))

- Enhance command handling and add shell pipeline support
  ([`7d74bf5`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7d74bf5))

- Remove test for nonexistent read file and update compression level flag
  ([`08b6e84`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/08b6e84))


## v1.5.0 (2025-10-03)

### Features

- Add bcftoolsconsequencescaller for consequence calling in pipeline
  ([`81fe6f8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/81fe6f8))


### Bug Fixes

- Update reference genome help text to clarify accepted formats
  ([`2003905`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2003905))

- Enhance gff3 output generation in process_reference method
  ([`1ce4319`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1ce4319))

- Remove bcbio-gff from runtime requirements in recipe.yaml
  ([`bdba584`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/bdba584))

- Remove bcbio-gff dependency from project and pixi configuration
  ([`7d33066`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7d33066))

- Improve error message for missing output files in pipeline execution
  ([`f937082`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f937082))

- Add support for python 3.13 in classifiers
  ([`c0434c0`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c0434c0))


## v1.4.0 (2025-10-02)

### Features

- Add development setup instructions to readme
  ([`3b5bdfd`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3b5bdfd))

- Refactor imports and add bug catching utility class
  ([`6bf6f82`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6bf6f82))


### Bug Fixes

- Remove vt dependency and replace with bcftools norm
  ([`0e8532d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0e8532d))

- Update example command in short function documentation
  ([`1ad69ef`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1ad69ef))

- Correct typo in filter vcf filename in freebayescaller output
  ([`dddd83e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/dddd83e))


## v1.3.0 (2025-08-06)

### Features

- Update benchmark tests to include 'og' and 'short' commands
  ([`db21824`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/db21824))

- Remove equality and inequality checks from atruntime class
  ([`86dfdeb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/86dfdeb))

- Update output_prefix to prefix in seqkitreadstats and related tests
  ([`e2abdee`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e2abdee))

- Add tests for rasusa downsampling stages and at_run_time functionality
  ([`b337ccb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b337ccb))

- Add cspell configuration for custom words in settings
  ([`2bf1980`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2bf1980))

- Add downsampling option to short cli for coverage specification
  ([`4cf75c8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4cf75c8))

- Remove file existence check for read files in seqkitreadstats
  ([`d07d1ea`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d07d1ea))

- Add metadata output in json format for reference processing
  ([`726adf7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/726adf7))

- Implement rasusa downsampling classes and atruntime wrapper for deferred execution
  ([`f405d50`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f405d50))

- Add rasusa dependency for read downsampling with citation
  ([`a2586f2`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a2586f2))

- Refactor logging in bugcatchinggroup and pipeline classes to use new logger methods
  ([`4ba6896`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4ba6896))

- Add tests for seqkitreadstats, seqkitreadstatsbasic, and seqkitreadstatsdetailed stages
  ([`c49a8d8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c49a8d8))

- Implement seqkitreadstats classes for generating read statistics
  ([`e9a5ff9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e9a5ff9))

- Add default output path and update reference option for short cli
  ([`9a37898`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/9a37898))

- Add seqkit dependency
  ([`3644aa0`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3644aa0))


### Refactoring

- Rename cleanreadsfastp to fastpcleanreads for consistency
  ([`7343465`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7343465))


### Chores

- Refactor cli into flat structure
  ([`9c83ed6`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/9c83ed6))


## v1.2.0 (2025-08-01)

### Features

- Add read cleaning option with fastp and implement cleanreadsfastp stage
  ([`c1aaafe`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c1aaafe))

- Add fastp dependency for read cleaning
  ([`ed424f1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/ed424f1))

- Add ssh remote configuration for glasscandle deployment
  ([`89b18bb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/89b18bb))


### Bug Fixes

- Restore import of alignmentfilter in snp calling pipeline
  ([`a37298a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a37298a))

- Update samtools version and add glasscandle dependency
  ([`6d8e4cf`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6d8e4cf))

- Add minimap2 to conda dependencies in watcher
  ([`ec78cc2`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/ec78cc2))

- Update minimap2 dependency format in pyproject.toml
  ([`c3cf074`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c3cf074))

- Enable slack notifications in watcher initialization
  ([`b980bc3`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b980bc3))

- Update ssh key usage and correct github actions user name in glasscandle workflow
  ([`035b0f1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/035b0f1))

- Comment out ssh remote configuration and adjust watcher initialization in glasscandle workflow
  ([`7e4afce`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7e4afce))

- Add permissions for id-token and contents in glasscandle workflow
  ([`d0214f5`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d0214f5))

- Update slack webhook url and adjust watcher script path in glasscandle workflow
  ([`6d62f55`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6d62f55))

- Update permissions to include pull-requests in glasscandle workflow
  ([`a23cd80`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a23cd80))

- Update git configuration and remote setup in glasscandle workflow
  ([`fa26b8d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fa26b8d))


### Refactoring

- Dynamically load dependencies from pyproject.toml in watcher
  ([`5036a3a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/5036a3a))


## v1.1.0 (2025-07-31)

### Features

- Implement alignmentfilter class for bam file filtering and add stats module
  ([`740346e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/740346e))

- Add alignfilter class for bam file filtering using samtools
  ([`97bed69`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/97bed69))

- Add options to skip and check dependency checks in cli
  ([`f05bab4`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f05bab4))

- Update glasscandle workflow to use watcher.py and add watcher script for dependency management
  ([`8caff61`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8caff61))

- Add glasscandle workflow and initial script for dependency management
  ([`34ac68b`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/34ac68b))

- Add minimap2 support and new short command for snp calling pipeline
  ([`f3959d1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f3959d1))


### Bug Fixes

- Change bam file type from str to path in aligneroutput and alignmentfilteroutput classes
  ([`1a7663e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1a7663e))

- Improve dependency validation by skipping already checked dependencies
  ([`e0ee5a2`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e0ee5a2))

- Add additional context section to bug report template in bug_catcher.py
  ([`33416c4`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/33416c4))

- Remove commented line for triggering github test in watcher.py
  ([`7843c6c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7843c6c))

- Adjust minimap2 command to remove unnecessary option for short reads
  ([`e0ecdcf`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/e0ecdcf))

- Update import statement to use glasscandle module in watcher.py
  ([`40f3c34`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/40f3c34))

- Enhance file format detection for reference processing
  ([`45d18ad`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/45d18ad))

- Remove unused 'outdir' field from basestage class
  ([`6dc0dc9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6dc0dc9))

- Update command example in documentation to use 'run' instead of 'snps'
  ([`6d44d2b`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6d44d2b))

- Update help description for --tmpdir option to remove redundancy
  ([`d2dca53`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d2dca53))

- Improve logging messages for dependency checks and working directory setup
  ([`6d0a762`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6d0a762))


### Refactoring

- Clean up glasscandle workflow configuration and improve comments
  ([`0860d79`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0860d79))


## v1.0.1 (2025-06-05)

### Bug Fixes

- Reorder global options for clarity and improve help descriptions
  ([`dea245a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/dea245a))

- Update version command to simplify output and improve help option
  ([`455af9e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/455af9e))

- Enhance help option settings and improve command-line options formatting
  ([`7368d87`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7368d87))


## v1.0.0 (2025-06-05)

### Features

- Implement bug-catching group with pre-filled bug report template
  ([`09ed70f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/09ed70f))

- Add horizontal rule function and integrate into pipeline logging
  ([`53821de`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/53821de))

- Update global options for cpu and ram configuration, enhance output directory handling
  ([`7852dad`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7852dad))

- Enhance logging in pipeline run method to improve stage visibility
  ([`6184df3`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6184df3))

- Add global options for cli commands
  ([`efecb87`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/efecb87))


### Bug Fixes

- Add color to horizontal rule in bug report template output
  ([`507deb7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/507deb7))

- Remove unused import of get_terminal_size in bugcatchinggroup
  ([`02749b5`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/02749b5))

- Update url references to use github_url and docs_url in pipeline and about files
  ([`2dc1e4c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2dc1e4c))

- Update bug report and feature request templates
  ([`f051223`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f051223))

- Remove conflicting field outdir
  ([`850d512`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/850d512))

- Update benchmark command to use 'run' instead of 'snps'
  ([`057c4a3`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/057c4a3))


### Chores

- Update readme
  ([`62e9dde`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/62e9dde))

- Add test for bugcatchinggroup to verify bug report output on exception
  ([`19de0df`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/19de0df))

- Rename tests
  ([`34d0a75`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/34d0a75))

- Add __pycache__ to .gitignore
  ([`9a9dc2a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/9a9dc2a))

- Remove obsolete google site verification file
  ([`ad75b52`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/ad75b52))

- Refactor cli
  ([`ed75edc`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/ed75edc))

## v0.4.0 (2025-04-10)

### Bug Fixes

- Update PyPI badge links to reflect the correct project name
  ([`187c3e9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/187c3e9f9f05c7af73ee5845128075ad4bbf26a5))

- Update README and .gitignore for consistency and clarity
  ([`1c6b0ed`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1c6b0ed938beb3b021f7ea1193052f39431c6e21))

### Features

- Add pypi job for building and publishing packages
  ([`c50cc9f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c50cc9f6ecabcf44268e13f56339884b743524c5))

- Update project name and enhance metadata in pyproject.toml
  ([`3be44bb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3be44bbbe21ef9ab89b4a6cb45aa2cdae9ce3f83))


## v0.3.9 (2025-04-10)

### Bug Fixes

- Update pixi-pack-install-script version to v2
  ([`6489690`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/64896908c9eaf500b9091e14fd5160dedeadf725))


## v0.3.8 (2025-04-10)

### Bug Fixes

- Update checkout step to use release tag and bump pixi-pack-install-script version
  ([`4d481a7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4d481a70fa77209663ff04bc968884b170b9976a))


## v0.3.7 (2025-04-10)

### Bug Fixes

- Add support for linux-aarch64 platform in release workflow
  ([`cb021ce`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/cb021cea133dff5eacd559562cdda53a3cc2f0f7))

- Update pixi-pack action version and adjust release workflow outputs
  ([`72c7fac`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/72c7fac91a9c874f8d89964b7f1c69f796a1406a))

### Documentation

- Add installation script command to README
  ([`71d7cc7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/71d7cc700dacfcdce47282181928dd58e200985c))


## v0.3.6 (2025-04-10)

### Bug Fixes

- Update release workflow to capture and use the release tag output
  ([`0738d37`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0738d37bd494d71cd9b6afe2f8351c9f2d5f7a35))


## v0.3.5 (2025-04-10)

### Bug Fixes

- Remove 'force' input from release workflow and adjust job conditions
  ([`a4ebebc`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a4ebebc4c5239d36c7dbb66f472b646c8980c688))

- Update pixi-pack action to v5 and modify build channels
  ([`0bb02a1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0bb02a1f93ece23de566799df057582d95fa83e3))


## v0.3.4 (2025-04-10)

### Bug Fixes

- Add 'force' input to release workflow to skip tests and benchmarks
  ([`c44bbd7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c44bbd7de66e53dbb8006e7ddf7142b2267afc92))

- Remove unused release commands and update semantic release configuration
  ([`8cdb251`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8cdb2517fd5ecbdfcb6e4e6528976e0275acec41))

- Update release job conditions to include 'force' input for immediate release
  ([`212ad06`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/212ad06e1d5ab72e1f2918b509383ba7f23c8356))

- Update SSH remote configuration and add repository URL for semantic release
  ([`6adbfb8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6adbfb89808fa87b43c1c77335df6e5d0dd85dba))


## v0.3.3 (2025-04-09)

### Bug Fixes

- Correct file extensions from .yml to .yaml in workflow configuration
  ([`af66615`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/af6661544c1a7ca6aded178eab87fcc2ddce1ca9))

- Correct file extensions in workflow configuration
  ([`00c9495`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/00c9495624b7542e2df5f0d5509ee3242e656b98))

- Correct permissions configuration in benchmark workflow
  ([`30430bf`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/30430bfeffa554c91e63c5823a1c8447196d510d))

- Enhance SSH configuration for deployment in release workflow
  ([`7a180e9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7a180e9ac564c4b1d4486f73d407c1a6d66cd037))

- Manual release workflow
  ([`2a6d406`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2a6d4065cd61ed90f973e7275411f4dfcefcd57b))

- Remove unnecessary checkout options in workflow configurations
  ([`a0e065d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a0e065dfb26cc404e4fbef24fc55909d7e678c3c))

- Remove unnecessary git add and commit commands from release workflow
  ([`b2dbead`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b2dbead5319eeed504a2f1e770fabd13eb50c93b))

- Remove unnecessary Git ref input from workflow configuration
  ([`50dba48`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/50dba48a7389ccdcdb8f46cfdadfb81a38f3943d))

- Remove unnecessary secrets from workflow configuration
  ([`054c62c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/054c62ccad2a8ab9668a99592bbe4de2b6eb05ce))

- Remove unnecessary secrets inheritance from benchmark workflow
  ([`f9d374a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f9d374a051aae6e96423e8512d4da84139ae2d33))

- Update job dependencies in release workflow
  ([`8cba440`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8cba4408ccaa0a63d61cca20cf7215b2fead26ef))

- Update permissions and secrets configuration in workflows
  ([`2cf045c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2cf045ce70a776152368ea3e818abddbdaf5160c))

- Update permissions to write for contents in workflows
  ([`f423c10`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f423c1072ea33c1b390095c2f5849e6458d586bb))

- Update release command to prevent pushing during versioning
  ([`904903e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/904903e7b78f27e4b818ee21cc77960b7431b8d7))

- Update release workflow to use semantic-release and rename jobs
  ([`6bb0728`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6bb0728be26b4e4c71307ec74e464a59440307ec))

### Chores

- Update commit message format for release step
  ([`ea4ced7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/ea4ced70efb1c46269f952d0f6aff9bcae70751c))

### Continuous Integration

- Add build and release jobs to workflow
  ([`7dd8d9d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7dd8d9d24588f6d76eca823173c83a538f6975db))

### Refactoring

- Streamline CI workflows by modularizing tests and benchmarks
  ([`87ba582`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/87ba582965c7454293dd5f9a6d7c5377b5769465))


## v0.3.2 (2025-04-09)

### Bug Fixes

- Push changes
  ([`11816a9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/11816a9e9971251c747a85053b12aab89a1cdf04))


## v0.3.1 (2025-04-09)

### Bug Fixes

- Commit release
  ([`4d494b8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4d494b8998b27c675db78af3e87b8b72f548d39c))


## v0.3.0 (2025-04-09)

### Chores

- Remove benchmark workflow
  ([`c835035`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c8350356cc3a70d7b1636173c29018d58a72dffa))

### Features

- Add benchmark to release workflow
  ([`0b9c3ff`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0b9c3ff52b09dce9f4ff348245a250bc4e533579))

- Update release
  ([`6c3370d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6c3370d4116b345a7c5e53212b549a586819c0d6))


## v0.2.1 (2025-04-09)

### Bug Fixes

- Update release.yaml
  ([`36ffa42`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/36ffa42f5c47e92b98d5a48a2156f458f3c5e4c1))


## v0.2.0 (2025-04-09)

### Features

- Add ssh key
  ([`814c2d8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/814c2d8c80aee7006c07b2e243cfcc8d704eed18))

- Add token to release workflow
  ([`2e3b879`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2e3b8796eca56869580b45d8716c774e93f88ac9))


## v0.1.0 (2025-04-09)

### Bug Fixes

- Replace semantic-release action
  ([`005dab1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/005dab184ff80719f55ca51b40224a90a2ad9c06))

### Features

- Add manual release workflow
  ([`d01c77a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d01c77a97b1f772a605b8237f2efb5d34774113d))


## v0.1.0-rc.1 (2025-04-09)

### Bug Fixes

- Rename mkdocs.yml to mkdocs.yaml
  ([`fbbf84f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fbbf84fb0a811bec5dc61f2fa8b94c7d519b8191))

- Rename workflow file for PRs
  ([`7ff1490`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7ff1490310248fb99f3643e6873c63f04b7c9cd4))

### Features

- Add release workflow
  ([`dba2b3a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/dba2b3a790a67a68e4019cc7c6b61610575a8bda))
