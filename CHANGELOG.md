# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased] 
### Added
- :arrow_up: Added Project.toml
- Automa v1 compatibility: Upgraded the Automa dependency to "1", enabling the new Automav1 API.
- BioGenerics support: Imported metadata functions from BioGenerics to unify VCF/BCF header handling.
- TranscodingStreams integration: Added using TranscodingStreams for more efficient stream transformations in VCF/BCF readers.
- New VCF record reader: Introduced `src/vcf/readrecord.jl` to encapsulate record parsing logic.

### Changed
- Streamlined imports: Limited BioSequences imports, upgraded BGZFStreams and BufferedStreams usage, and replaced BioCore I/O types with BioGenerics abstractions.
- BCF reader refactoring: Transitioned Reader to subtype BioGenerics.IO.AbstractReader, centralized exception handling, and cleaned up parse logic.
- VCF header & metainfo: Fixed header parsing and improved metainfo tag/value functions to use BioGenerics APIs.
- Project.toml targets: Reorganized `[extras]` and `[targets]` sections.

### Removed
- Deprecated dependencies: Dropped older Automa versions (0.7, 0.8) and obsolete IO imports from BioCore.


## [0.4.0] - 2018-11-22
### Added
- :arrow_up: Support for julia v0.7 / v1.0

### Changed
- Fixed an issue with parsing VCF files from recent GATK releases.

### Removed
- :exclamation: Dropped support for julia v0.6 and v0.7

## [0.3.2] - 2018-07-25
### Added
- Project files: HUMANS.md & CHANGELOG.md

### Changed
- Fixed an error in the `BCF.n_samples` getter method that was causing it to return incorrect values.
- Updated to README.md to latest style. 

## [0.3.1] - 2017-10-10
### Added
- Method computing the average number of mutations, that should have gone into
  version 0.3.0.

## [0.3.0] - 2017-10-05
### Added
- Number of segregating sites method.
- Methods for computing allele frequencies and nucleotide diversity from sequences.

### Changed
- Documentation updated.
- API changes to dN/dS methods.

## [0.2.0] - 2017-07-31
### Added
- MASH distances.
- NG86 method of dN/dS computation.
- Proportion distance.

### Removed
- :exclamation: Dropped support for julia v0.5.

## [0.1.0] - 2017-06-23

Initial release, split from [Bio.jl](https://github.com/BioJulia/Bio.jl).

[Unreleased]: https://github.com/BioJulia/GeneticVariation.jl/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/BioJulia/GeneticVariation.jl/compare/v0.3.1...v0.4.0
[0.3.2]: https://github.com/BioJulia/GeneticVariation.jl/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/BioJulia/GeneticVariation.jl/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/BioJulia/GeneticVariation.jl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/BioJulia/GeneticVariation.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/BioJulia/GeneticVariation.jl/tree/v0.1.0