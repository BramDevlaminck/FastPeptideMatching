# Suffix Array Mappings

![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/test.yml?logo=github)
![Codecov](https://img.shields.io/codecov/c/github/unipept/unipept-index?token=IZ75A2FY98&flag=sa-mappings&logo=codecov)
![Static Badge](https://img.shields.io/badge/doc-rustdoc-blue)

A suffix array returns a range of matches. This range has to be mapped to some informational data. The `sa-mappings` library 
offers a few utilities to build these mapping tables. The library offers functionality to map each SA output onto its 
proteins, taxonomy and functional annotations.

- `sa_mappings::taxonomy::TaxonAggregator` can aggregate a list of taxa.
- `sa_mappings::functionality::FunctionAggregator` can aggregate a list of functional annotations.
- `sa_mappings::proteins::Proteins` can map an SA entry to a protein and all its information
