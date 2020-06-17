
# codon-degeneracy
[![Python application](https://github.com/nickmachnik/codon-degeneracy/workflows/Python%20application/badge.svg)](https://pypi.org/project/codon-degeneracy)
![License](https://img.shields.io/github/license/nickmachnik/codon-degeneracy)

This python package provides routines for the extraction of [degenerate sites](https://en.wikipedia.org/wiki/Codon_degeneracy) from sequences and alignments. The latter is particularly useful for estimations of rates of neutral evolution.

## Dependencies

This code uses [biopython](https://biopython.org/) and [scikit-bio](http://scikit-bio.org/) internally. In order to installl via pip, [numpy](https://numpy.org/) has to be installed.

## Installing

Simply clone this repo:

```
git clone https://github.com/nickmachnik/codon-degeneracy.git [TARGET DIR]
```

and then install  via pip
```
pip install [TARGET DIR]
```

or install directly from PyPI (this won't include unreleased changes as specified in the [changelog](CHANGELOG.md)):
```
pip install codon-degeneracy
```

## Testing

Test the cloned package:
```
cd [TARGET DIR]
python -m unittest
```

## Usage

There are more useful and well documented functions under the hood than shown here, which I enourage to explore by browsing the code.

### Counting substitutions per four fold degenerate site

One of the main features of the package is the counting of neutral substitutions at four fold degenerate sites.
This is best done with known orthologue pairs between species.
`substitution_rate_at_ffds` provides that functionality and is easy to use like so:
```python
from codon_degeneracy import substitution_rate_at_ffds as nsr
seq_a = (
    "ATACCCATGGCCAACCTCCTACTCCTCATTGTACCCATTC"
    "TAATCGCAATGGCATTCCTAATGCTTACCGAACGA")
seq_b = (
    "ATGACCACAGTAAATCTCCTACTTATAATCATACCCACAT"
    "TAGCCGCCATAGCATTTCTCACACTCGTTGAACGA")
(number_of_substitutions, number_of_sites), (orf_a, orf_b) = nsr(
    # the input sequences
    seq_a,
    seq_b,
    # NCBI codon table names as used in Bio.Data.CodonTable
    "Vertebrate Mitochondrial",
    "Vertebrate Mitochondrial")
```
The ORFs returned are there for sanity checks. The default behaviour is to select the first ATG codon
as start.

> NOTE: The numbers of neutral substitutions per site reported by this function are merely a lower bound,
> as they do not include the possibility of multiple substitutions per site.

### Substitutions at four fold degenerate sites separated by CpG context

In certain situations, it may be useful to differentiate between four fold degenerate sites
that could potentially exist in a CpG context and could therefore exhibit an elevated
mutation rate and those that do not. `substitutions_per_ffds_by_cpg_context` provides that
functionality.
It differentiates between four CpG contexts. Sites that are:
    - preceded by C and not followed by G (nonCpG)
    - preceded by C but not followed by G (postC)
    - followed by G but not preceded by C (preG)
    - preceded by C and followed by G (postCpreG)

> Note: the number of sites considered here may be substantially lower than
> in `substitutions_per_ffds`, as this function requires the sites
> preceeding and following site of a four fold degenerate site
> to be identical in the two aligned sequences.

The function can be used in exactly the same that is shown for `substitutions_per_ffds` above.

## License

MIT license ([LICENSE](LICENSE.txt) or https://opensource.org/licenses/MIT)

<!-- 
End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

 -->