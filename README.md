[![Licence](https://img.shields.io/badge/License-MIT-9cf.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/csvtag/pytest.yml?branch=main&label=Test&color=brightgreen)](https://github.com/akikuno/csvtag/actions)


# csvtag

`csvtag` is a toolkit for `csv` tag, a format of cs tag that supports inversion.

# Specification

This is essentially the same encoding as the [minimap2 cs tag](https://lh3.github.io/minimap2/minimap2.html#10), but with the one difference that **lowercase letters represent inversions**:

| Prefix  | Sequence                   | Description                     |
| ------- | -------------------------- | ------------------------------- |
| =       | [ACGTN]+                   | Identical sequence (long form)  |
| :       | [0-9]+                     | Identical sequence length       |
| *       | [ACGTN][ACGTN]             | Substitution: ref to query      |
| +       | [ACGTN]+                   | Insertion to the reference      |
| -       | [ACGTN]+                   | Deletion from the reference     |
| ~       | [ACGTN]{2}[0-9]+[ACGTN]{2} | Intron length and splice signal |
| [=+-*~] | [acgtn]                    | Inversion                       |

> [!IMPORTANT]
> All csv tags are based on the forward strand of the reference sequence (SAM FLAG is 0). The reverse strand is entirely reverse complemented.


# Definision of Inversion

- Inversion detection uses RNAME, POS, and FLAG from SAM files.
  - Sort alignments by RNAME and POS.
  - If there are 2 or fewer reads for a QNAME, there is no Inversion, so output the cstag in uppercase.
  - If there are 3 or more reads for a QNAME, detect Inversion.
    - Extract three alignments in order of ascending POS (first, second, third).
    - (1) If the reads of first, second, and third are within 50 bp of each other, and only the second is reverse-oriented, then the second is determined to be an Inversion.
    - Reverse complement the cs tag of the second and output it as a csv tag in lowercase.
    - If there are gaps between first, second, and third, fill them with `N`.
  - Apply the same process to any adjacent reads.


# Functions

- `csvtag.call()`: Generate a csv tag
- `csvtag.to_sequence()`: Reconstruct a reference subsequence from the alignment
<!-- - `csvtag.to_vcf()`: Generate an VCF representation -->
<!-- - `csvtag.to_html()`: Generate an HTML representation -->


# Usage

