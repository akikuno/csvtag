[![Licence](https://img.shields.io/badge/License-MIT-9cf.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)

# csv tag

# Operators

This is essentially the same encoding as the [minimap2 cs tag](https://lh3.github.io/minimap2/minimap2.html#10), but with the one difference that **lowercase letters represent Inversions**:

| Prefix  | Sequence          | Description                  |
| --- | -------------- | ---------------------------- |
| =   | [ACGTN]+        | Identical sequence (long form)          |
| :   | [0-9]+        | Identical sequence length          |
| *   | [ACGTN][ACGTN] | Substitution: ref to query                 |
| +   | [ACGTN]+        | Insertion to the reference   |
| -   | [ACGTN]+        | Deletion from the reference  |
| ~	  | [ACGTN]{2}[0-9]+[ACGTN]{2} | Intron length and splice signal |
| [=+-*~]    | [acgtn]        | Inversion                    |

# Features

- `csvtag.call()`: Generate a csv tag
- `csvtag.to_vcf()`: Generate an VCF representation
- `csvtag.to_html()`: Generate an HTML representation
- `csvtag.to_pdf()`: Produce an PDF file

# Usage

