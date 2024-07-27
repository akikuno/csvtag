[![Licence](https://img.shields.io/badge/License-MIT-9cf.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)

# csv tag

csvtag is a toolkit for `csv` tag, a format of cs tag that supports inversion.

# Specification of `csv` tag

This is essentially the same encoding as the [minimap2 cs tag](https://lh3.github.io/minimap2/minimap2.html#10), but with the one difference that **lowercase letters represent inversions**:

| Prefix  | Sequence          | Description                  |
| --- | -------------- | ---------------------------- |
| =   | [ACGTN]+        | Identical sequence (long form)          |
| :   | [0-9]+        | Identical sequence length          |
| *   | [ACGTN][ACGTN] | Substitution: ref to query                 |
| +   | [ACGTN]+        | Insertion to the reference   |
| -   | [ACGTN]+        | Deletion from the reference  |
| ~	  | [ACGTN]{2}[0-9]+[ACGTN]{2} | Intron length and splice signal |
| [=+-*~]    | [acgtn]        | Inversion                    |

# Definision of Inversion

- Inversionの検出には、SAM/BAMファイルのPOSとFLAGを用います。
  - FLAGによりprimary readsとそのstrandnessを判定します。
  - POSにより、primary readsの10bp以内に隣接したreadsを抽出します
  - 隣接したreadsのstrandnessが異なる場合、Inversionと判定します。
  - 隣接したreadsの、さらに隣接したreadsに対しても、同様に処理します。
- readsの重複については、QUALが高い塩基を優先します
  - 同じQUALの場合には、参照ゲノムと合っている塩基を優先します


# Functions

- `csvtag.call()`: Generate a csv tag
- `csvtag.to_vcf()`: Generate an VCF representation
- `csvtag.to_html()`: Generate an HTML representation
- `csvtag.to_pdf()`: Produce an PDF file


# Usage

