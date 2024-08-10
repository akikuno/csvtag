from __future__ import annotations

import re

from csvtag.splitter import split_by_tag
from csvtag.template.html import HTML_FOOTER, HTML_HEADER, HTML_LEGEND


def _append_mark_to_n(csv_tag: str) -> str:
    """Process each csv tag by adding specific markers `@` to `N` or `n`."""

    def _append_mark(tag: str) -> str:
        if tag.startswith("N") or tag.startswith("n"):
            return "@" + tag
        elif re.match(r"^[ACGTacgt]", tag):
            return "=" + tag
        return tag

    csv_tag = csv_tag.replace("=N", "N").replace("=n", "n")
    csv_tag_processed = [_append_mark(tag) for tag in re.split(r"(N+|n+)", csv_tag) if tag and tag != "="]

    return "".join(csv_tag_processed)


def _apply_css(cs: str, css_class: str) -> str:
    return f"<span class='{css_class}'>{cs.upper()}</span>"


def _make_html_body(csv_tag: str) -> str:
    # Format csv_tag
    csv_tag_marked = _append_mark_to_n(csv_tag)
    csv_tag_split = list(split_by_tag(csv_tag_marked))

    # Build html
    html_body = []
    idx = 0
    while idx < len(csv_tag_split):
        tag = csv_tag_split[idx]
        operand, nucleotides = tag[0], tag[1:]
        if operand.startswith("="):
            html_body.append(nucleotides)
        elif operand.startswith("@"):
            html_body.append(_apply_css(nucleotides, "Unknown"))
        elif operand.startswith("*"):
            substitutions = [nucleotides[1]]
            while idx < len(csv_tag_split) - 1 and csv_tag_split[idx + 1].startswith("*"):
                substitutions.append(csv_tag_split[idx + 1][2])
                idx += 1
            html_body.append(_apply_css("".join(substitutions), "Sub"))
        elif operand.startswith("+"):
            html_body.append(_apply_css(nucleotides, "Ins"))
        elif operand.startswith("-"):
            html_body.append(_apply_css(nucleotides, "Del"))
        elif operand.startswith("~"):
            left, right = tag[1:3], tag[-2:]
            splice = "-" * (int(tag[3:-2]) - 4)
            html_body.append(_apply_css(f"{left}{splice}{right}", "Splice"))
        idx += 1

    return f"<p class='p_seq'>{''.join(html_body)}</p>"


def to_html(csv_tag: str, description: str = "") -> str:
    """Output HTML string showing a sequence with mutations colored
    Args:
        csv_tag (str): csv tag in the **long** format
        description (str): (optional) header information in the output string
    Return:
        HTML string
    Example:
        >>> import cstag
        >>> csv_tag = "=AC+ggg=T-acgt*at~gt10cg=GNNN"
        >>> description = "Example"
        >>> html_string = cstag.to_html(csv_tag, description)
    """

    description_str = f"<h1>{description}</h1>" if description else ""
    html_body = _make_html_body(csv_tag)
    report = "\n".join(
        [
            HTML_HEADER,
            description_str,
            HTML_LEGEND,
            html_body,
            HTML_FOOTER,
        ]
    )
    return report
