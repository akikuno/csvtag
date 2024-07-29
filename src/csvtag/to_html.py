from __future__ import annotations

import re

from csvtag.splitter import split_by_tag
from csvtag.template.html import HTML_FOOTER, HTML_HEADER, HTML_LEGEND


def append_mark_to_n(csv_tag: str) -> str:
    """Process each csv tag by adding specific markers `@` to `N`."""

    def append_mark(cs: str) -> str:
        if cs.startswith("N"):
            return "@" + cs
        elif re.match(r"^[ACGT]", cs):
            return "=" + cs
        return cs

    csv_tag = csv_tag.replace("=N", "N")
    csv_tag_processed = [append_mark(cs) for cs in re.split(r"(N+)", csv_tag) if cs and cs != "="]

    return "".join(csv_tag_processed)


def apply_css(cs: str, css_class: str) -> str:
    return f"<span class='{css_class}'>{cs.upper()}</span>"


def process_csv_tag(csv_tag: str) -> str:
    # Format csv_tag
    csv_tag = csv_tag.replace("cs:Z:", "")
    csv_tag_marked = append_mark_to_n(csv_tag)
    csv_tag_split = split_by_tag(csv_tag_marked)

    # Build html
    html_body = []
    idx = 0
    while idx < len(csv_tag_split):
        cs = csv_tag_split[idx]
        if cs.startswith("="):
            html_body.append(cs[1:])
        elif cs.startswith("@"):
            html_body.append(apply_css(cs[1:], "Unknown"))
        elif cs.startswith("*"):
            substitutions = [cs[2]]
            while idx < len(csv_tag_split) - 1 and csv_tag_split[idx + 1].startswith("*"):
                substitutions.append(csv_tag_split[idx + 1][2])
                idx += 1
            html_body.append(apply_css("".join(substitutions), "Sub"))
        elif cs.startswith("+"):
            html_body.append(apply_css(cs[1:], "Ins"))
        elif cs.startswith("-"):
            html_body.append(apply_css(cs[1:], "Del"))
        elif cs.startswith("~"):
            left, right = cs[1:3], cs[-2:]
            splice = "-" * (int(cs[3:-2]) - 4)
            html_body.append(apply_css(f"{left}{splice}{right}", "Splice"))
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
    html_body = process_csv_tag(csv_tag)
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
