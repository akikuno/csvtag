HTML_HEADER = """<!DOCTYPE html>
    <html>
    <head>
    <style>
    h1 {
    font-family: Consolas, 'Courier New', monospace;
    color: #333;
    padding: 0.1em 0;
    border-top: solid 3px #333;
    border-bottom: solid 3px #333;
    }
    .p_seq {
    font-family: Consolas, 'Courier New', monospace;
    color: #585858;
    word-wrap: break-word;
    letter-spacing: 0.15em;
    }
    .p_legend {
    font-family: Consolas, 'Courier New', monospace;
    color: #585858;
    word-wrap: break-word;
    }

    .Ins {
    color: #333;
    background-color: #ee827c;
    font-weight: bold;
    }
    .Del {
    color: #333;
    background-color: #a0d8ef;
    font-weight: bold;
    }
    .Sub {
    color: #333;
    background-color: #98d98e;
    font-weight: bold;
    }
    .Splice {
    color: #333;
    background-color: #f8e58c;
    font-weight: bold;
    }
    .Unknown {
    color: #333;
    background-color: #c0c6c9;
    font-weight: bold;
    }
    .Inv {
    color: #333;
    border-top: 0.1em solid;
    border-bottom: 0.1em solid;
    font-weight: bold;
    }

    </style>
    </head>
    <body>
"""

HTML_LEGEND = """
<p class = "p_legend">
Labels:
<span class="Ins">Insertion</span>
<span class="Del">Deletion</span>
<span class="Sub">Substitution</span>
<span class="Splice">Splice</span>
<span class="Unknown">Unknown</span>
<span class="Inv">Inversion</span>
</p>
<hr>
"""

HTML_FOOTER = """
</body>
</html>
"""
