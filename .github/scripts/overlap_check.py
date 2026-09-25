#!/usr/bin/env python3
"""Report word sequences shared between documentation pages and source texts.

Prose in the documentation may draw on books and articles for what to explain,
never for how to say it. This script is the mechanical half of that rule: it
lists every run of N consecutive words (8 by default) that a page shares with
any of the source texts, after lowercasing and dropping punctuation, math and
code, so that a sentence carried over by accident is found before it is
committed.

The sources are plain-text files, typically extracted from the PDFs of the
cited works (for instance with PyMuPDF), which are not distributed with the
package. Extraction mangles ligatures, so a clean report is necessary rather
than sufficient: it catches what was copied, not what was paraphrased.

    python3 .github/scripts/overlap_check.py --sources a.txt b.txt -- page.md ...
    python3 .github/scripts/overlap_check.py -n 7 --sources a.txt -- docs/src/theory/*.md

The exit status is 1 when a shared sequence is found, 0 otherwise.
"""

import argparse
import re
import sys

WORD = re.compile(r"[a-z]+(?:'[a-z]+)?")


def strip_markdown(text):
    text = re.sub(r"```.*?```", " ", text, flags=re.S)       # fenced blocks
    text = re.sub(r"``.*?``", " ", text)                      # inline math
    text = re.sub(r"`[^`]*`", " ", text)                      # inline code
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)      # links -> text
    return text


def words(text):
    return WORD.findall(text.lower())


def ngrams(ws, n):
    return {tuple(ws[i:i + n]) for i in range(len(ws) - n + 1)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-n", type=int, default=8, help="sequence length in words")
    ap.add_argument("--sources", nargs="+", required=True)
    ap.add_argument("pages", nargs="+")
    args = ap.parse_args()

    source_grams = {}
    for path in args.sources:
        with open(path, encoding="utf-8", errors="ignore") as f:
            for g in ngrams(words(f.read()), args.n):
                source_grams.setdefault(g, path)

    found = 0
    for page in args.pages:
        with open(page, encoding="utf-8") as f:
            ws = words(strip_markdown(f.read()))
        hits = sorted({g for g in ngrams(ws, args.n) if g in source_grams})
        for g in hits:
            found += 1
            print(f"{page}: \"{' '.join(g)}\"  <- {source_grams[g]}")
    print(f"{found} shared sequence(s) of {args.n} words", file=sys.stderr)
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
