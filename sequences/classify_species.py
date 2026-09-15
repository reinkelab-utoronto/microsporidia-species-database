#!/usr/bin/env python3
"""classify_species.py -- THE canonical named/provisional rule (Python twin).

Single source of truth for the Python side. Import it everywhere:

    from classify_species import classify_microsporidia

Do not keep local copies of this function in any script. The R twin is
classify_species.R; check_classifier_parity.py verifies the two agree on
every Species Name in the database.

The rule: provisional if the name is empty; starts with an uninformative
prefix (uncultured/unnamed/unknown/undescribed/mutated/...); has fewer than
two informative tokens; carries sp./spp. among the first three tokens; or
has a purely numeric second token. Otherwise named. Nomenclatural acts
("sp. n.", "n. gen. sp.") are stripped BEFORE the sp. test, so
"Ameson earli sp. n." is named; "Microsporidium" remains a valid
collective-group genus. Duplicate tokens are KEPT (matching the R twin).
"""
import re
import unicodedata

__all__ = ["classify_microsporidia", "is_provisional",
           "RX_NOM_1", "RX_NOM_2", "UNINFORMATIVE_RX"]

RX_NOM_1 = r"\b(?:nov|n)\.(?:\s*,?\s*(?:n\.\s*)?(?:gen|sp|spp|comb)\.?)+"
RX_NOM_2 = (r"\b(?:gen|sp|spp|ssp|subsp|var|comb|stat)\.?\s*,?\s*"
            r"(?:et\s*|and\s*|&\s*)?(?:nov\.?(?![a-z])|n\.)")
UNINFORMATIVE_RX = (r"^(?:uncultured|unidentified|unclassified|environmental|"
                    r"microsporidian?|fungal|eukaryot\w*|unnamed|unknown|"
                    r"undescribed|mutated)\b")

_STOP = ("nov", "gen", "et", "and", "comb", "stat")


def _deaccent(text):
    s = unicodedata.normalize("NFKD", str(text))
    return "".join(c for c in s if not unicodedata.combining(c))


def classify_microsporidia(name):
    """Return 'named' or 'provisional' for one microsporidia species name."""
    s = _deaccent(str(name)).replace("\u00a0", " ").lower()
    s = re.sub(r"\(.*?\)", " ", s)
    s = re.sub(RX_NOM_1, " ", s)
    s = re.sub(RX_NOM_2, " ", s)
    s = re.sub(r"\b(?:cf|aff|nr)\.?\s+", " ", s)
    s = re.sub(r"[^a-z0-9\s]", " ", s)
    tok = [t for t in s.split() if t not in _STOP]
    if not tok:
        return "provisional"
    if re.match(UNINFORMATIVE_RX, " ".join(tok)):
        return "provisional"
    if len(tok) < 2 or any(t in ("sp", "spp") for t in tok[:3]):
        return "provisional"
    if re.fullmatch(r"[0-9]+", tok[1]):
        return "provisional"
    return "named"


def is_provisional(name):
    return classify_microsporidia(name) == "provisional"
