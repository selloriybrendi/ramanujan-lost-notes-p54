# The two lost notes on page 54 of Ramanujan's lost notebook

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23030496.svg)](https://doi.org/10.5281/zenodo.23030496)

On page 54 of his lost notebook Ramanujan states two families of cubic theta-function identities and, in each family,
writes "see note" next to one further claim (A³ − B³ in the first family, A³ + B³ in the second). Both notes are lost
(Andrews–Berndt, *Ramanujan's Lost Notebook, Part II*, pp. 177–178).

We show that both missing identities are special cases of another claim on the same page, Ramanujan's cubic circular
summation (Entry 8.2.3), and that with the Borweins' cubic theta function a(q) = Σ_{m,n} q^(m²+mn+n²) they read

- A = f(−q⁷,−q⁸), B = q f(−q²,−q¹³):  **A³ − B³ = q² f³(−q³,−q¹²) + f(−q²,−q³) a(q⁵)**
- A = f(−q⁴,−q¹¹), B = q f(−q,−q¹⁴):  **q(A³ + B³) = f³(−q⁶,−q⁹) − f(−q,−q⁴) a(q⁵)**

The paper also gives an elementary second form, and computational evidence that the two missing combinations have no
one- or two-term theta-product representation in a large explicit search space.

- **Interactive checker:** https://otakhonkenjaev.com/ramanujan-lost-notes-p54/ — both identities coefficient by coefficient, numerical evaluation, Entry 8.2.3 for any a, b, and the product obstruction (runs in the browser)
- **Paper:** [`paper/lost_notes_p54.pdf`](paper/lost_notes_p54.pdf) (7 pages, version 2) and its LaTeX source
- **Programs:** `src/`

## Re-check every identity

```bash
cd src
python3 verifier.py 2000      # JAMI: 12/12 MOS  — all 12 identities, exact integer arithmetic, 2000 coefficients
python3 sonli_tekshir.py      # 50-digit check of every identity as printed (needs mpmath)
python3 ildiz.py              # ILDIZ: 17/17 MOS — theta series vs Jacobi triple product, Entry 8.2.3 as a
                              #   two-variable identity, a³ = b³ + c³
gp -q ildiz_gp.gp             # GP-ILDIZ: 15/15 MOS — independent PARI/GP re-implementation
python3 audit_v2.py           # AUDIT: 63/66 MOS, the 3 FARQ lines are deliberate negative controls
```

Outputs are in Uzbek: `MOS` = match, `FARQ` = mismatch, `JAMI` = total.

The remaining scripts reproduce the search of Section 4 (`mosla*.py`, `faktor*.py`, `davr_jadval.py`) and the
exploration that led to the circular-summation form (`ijod*.py`, `yoyish.py`, `tekshir.py`, `rekon.py`).

## Cite

Kenjaev, O. U. (2026). *The two lost notes on page 54 of Ramanujan's lost notebook* (Version 2). Zenodo.
https://doi.org/10.5281/zenodo.23030496

## Use of generative AI

Claude (Anthropic) via Claude Code was used to write and run the programs and to help prepare the text, under the
author's direction; every identity was checked by the programs above. The author takes full responsibility.

## License

Code: MIT. Paper: CC BY 4.0.
