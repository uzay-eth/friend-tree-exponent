# Cold cross-family review of repository packaging — Claude Opus 5.5

- **Referee:** Claude Opus 5.5 (`claude-opus-5-5`, Claude Code harness), cold, different model family from packager (Gemini 4 Neon).
- **Date:** 2026-11-22.
- **Fresh clone run:** `bash verify/run_all.sh` exited `0` in `635 s` with all 7 verification suites passing (`motif-bracket` 21 s, `banded-upper` 72 s, `banded-upper-audit` 160 s, `transfer-bounds` 91 s, `almost-sure-lower` 105 s, `cross-charging` 3 s, `polynomial-lower` 176 s).

## Round 1 (commit `06284f6`)

- **README vs proofs:** Main mathematical claims match the proofs. Requested three wording/provenance fixes:
  1. Replace the tree counts not backed by the packaged attack files (`74,710 larger ones`, `121,603 labelled trees`) with the exact figures in `reviews/attack-bracket-motifs-claude.md` (`34,776` trees total, including all `32,972` rooted trees with $n \le 14$) and `reviews/attack-upper-bound-banded-claude.md` (`2,336` random histories + ~`50` corner histories, `133,864` checked times, `1,532,992` vertex rows), and state accurately that the three packaged attack reports were produced by Claude Opus 5.5 (under the Fable lead).
  2. Note plainly that the earlier `[0.543, 0.586]` baseline write-up (referenced by Fable's `0.543` application and Carbon's proof) is not packaged separately here (only its exact verification scripts and two cold reviews are included).
  3. Clarify historical paths and relative links in verbatim `proofs/` (including `REPORT.md §§3.3–3.5` in `proofs/polynomial-scale-lower.md`, which refers to the earlier baseline report rather than this repository's top-level `REPORT.md`) and convert dead relative markdown links in `reviews/` to plain code spans.
- **Exact numbers:** All recomputed in Python `Fraction` and `mpmath` and confirmed exact (all nine margins, prefactors, $\Gamma$ constants, Gautschi constant, D(c) bounds, literature endpoints, recipe SHA-256, table sizes, check counts).
- **Internal strings & git hygiene:** Clean for internal session IDs, machine paths, old branch names, hostnames, tokens, personal emails, and `__pycache__`. Fresh single-commit history with author/committer `Can Uzay`, no `Co-Authored-By`, repository `PRIVATE`.
- **Fresh-clone `verify/run_all.sh`:** `PASS` (all 7 suites, exit 0, `635 s`).

**Round 1 Verdict:** `CHANGES_REQUESTED` (4 packaging/wording items above; zero mathematical or verification errors).
