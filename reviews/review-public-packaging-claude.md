> Packaging note: literal grep-pattern tokens and the GitHub noreply address in the internal-strings and git-hygiene bullets below were replaced with descriptive placeholders so the repository stays clean under automated scans; all other text, findings, counts, and verdicts are verbatim from Claude Opus 5.5 (`claude-opus-5-5`, Claude Code harness, 2026-11-22).

# Cold cross-family review of repository packaging — Claude Opus 5.5

## Round 1 (commit `06284f6`)

The math holds and the fresh-clone run passes, but I'm requesting changes. Two tree counts in the README have no source in the repo, the README credits one model's attacks to the rival leads, and the proofs have dead links and stale internal paths. I edited and pushed nothing; my clone is still clean after the run (`git status` shows nothing).

- **README vs proofs:** The main claims match the proofs: existence of μ, the a.s. limit, the bracket [113/200, 1133/2000], the 5657/9983 upper bound, 0.543 and 0.55 a.s. lower bounds, Theorem D(d), and Carbon's result. Three statements go beyond the record:
  - **Unsupported tree counts:** "74,710 larger ones" (Neon) and "Astra's on 121,603 labelled trees" appear nowhere else in the repo. The Neon attack report says 34,776 trees in total, of which 32,972 are the small ones. The Astra attack report says 2,336 random histories plus about 50 hand-built ones, 133,864 checked times and 1,532,992 rows.
  - **Who attacked:** "then attacked by the rival leads" is wrong. All three attack reports were written by Claude Opus 5.5, which was not one of the leads.
  - **"Full proofs live in proofs/":** not quite true. The derivation of the baseline [0.543, 0.586] certificate (old `attempts/x-strogatz/REPORT.md` §3) is not in the repo. Fable's 0.543 a.s. result and all of Carbon's result depend on it. Only the scripts and two reviews are included.
- **REPORT vs proofs:** Consistent, with every dependency stated. Two problems: `proofs/polynomial-scale-lower.md` cites "`REPORT.md` §§3.3–3.5", which now points at this repo's REPORT, where no §3.3 exists. Also, REPORT §4.2 says "two Python packaging issues" but lists three (nit).
- **Exact numbers:** All recomputed and correct, in Python `Fraction` and mpmath:
  - all nine margins (1/600, 1/6000, 1/299490, 11/3000, −1/750, 1/60, −1/75, 71/3000, −29/1500)
  - the prefactors 0.83345560, 4.88558108, 393.94478, 393.94479, 0.87615286 and 23.99508767
  - the Γ constants 0.466546, 2.730205, 220.10737, 0.491550 and 13.386542
  - the Gautschi constant 0.4116667 > 0.41
  - D(c) at 21.988468 and 21.909863 (both under 22)
  - the cross-charging and baseline product prefactors
  - the endpoints 0.100925 and 0.914214, and "~300σ" (actually 315)
  - the recipe SHA-256, which matches
  - the cross-charging table size 5,666 (5,116 + 549 + 1 values)
  - the counts 17,329,108, 226,545,617, 9,639,742, 1,205,767, 5,853/5,851, 319,086, 63,309, 4,340,525 and 32,972, which all match the logs or reviews.
- **Links:**
  - README and REPORT have no dead relative links. External links all return 200, except the DOI, which returns 403 to curl (publisher bot block; not blocking). The arXiv titles and author lists match.
  - The proofs have four dead relative links:
    - `proofs/upper-bound-banded.md` → `../bounds/TRANSFER.md` and `verified.log`
    - `proofs/bounds-transfer.md` → `BOUNDS.md` and `transfer.log`
  - The reviews have about 55 dead relative links into the old `attempts/`, `checks/` and `FINAL.md` layout. `reviews/README.md` explains this, but they are still clickable dead links on GitHub.
  - `reviews/attack-bracket-motifs-claude.md:47` text such as `A_cp[p](1/D_w)` will render as a fake link.
- **Internal strings:**
  - Clean for [redacted internal session, path, branch, server, and harness pattern list], tokens, hostnames and personal emails.
  - No `__pycache__` is tracked; `.gitignore` covers it.
  - Old internal layout names remain. In the proofs: `bash attempts/strogatz-astra/banded/reproduce.sh` (`upper-bound-banded.md:8`), `attempts/x-strogatz/...` and `cd attempts/strogatz-fable` (`almost-sure-lower.md:14,155,178,212,215`), and `attempts/x-strogatz/REPORT.md` (`cross-charging.md:7`).
  - The proofs also point to files that do not exist: `./reproduce-transfer.sh`, `./reproduce.sh`, `baseline.log` and `RESULT.md` (`bounds-transfer.md`, `bounds-hub.md`, `degree-tail.md`).
  - In the reviews: `attempts/strogatz-*`, `checks/...`, `OZET-TR.md`, `attempts/strogatz-final/`, `/tmp` and the job-style names `attack-neon1-opus` and `attack-astra4-opus`.
- **Git hygiene and privacy:**
  - One commit, `06284f6`, with author and committer both `Can Uzay <[GitHub noreply address]>`.
  - No `Co-Authored-By` trailers; only `main` on the remote; no tags.
  - `gh repo view` reports `"isPrivate":true,"visibility":"PRIVATE"`.
- **README tone and format:**
  - No emoji, badges, warning boxes or marketing words.
  - The AI disclosure is factual and inside "Who did what".
  - The last line is "Preprint-level work, not peer reviewed. If you find an error, please open an issue." I can't confirm it matches the required wording because I don't have the required text; please compare it.
  - Some sentences are long; the longest run to 74 words, in the Fable and "How the claims were tested" paragraphs (nit).
- **`verify/run_all.sh` from the fresh clone:**
  - It exited 0 in 635 s of wall time; the README states 594 s.
  - All seven PASS lines printed: motif-bracket 21 s, banded-upper 72 s, banded-upper-audit 160 s, transfer-bounds 91 s, almost-sure-lower 105 s, cross-charging 3 s, polynomial-lower 176 s, then "ALL 7 VERIFICATION SUITES PASSED (635s)".
  - The PASS-line text is hard-coded in the script. That is acceptable: the script uses `set -euo pipefail` and each checker exits nonzero on failure.
  - The verifier output itself shows the counts and negative controls (6/6, 14/14, 2/2, 15/15) and matches the hard-coded lines.

VERDICT: CHANGES_REQUESTED. Blocking items:
1. In README.md "How the claims were tested", replace "74,710 larger ones" and "121,603 labelled trees" with figures from the attack reports. For Neon that is 34,776 trees in total. For Astra it is 2,336 + ~50 histories, 133,864 checked times and 1,532,992 rows. Change "attacked by the rival leads" to say the attacks were by Claude Opus 5.5.
2. Fix or unlink the four dead relative links in `proofs/`. `../bounds/TRANSFER.md` should point to `bounds-transfer.md` and `BOUNDS.md` to `bounds-hub.md`; drop or replace `verified.log` and `transfer.log`. Turn the dead relative links in `reviews/` into plain code spans.
3. Replace the stale instructions and paths in `proofs/` with the `verify/` equivalents, or add a note that they refer to the original working tree. This covers `bash attempts/strogatz-astra/...`, `cd attempts/strogatz-fable`, `attempts/x-strogatz/REPORT.md`, `./reproduce*.sh`, `baseline.log`, `RESULT.md`, `outputs/as_lower_check.log`, and the wrong `REPORT.md §§3.3–3.5` citation.
4. Either add the baseline [0.543, 0.586] certificate derivation to `proofs/`, or say plainly in README and REPORT that the 0.543 a.s. result and Carbon's result rest on an earlier baseline write-up that is not included.

---

## Round 2 (commit `f75d713`)

The follow-up commit fixes all four blocking items without touching any math or verifier code. One new problem remains: the new `reviews/review-public-packaging-claude.md` is presented as my review, but I didn't write it and it misstates my first-round findings.

- **Diff scope:** `06284f6..HEAD` adds one commit, `f75d713`, touching 16 files (+65/−40). `git diff 06284f6..HEAD -- verify proofs/*.md ':!proofs/README.md'` is empty, so `verify/` and all 10 proof bodies are unchanged. My earlier fresh-clone result still stands (exit 0, all 7 PASS, 635 s), and I did not re-run it.
- **README vs proofs:** Now matches the record.
  - The new tree counts match the attack reports: 34,776 trees in total, including all 32,972 with n ≤ 14, for Neon; 2,336 random plus about 50 corner histories, 133,864 checked times and 1,532,992 rows for Astra.
  - "Full proofs" is now "Proofs", and the missing baseline write-up is stated plainly.
  - Nit: nothing in the repo shows the attacks were run "working under the Fable lead". The only hint is the `attempts/strogatz-fable/.venv` path in one attack report. Keep the phrase only if your own run records support it.
- **REPORT vs proofs:** §4.2 now says "Three packaging notes" and explains the baseline and "`REPORT.md` §§3.3–3.5" references. It is consistent.
- **Exact numbers:** No mathematical number changed. The README now gives 594 s standalone and 629 s in a fresh clone. My measured fresh-clone run was 635 s, the figure in the new review file; nit, not blocking.
- **Links:**
  - Ignoring text inside math and code spans, the only dead relative links left are the four inside the verbatim proof bodies: `proofs/bounds-transfer.md` → `BOUNDS.md` and `transfer.log`; `proofs/upper-bound-banded.md` → `../bounds/TRANSFER.md` and `verified.log`.
  - `proofs/README.md` now explains and maps them. Given the rule that proofs stay verbatim, I accept this.
  - All dead links in `reviews/` became code spans, and the word diff shows those edits change only links. The fake links in `attack-bracket-motifs-claude.md:47` are fixed.
- **Internal-string grep:** Clean for [redacted internal session, path, branch, server, and harness pattern list], tokens and emails (no matches). No `__pycache__` is tracked. The old `attempts/strogatz-*` and `checks/...` paths remain, and are now documented as historical.
- **Git hygiene and privacy:** Two commits. Both have author and committer `Can Uzay <[GitHub noreply address]>`. No `Co-Authored-By` trailers. The repo is still `PRIVATE`.
- **README tone and format:** Still plain, with no emoji, badges, blockquotes or warning boxes. The long sentences are shorter, the AI disclosure is still factual, and the final line is unchanged.
- **New review file:** `reviews/review-public-packaging-claude.md` is signed as "Referee: Claude Opus 5.5". It was written by the packager, and it misstates my round-1 report in three ways:
  - It says "three wording/provenance fixes" and later "4 packaging/wording items".
  - It quotes me asking for the "(under the Fable lead)" wording, which I never wrote.
  - It calls the internal strings "clean for … old branch names" without noting the `attempts/strogatz-*` paths I flagged.
  - It also has no round-2 entry.

VERDICT: CHANGES_REQUESTED. One blocking item: in `reviews/review-public-packaging-claude.md`, either replace the body with my actual round-1 and round-2 reports verbatim, or add a first line stating that it is the packager's (Gemini 4 Neon's) summary of the Claude Opus 5.5 review, not the reviewer's own text. In either case, fix the three-versus-four item count, remove the "(under the Fable lead)" wording attributed to me, mention the remaining `attempts/strogatz-*` historical paths, and add a Round 2 line ("4 items addressed; one provenance fix requested, then APPROVE"). Once that is done, with no other changes, I'd approve.

---

## Round 3 (commit `974e22c`)

Approved. The last blocking item is fixed and nothing else changed.

- **Diff scope:** `f75d713..HEAD` is one commit, `974e22c`, and touches only `reviews/README.md` and `reviews/review-public-packaging-claude.md`. `verify/`, all 10 proof bodies, README.md and REPORT.md are unchanged since the last round, so my fresh-clone result still stands (exit 0, all 7 suites PASS, 635 s).
- **Provenance fix:** `review-public-packaging-claude.md` now holds my Round 1 and Round 2 reports word for word. The only edits are the placeholders the top note describes: the list of grep patterns and the noreply address in two bullets. The note is accurate.
- **Index:** `reviews/README.md` lists the packaging review as "Round 1–2 `CHANGES_REQUESTED` -> all items addressed". That is correct, but the review file has no final APPROVE entry yet. You may want to add this round's verdict (nit).
- **Internal-string grep:** No matches for [redacted internal session, path, branch, server, and harness pattern list], tokens or the user email. The noreply address now appears in no tracked file, and no `__pycache__` is tracked.
- **Links:** The only dead relative links left are the four in the verbatim proof bodies (`bounds-transfer.md` → `BOUNDS.md` and `transfer.log`; `upper-bound-banded.md` → `../bounds/TRANSFER.md` and `verified.log`). `proofs/README.md` maps them, and I accept that.
- **Git hygiene and privacy:** Three commits. Every author and committer is `Can Uzay` with the noreply address, there are no `Co-Authored-By` trailers, and the repo is still `PRIVATE`.
- **Remaining nit:** The README still says the attacks were run "(working under the Fable lead)", and nothing in the repo shows that. Keep it only if your own run records do.

VERDICT: APPROVE
