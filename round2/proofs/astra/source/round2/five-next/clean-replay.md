# Cold production receipt for the frozen five passes

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 15: Claude on frozen composition and selected inequalities, HY4 Flash on components, GPT on full production; inputs include Claim 9, the five-pass replacement, and the unsigned K7/144 dual. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


The fresh production replay passed with **36,290,118 lower** and **32,703,071 upper** newborn ticks at scale $`2^{40}`$. All 25 input digests and all 50 rank counts matched. The two signed runs enumerated 4,431,601,576 rows in total. Elapsed time was 2,137.56 seconds and peak resident memory was 2,866,080 KiB. The proof SHA-256 is `dc81f75e223e17dee63154e9bf4d98a55c46a1583011680d0e8e64dcbaa67508`. Signed receipt and projection digests · Production stdout · Full command output

From the repository root, this command builds new binaries, tests the slopes and the small physical trees, and computes both signs in a fresh working directory. It also compares the integer transcript to the signed frozen reference.

~~~sh
bash round2/five-next/replay-split.sh
~~~

After the cold run, `python3 -B round2/five-next/cold_receipt.py` reconstructs its proof digest from the raw output and the fresh projections. The production receipt records the manifest digest, signed integer transcript digest, both primary projection digests, raw stdout digest, every rank count and the individual newborn integers. Its timing-independent proof digest is SHA-256 of the following ASCII records, each terminated by a line feed. The four digest values are lowercase hexadecimal. The first line distinguishes this receipt format from other calculations. The raw stdout digest is reported separately because the raw log includes timing and peak memory. Receipt generator

~~~text
FRIEND-TREE-Q4-FIVE-PASS-COLD-1
source_manifest <sha256>
integer_transcript <sha256>
primary_lower_projection <sha256>
primary_upper_projection <sha256>
~~~

The source manifest pins the two initial four-vertex tables, cutoff profiles, Bellman solver, exact checks and original reference. The integer transcript contains 72 lines of signed pass labels, rank counts, pass values and replacements. Each binary projection stores the complete primary five-set weight at 17,776,584 exact configurations across 24 shapes. The arithmetic mesh counts include zero rows where the specialized internal support stops. [Source manifest](split-inputs.sha256) · [Signed reference](split-reference.tsv) · Independent mathematical review

| Pass | Five sites | Four sites | Three sites | Two sites | One site | Rows per sign |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 2,021,486,872 | 56,279,691 | 70,491,064 | 32,031,570 | 1,061,159 | 2,181,350,356 |
| 1 | 7,312,312 | 6,449,232 | 2,975,049 | 976,140 | 52,278 | 17,765,011 |
| 2 | 1,177,302 | 381,038 | 2,975,049 | 976,140 | 52,278 | 5,561,807 |
| 3 | 1,177,302 | 381,038 | 2,975,049 | 976,140 | 52,278 | 5,561,807 |
| 4 | 1,177,302 | 381,038 | 2,975,049 | 976,140 | 52,278 | 5,561,807 |
