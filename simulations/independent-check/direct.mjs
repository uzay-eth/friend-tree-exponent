import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createInterface } from 'node:readline';

const WORD = 2 ** 32;
const rotate = (x, bits) => (x << bits) | (x >>> (32 - bits));

class Rng {
  constructor(seed) {
    const nextSeedWord = () => {
      seed = (seed + 0x9e3779b9) >>> 0;
      let word = Math.imul(seed ^ (seed >>> 16), 0x85ebca6b);
      word = Math.imul(word ^ (word >>> 13), 0xc2b2ae35);
      return word ^ (word >>> 16);
    };
    [this.a, this.b, this.c, this.d] = Array.from({ length: 4 }, nextSeedWord);
  }

  word() {  // xoshiro128**
    const result = Math.imul(rotate(Math.imul(this.b, 5), 7), 9) >>> 0;
    const shift = this.b << 9;
    this.c ^= this.a;
    this.d ^= this.b;
    this.b ^= this.c;
    this.a ^= this.d;
    this.c ^= shift;
    this.d = rotate(this.d, 11);
    return result;
  }

  below(n) {
    const limit = WORD - WORD % n;
    let word;
    do { word = this.word(); } while (word >= limit);
    return word % n;
  }

  unit() {  // 53 random bits, yielding U in (0, 1]
    return ((this.word() >>> 6) * 2 ** 27 + (this.word() >>> 5) + 1) / 2 ** 53;
  }
}

class Tree {
  extra = [];

  grow(size, rng, observe) {
    assert(size >= 2 && size < WORD, `unsupported size ${size}`);
    if (!this.degree || size > this.degree.length) {
      const capacity = 2 ** Math.ceil(Math.log2(size));
      this.degree = new Uint32Array(capacity);
      this.first = new Uint32Array(capacity);
    }
    const { degree, first, extra } = this;
    first[0] = 1;
    first[1] = 0;
    degree[0] = degree[1] = 1;
    let core = 0;
    observe?.(2, core);

    for (let n = 2; n < size; n++) {
      const vertex = rng.below(n);
      const index = rng.below(degree[vertex]);
      const target = index === 0 ? first[vertex] : extra[vertex][index - 1];
      first[n] = target;
      degree[n] = 1;
      (extra[target] ??= [])[degree[target] - 1] = n;
      if (degree[target]++ === 1) core++;
      observe?.(n + 1, core);
    }
    return core;
  }
}

class Mean {
  count = 0;
  mean = 0;
  residualSquares = 0;

  add(x) {
    const delta = x - this.mean;
    this.mean += delta / ++this.count;
    this.residualSquares += delta * (x - this.mean);
  }

  get se() {
    return Math.sqrt(this.residualSquares / (this.count * (this.count - 1)));
  }
}

function exact(table, runs, seed) {
  const reference = [0];  // n=2: both endpoints have degree one
  for (const row of readFileSync(table, 'utf8').trim().split(/\r?\n/)) {
    const [size, , fraction] = row.trim().split(/\s+/);
    assert.equal(Number(size), reference.length + 2);
    const [numerator, denominator = 1] = fraction.split('/').map(Number);
    reference.push(numerator / denominator);
  }
  assert.equal(reference.length, 20, 'expected all reference sizes 2 through 21');
  const estimates = reference.map(() => new Mean());
  const tree = new Tree();
  const rng = new Rng(seed);
  const observe = (size, core) => estimates[size - 2].add(core);
  for (let run = 0; run < runs; run++) tree.grow(21, rng, observe);

  console.log('n\texact\tdirect\tSE\tz');
  for (const [index, estimate] of estimates.entries()) {
    const expected = reference[index];
    const delta = estimate.mean - expected;
    const z = estimate.se ? delta / estimate.se : delta === 0 ? 0 : Infinity;
    console.log([index + 2, fmt(expected), fmt(estimate.mean), fmt(estimate.se), fmt(z, 3)].join('\t'));
  }
}

async function compare(time, runs, treeSeed, sizeSeed) {
  const fast = [new Mean(), new Mean()];
  let count = 0;
  for await (const row of createInterface({ input: process.stdin, crlfDelay: Infinity })) {
    if (!row.trim() || row.startsWith('#')) continue;
    const [id, t, size, core] = row.trim().split(/\s+/, 4).map(Number);
    assert(id === count && t === time, `unexpected checkpoint on data row ${count + 1}`);
    assert(Number.isSafeInteger(size) && Number.isSafeInteger(core) && core >= 1 && core <= size - 2,
      `invalid size or core on data row ${count + 1}`);
    fast[0].add(size);
    fast[1].add(core);
    count++;
  }
  assert.equal(count, runs, 'wrong number of CT observations');

  const direct = [new Mean(), new Mean()];
  const tree = new Tree();
  const treeRng = new Rng(treeSeed);
  const sizeRng = new Rng(sizeSeed);
  const p = Math.exp(-time);
  const logFailure = Math.log1p(-p);
  for (let run = 0; run < runs; run++) {
    let size = 0;
    for (let founder = 0; founder < 3; founder++) {
      size += 1 + Math.floor(Math.log(sizeRng.unit()) / logFailure);
    }
    direct[0].add(size);
    direct[1].add(tree.grow(size, treeRng));
  }

  console.log(`t=${time}, runs per method=${runs}, Yule mean=${fmt(3 / p)}, Yule size SE=${fmt(Math.sqrt(3 * (1 - p) / runs) / p)}`);
  console.log('quantity\tdirect\tSE\tfast_CT\tSE\tdirect-fast\tSE_diff\tz');
  for (const [index, name] of ['N', 'C'].entries()) {
    const a = direct[index];
    const b = fast[index];
    const diff = a.mean - b.mean;
    const se = Math.hypot(a.se, b.se);
    console.log([name, fmt(a.mean), fmt(a.se), fmt(b.mean), fmt(b.se), fmt(diff), fmt(se), fmt(diff / se, 3)].join('\t'));
  }
}

const fmt = (value, digits = 6) => value.toFixed(digits);
const [mode, input, ...args] = process.argv.slice(2);
const [runs, ...seeds] = args.map(Number);
assert(input && Number.isSafeInteger(runs) && runs >= 2 &&
  seeds.every(seed => Number.isInteger(seed) && seed >= 0 && seed < WORD),
  'usage: node direct.mjs exact TABLE RUNS SEED | compare TIME RUNS TREE_SEED SIZE_SEED');
if (mode === 'exact' && seeds.length === 1) {
  exact(input, runs, seeds[0]);
} else if (mode === 'compare' && seeds.length === 2) {
  const time = Number(input);
  assert(Number.isInteger(time) && time > 0 && seeds[0] !== seeds[1]);
  await compare(time, runs, ...seeds);
} else {
  throw new Error('expected exact with one seed or compare with two seeds');
}
