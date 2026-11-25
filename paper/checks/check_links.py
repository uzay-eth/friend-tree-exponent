#!/usr/bin/env python3
"""Check local links in tracked Markdown, plus the paper's pinned source links."""
import argparse
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
import posixpath
import re
import subprocess
import unicodedata
from urllib.parse import unquote, urlsplit

from markdown_it import MarkdownIt

ROOT = Path(__file__).resolve().parents[2]
PIN = 'a0972085ad72ad5685ae7b424988648d1965e73c'
# The stochastic-order sources are cited through \\ordersrc at the commit that froze them.
ORDER_PIN = '62081cad01112c42c2011ac3eed11aa3aa665e09'
# Markdown may also link two archived trees; such links are checked against those trees.
ARCHIVE_PIN = 'c250a1f2fdd53ced91b7b485e94ac8a16106b9c7'
BASE_PIN = '56f2d807ba63d19e9cb0e56d7678ae50843d3f6a'
PAPER = 'paper/friend-tree-exponent.tex'
REPO_URL = 'https://github.com/uzay-eth/friend-tree-exponent/'
PINNED_URL = re.compile(re.escape(REPO_URL) + r'(blob|tree)/' + PIN + r'/(.+)\Z')
OTHER_PINNED_URL = re.compile(re.escape(REPO_URL) + r'(blob|tree)/([0-9a-f]{40})/(.+)\Z')
MD = MarkdownIt('commonmark', {'inline_definitions': True, 'store_labels': True}).enable(['table', 'strikethrough'])


def git(*args, data=None):
    return subprocess.run(['git', *args], cwd=ROOT, input=data, capture_output=True, check=True).stdout


def manifest(revision=None):
    cmd = ('ls-tree', '-rz', '--name-only', revision) if revision else ('ls-files', '-z')
    return set(git(*cmd).decode().strip('\0').split('\0'))


def blobs(revision, paths):
    query = ''.join(f'{revision}:{p}\n' for p in paths).encode()
    output = git('cat-file', '--batch', data=query)
    cursor = 0
    for path in paths:
        end = output.index(b'\n', cursor)
        _, kind, size = output[cursor:end].split()
        if kind != b'blob':
            raise ValueError(f'{revision}:{path} is not a blob')
        cursor = end + 1
        yield path, output[cursor:cursor + int(size)].decode('utf-8-sig')
        cursor += int(size) + 1
    if cursor != len(output):
        raise ValueError('unexpected git cat-file output')


class HTML(HTMLParser):
    def __init__(self):
        super().__init__()
        self.anchors = set()
        self.links = []

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if value is None:
                continue
            if key == 'id' or (tag == 'a' and key == 'name'):
                self.anchors.add(value)
                self.anchors.add('user-content-' + value)
            if (tag, key) in {('a', 'href'), ('img', 'src'), ('source', 'src')}:
                self.links.append((self.getpos()[0] - 1, value))


def slug(title):
    return ''.join(
        '-' if c == ' ' else c for c in title.lower()
        if c in ' -_' or unicodedata.category(c)[0] in 'LMN'
    )


def parse(text):
    headings, anchors, links = set(), set(), []
    tokens = MD.parse(text)
    line = 1
    for i, token in enumerate(tokens):
        if token.map:
            line = token.map[0] + 1
        if token.type == 'heading_open':
            title = ''.join(t.content for t in tokens[i + 1].children
                            if t.type in {'text', 'code_inline', 'image'})
            base = slug(title)
            anchor = base
            suffix = 0
            while anchor in headings:
                suffix += 1
                anchor = f'{base}-{suffix}'
            headings.add(anchor)
        if token.type == 'definition':
            links.append((line, token.meta['url'], 'reference definition'))
        for child in [token, *(token.children or [])]:
            if child.type in {'link_open', 'image'}:
                # The definition is checked once at its own line, including unused ones.
                if 'label' not in child.meta:
                    dest = child.attrGet('href' if child.type == 'link_open' else 'src')
                    if dest is not None:
                        links.append((line, dest, 'inline'))
            if child.type in {'html_block', 'html_inline'}:
                parser = HTML()
                parser.feed(child.content)
                anchors.update(parser.anchors)
                links.extend((line + offset, url, 'HTML') for offset, url in parser.links)
    return headings | anchors, links


def check(path, dest, files, directories, anchors, roots=False):
    try:
        url = urlsplit(dest)
    except ValueError:
        return 'malformed URL'
    if url.scheme or url.netloc:
        return None
    raw = unquote(url.path)
    fragment = unquote(url.fragment)
    if raw.startswith('/') or roots:
        target = posixpath.normpath(raw.lstrip('/'))
    elif raw:
        target = posixpath.normpath(posixpath.join(posixpath.dirname(path), raw))
    else:
        target = path
    if target not in files and target not in directories:
        return f'untracked or missing target {target!r}'
    # Non-Markdown targets are path-checked only (e.g. GitHub's code-line or PDF fragment syntax).
    if fragment and target in anchors and fragment not in anchors[target]:
        return f'missing Markdown anchor in {target!r}'
    return None


def scan(files, docs, pinned=None, pinned_anchors=None, others=None):
    # ``others`` maps further pinned commits to (manifest, anchors) for Markdown links.
    anchors, entries, errors = {}, [], []
    dirs = {str(d) for path in files for d in PurePosixPath(path).parents}
    for path, text in docs:
        if not path.lower().endswith('.md'):
            continue
        anchors[path], links = parse(text)
        entries.extend((path, *link) for link in links)
    local = 0
    pin_dirs = {str(d) for p in (pinned or ()) for d in PurePosixPath(p).parents}
    for path, line, dest, kind in entries:
        match = PINNED_URL.fullmatch(dest) if pinned is not None else None
        other = OTHER_PINNED_URL.fullmatch(dest) if others and not match else None
        if match:
            problem = check(path, match[2], pinned, pin_dirs if match[1] == 'tree' else set(),
                            pinned_anchors, roots=True)
            kind = 'pinned Markdown'
        elif other and other[2] in others:
            tree, tree_anchors = others[other[2]]
            tree_dirs = {str(d) for q in tree for d in PurePosixPath(q).parents} if other[1] == 'tree' else set()
            problem = check(path, other[3], tree, tree_dirs, tree_anchors, roots=True)
            kind = 'pinned Markdown'
        else:
            if not re.match(r'^(?:[A-Za-z][A-Za-z0-9+.-]*:|//)', dest):
                local += 1
            problem = check(path, dest, files, dirs, anchors)
        if problem:
            errors.append((path, line, dest, kind, problem))
    return errors, len(anchors), local


def paper_links(text, pinned, pinned_anchors, extra=None):
    r"""Check \proofsrc calls against PIN and each extra macro's calls against its own pin.

    ``extra`` maps a macro name to ``(pin, manifest, anchors)``.
    """
    errors = []
    count = 0
    text = re.sub(r'(?<!\\)(?:\\\\)*(%[^\n]*)', lambda m: m[0][:m.start(1) - m.start()], text)
    macros = {'proofsrc': (PIN, pinned, pinned_anchors), **(extra or {})}
    allowed = set()
    for name, (pin, files, anchors) in macros.items():
        definition = r'\\(?:newcommand|renewcommand|DeclareRobustCommand)\s*\{\s*\\' + name + r'\s*\}\s*\[2\]\s*\{'
        expected = re.escape(REPO_URL + 'blob/' + pin + '/#1')
        macro = re.findall(definition + r'\s*\\href\s*\{' + expected + r'\}\s*\{\s*#2\s*\}\s*\}', text)
        calls = 0
        for match in re.finditer(r'(?<!\\)(\\+)' + name + r'\b(?!\s*\})', text):
            line = text.count('\n', 0, match.start()) + 1
            if match[1] != '\\':
                errors.append((PAPER, line, '', 'pinned source', f'{name} must start with a single backslash'))
                continue
            args = re.match(r'\s*\{([^{}\s]+)\}\s*\{', text[match.end():])
            if not args:
                errors.append((PAPER, line, '', 'pinned source', 'malformed source link'))
                continue
            calls += 1
            dest = args[1].replace(r'\#', '#').replace(r'\_', '_')
            url = urlsplit(dest)
            if (url.scheme or url.netloc or url.query or not url.path or url.path.startswith('/')
                    or posixpath.normpath(url.path) != url.path or url.path.startswith('../')):
                problem = 'expected repository-relative pinned blob path'
            else:
                problem = check(PAPER, dest, files, set(), anchors, roots=True)
            if problem:
                errors.append((PAPER, line, dest, 'pinned source', problem))
        if calls and (len(macro) != 1 or len(re.findall(definition, text)) != 1):
            errors.append((PAPER, 1, '', 'pinned source', f'{name} must expand to the pinned GitHub blob URL'))
        if calls or name == 'proofsrc':
            allowed.add(REPO_URL + 'blob/' + pin + '/#1')
        count += calls
    if count:
        for m in re.finditer(r'\\href\s*\{(' + re.escape(REPO_URL) + r'[^{}]*)\}', text):
            if m[1] not in allowed:
                line = text.count('\n', 0, m.start()) + 1
                errors.append((PAPER, line, m[1], 'pinned source', 'repository sources must use proofsrc or ordersrc'))
    return errors, count


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--allow-existing', action='store_true', help='report pinned Markdown defects but fail only on regressions')
    ap.add_argument('--include', action='append', default=[], metavar='PATH', help='also check a newly written, unstaged Markdown document')
    ap.add_argument('--tex', type=Path, default=ROOT / PAPER, help='check links in another TeX input')
    args = ap.parse_args()
    files, pinned = manifest(), manifest(PIN)
    for p in args.include:
        candidate = (ROOT / p).resolve()
        if not candidate.is_relative_to(ROOT) or not candidate.is_file() or candidate.suffix.lower() != '.md':
            ap.error(f'expected an existing Markdown document inside the repository: {p}')
        files.add(candidate.relative_to(ROOT).as_posix())
    old_docs = dict(blobs(PIN, sorted(p for p in pinned if p.lower().endswith('.md'))))
    pinned_anchors = {p: parse(s)[0] for p, s in old_docs.items()}
    baseline, _, _ = scan(pinned, old_docs.items(), pinned, pinned_anchors)
    old = Counter((p, d, k, e) for p, _, d, k, e in baseline)
    paths = sorted(p for p in files if p.lower().endswith('.md') and (ROOT / p).is_file())
    working = {p for p in files if (ROOT / p).is_file()}
    others = {}
    for rev in (ORDER_PIN, ARCHIVE_PIN, BASE_PIN):
        tree = manifest(rev)
        docs_at = dict(blobs(rev, sorted(p for p in tree if p.lower().endswith('.md'))))
        others[rev] = (tree, {p: parse(s)[0] for p, s in docs_at.items()})
    issues, docs, local = scan(working, ((p, (ROOT / p).read_text(encoding='utf-8-sig')) for p in paths),
                               pinned, pinned_anchors, others)
    order_tree, order_anchors = others[ORDER_PIN]
    tex_issues, sources = paper_links(args.tex.read_text(encoding='utf-8'), pinned, pinned_anchors,
                                      {'ordersrc': (ORDER_PIN, order_tree, order_anchors)})
    if not sources:
        tex_issues.append((PAPER, 1, '', 'pinned source', 'no pinned proofsrc or ordersrc calls found'))
    prior = new = 0
    for p, n, d, k, e in issues + tex_issues:
        key = p, d, k, e
        if old[key]:
            prior += 1
            old[key] -= 1
            status = 'existing'
        else:
            new += 1
            status = 'NEW'
        print(f'{p}:{n}: {status} [{k}] {d!r}: {e}')
    print(f'{docs} Markdown docs; {local} local destinations; {sources} pinned TeX links; '
          f'{prior} existing and {new} new issues; {sum(old.values())} baseline issues no longer present')
    return int(bool(new or (prior and not args.allow_existing)))


if __name__ == '__main__':
    raise SystemExit(main())
