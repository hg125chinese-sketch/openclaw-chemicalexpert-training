---
name: chem-literature
description: Search, retrieve, and critically read chemistry / AI4Chem papers from arXiv, Semantic Scholar, and PubChem. Produce structured notes for the research knowledge base.
homepage: https://arxiv.org
metadata: { "openclaw": { "emoji": "📑", "requires": { "bins": ["curl", "python3"] } } }
---

# Chemical Literature Search & Critical Reading

Find, retrieve, and deeply analyze chemistry and AI-for-chemistry papers, then produce structured notes that integrate into the `research/` knowledge base.

## When to Use

- User asks to find papers on a chemistry / ML-for-chemistry topic
- User asks to read, summarize, or critically analyze a paper
- User provides an arXiv ID, DOI, or title and wants a deep-read note
- User asks to update the reading list or literature survey

## Tools Available

### 1. arXiv API (free, no key)

Search:
```bash
curl -s "http://export.arxiv.org/api/query?search_query=all:molecule+generation+vae&start=0&max_results=10&sortBy=submittedDate&sortOrder=descending" | python3 -c "
import sys, xml.etree.ElementTree as ET
ns = {'a': 'http://www.w3.org/2005/Atom'}
root = ET.parse(sys.stdin).getroot()
for i, entry in enumerate(root.findall('a:entry', ns)):
    title = entry.find('a:title', ns).text.strip().replace('\n',' ')
    arxiv_id = entry.find('a:id', ns).text.strip().split('/')[-1]
    published = entry.find('a:published', ns).text[:10]
    authors = ', '.join(a.find('a:name', ns).text for a in entry.findall('a:author', ns)[:3])
    cats = entry.find('a:category', ns).attrib.get('term','')
    print(f'{i+1}. [{arxiv_id}] {published} | {authors} | {title} | {cats}')
"
```

Fetch single paper metadata:
```bash
curl -s "http://export.arxiv.org/api/query?id_list=2312.01234" | python3 -c "
import sys, xml.etree.ElementTree as ET
ns = {'a': 'http://www.w3.org/2005/Atom'}
entry = ET.parse(sys.stdin).getroot().find('a:entry', ns)
print('Title:', entry.find('a:title', ns).text.strip().replace('\n',' '))
print('Authors:', ', '.join(a.find('a:name', ns).text for a in entry.findall('a:author', ns)))
print('Abstract:', entry.find('a:summary', ns).text.strip()[:500])
print('Categories:', ' '.join(c.attrib['term'] for c in entry.findall('a:category', ns)))
print('Published:', entry.find('a:published', ns).text[:10])
print('PDF:', entry.find('a:id', ns).text.strip().replace('abs','pdf'))
"
```

Download PDF:
```bash
curl -L -o paper.pdf "https://arxiv.org/pdf/2312.01234"
```

### 2. Semantic Scholar API (free, no key, richer metadata)

Search:
```bash
curl -s "https://api.semanticscholar.org/graph/v1/paper/search?query=molecule+generation+diffusion&limit=10&fields=title,authors,year,citationCount,externalIds,abstract" | python3 -c "
import sys, json
data = json.load(sys.stdin)
for i, p in enumerate(data.get('data', [])):
    authors = ', '.join(a['name'] for a in (p.get('authors') or [])[:3])
    arxiv = (p.get('externalIds') or {}).get('ArXiv', '')
    print(f\"{i+1}. [{p.get('year','')}] {p['title']} | {authors} | cites:{p.get('citationCount',0)} | arXiv:{arxiv}\")
"
```

Get citation graph (who cites whom):
```bash
# citations of a paper
curl -s "https://api.semanticscholar.org/graph/v1/paper/ArXiv:1610.02415/citations?fields=title,year,citationCount&limit=20"
# references of a paper
curl -s "https://api.semanticscholar.org/graph/v1/paper/ArXiv:1610.02415/references?fields=title,year,citationCount&limit=20"
```

### 3. PubChem (compound/assay lookup)

```bash
# Search compound by name
curl -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/property/MolecularFormula,MolecularWeight,IsomericSMILES,IUPACName/JSON" | python3 -m json.tool

# Search by SMILES
curl -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/CC(=O)Oc1ccccc1C(=O)O/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
```

### 4. Browser (for paywalled or complex pages)

If a paper is not on arXiv or a full-text needs browser rendering:
```bash
npx agent-browser open "https://doi.org/10.1234/example"
npx agent-browser snapshot -i
npx agent-browser screenshot --full "paper-page.png"
npx agent-browser close
```

## Rate Limits & Etiquette

- **arXiv**: Max 1 request/3 seconds. Batch queries, don't hammer.
- **Semantic Scholar**: 100 requests/5 min without key. Add `x-api-key` header if you have one.
- **PubChem**: 5 requests/second. Be polite.

## Critical Reading Protocol

When doing a deep read of a paper, produce a note with ALL of these sections:

### Paper Note Template

Save to: `research/ai4chem/papers/<topic>/<NNN>-<first-author>-<year>-<short-slug>.md`

```markdown
# <Full Title>

| Field       | Value |
|-------------|-------|
| Authors     | ... |
| Year        | ... |
| Venue       | ... |
| arXiv       | ... |
| DOI         | ... |
| Code        | ... (repo URL if available) |
| Cites       | ... (Semantic Scholar count) |

## TL;DR (1-2 sentences)

## Problem & Motivation
What gap does this fill? Why now?

## Method
- Architecture / algorithm
- Key design decisions and why
- Training details (data, loss, optimizer, epochs, hardware)

## Results
- Main claims + evidence (tables/figures referenced)
- Baselines compared against
- Ablations (what was tested)

## Strengths
- (numbered list)

## Weaknesses & Red Flags
- (numbered list)
- Data leakage risk?
- Cherry-picked results?
- Missing baselines?
- Reproducibility issues?

## Key Equations / Figures
- Equation N: ... (what it means in plain language)

## Connections
- Builds on: [paper refs]
- Superseded by: [if known]
- Related in our reading list: [cross-refs to our notes]

## Reproducibility Assessment
- Code available? Runs?
- Data available? License?
- Can we replicate Table X with our setup?

## Open Questions
- Things to investigate or verify

## Raw Notes
- (free-form observations during reading)
```

## Workflow: Full Literature Search

When the user asks "find me papers about X":

1. **Search** arXiv + Semantic Scholar (2 queries each, different angles)
2. **Deduplicate** by arXiv ID or title similarity
3. **Rank** by: relevance to query > citation count > recency
4. **Present** top 5-10 as a table: `| # | Year | Title | Authors | Cites | arXiv | Why relevant |`
5. **Ask** user which to deep-read, or auto-select top 3 if user said "go deep"
6. **Deep-read** each selected paper using the template above
7. **Update** `research/ai4chem/papers/<topic>/reading-list.md` with new entries
8. **Cross-reference** with existing notes (check `research/ai4chem/index.md`)

## Workflow: Single Paper Deep Read

When given a specific paper (arXiv ID, DOI, or title):

1. **Fetch metadata** from arXiv + Semantic Scholar
2. **Download PDF** if available
3. **Read** the paper (use PDF text extraction or browser if needed)
4. **Write** the structured note using the template
5. **Save** to the correct location in `research/`
6. **Update** reading-list.md and index.md
7. **Report** TL;DR + strengths + weaknesses + open questions to user

## Workflow: Literature Survey / Comparison

When asked to compare papers or survey a subfield:

1. Gather all relevant notes from `research/ai4chem/papers/<topic>/`
2. Build a comparison table: `| Paper | Method | Data | Main Result | Weakness |`
3. Write a synthesis note at `research/ai4chem/notes/<topic>/<slug>.md`
4. Identify gaps: what hasn't been tried? What claims conflict?
5. Propose next steps or experiments

## Python Helper: Extract Text from PDF

If `pdftotext` is unavailable, use Python:
```python
#!/opt/conda/envs/chem/bin/python
"""Extract text from a PDF using PyMuPDF (fitz) if available, else basic method."""
import sys

def extract_with_fitz(path):
    import fitz
    doc = fitz.open(path)
    return '\n'.join(page.get_text() for page in doc)

def extract_basic(path):
    """Fallback: read raw bytes and extract printable text."""
    with open(path, 'rb') as f:
        raw = f.read()
    import re
    # Extract text between BT and ET operators (very rough)
    texts = re.findall(rb'\((.*?)\)', raw)
    return '\n'.join(t.decode('latin-1', errors='ignore') for t in texts[:200])

if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'paper.pdf'
    try:
        print(extract_with_fitz(path))
    except ImportError:
        print("[WARN] PyMuPDF not installed, using basic extraction")
        print(extract_basic(path))
```

Install PyMuPDF if needed (run once):
```bash
/opt/conda/envs/chem/bin/pip install PyMuPDF
```

## Integration with Existing Knowledge Base

- **Index**: Always update `research/ai4chem/index.md` after adding new material
- **Cross-links**: Use relative links between notes: `[see VAE training tricks](../notes/molecule-generation/vae-training-tricks-for-strings.md)`
- **Numbering**: Papers are numbered sequentially per topic: `001-`, `002-`, etc. Check existing files before assigning the next number.
- **Git**: After writing notes, `cd /home/node/.openclaw/workspace-chemicalexpert && git add -A && git commit -m "lit: add <short description>"`

## Checklist Before Finishing

- [ ] Note saved in correct location with correct numbering
- [ ] reading-list.md updated
- [ ] index.md updated (if new topic or significant addition)
- [ ] Git committed
- [ ] User given: TL;DR + key finding + any red flags
