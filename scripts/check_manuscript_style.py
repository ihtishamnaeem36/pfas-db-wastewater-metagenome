import docx
import re
from validate_human_style import analyze_paragraph

doc = docx.Document('../PFAS_Candidate_Gene_Manuscript.docx')
total_paras = 0
total_sentences = 0
short_sentences = 0
all_issues = []

for i, p in enumerate(doc.paragraphs):
    text = p.text.strip()
    if not text or len(text.split()) < 5:
        continue
    # skip references list entries
    if re.match(r'^[A-Z][a-zA-Z\s\.,\-]+,\s*\d{4}\.', text) or text.startswith('Methodological references:') or text.startswith('PFAS/dehalogenase-biochemistry references:'):
        continue
    total_paras += 1
    issues = analyze_paragraph(text, section_name=text[:30])
    if issues:
        all_issues.append((i, text[:50], issues))
    sents = [s for s in re.split(r'(?<=[.!?])\s+', text) if s.strip()]
    total_sentences += len(sents)
    short_sentences += sum(1 for s in sents if len(s.split()) < 15)

print(f"Total analyzed body paragraphs: {total_paras}")
print(f"Total sentences: {total_sentences}")
print(f"Short sentences (<15 words): {short_sentences} ({short_sentences/total_sentences*100:.1f}%)")
print(f"Paragraphs with flagged issues: {len(all_issues)}")
for idx, snippet, iss in all_issues:
    print(f"\n  Para {idx} (\"{snippet}...\"):")
    for item in iss:
        print(f"    - {item}")
