"""
Validator script to enforce the Human Academic Writing Style Guide rules:
1. Sentence-length rhythm: at least 30% of sentences < 15 words; jagged varied lengths.
2. Zero banned AI vocabulary words.
3. Zero copula-avoidance verbs (represents, constitutes, occupies, serves as, reflects).
4. Zero rhetorical em-dashes or double-dashes (excluding CLI arguments).
5. Zero trailing '-ing' significance clauses.
6. Zero misuse of 'enriched' (only permitted as 'enrichment culture').
7. British spelling compliance.
"""
import re
import sys
from pathlib import Path

BANNED_WORDS = [
    r'\badditionally\b', r'\bmoreover\b', r'\bfurthermore\b', r'\bnotably\b',
    r'\bimportantly\b', r'\bcrucial\b', r'\bcrucially\b', r'\bpivotal\b',
    r'\blandscape\b', r'\btapestry\b', r'\brealm\b',
    r'\bunderscore[s|d|ing]?\b',
    r'\bhighlight[s|ed|ing]?\b',  # checked separately if not in heading
    r'\bshowcase[s|d|ing]?\b',
    r'\bdelve[s|d|ing]?\b', r'\bfoster[s|d|ing]?\b', r'\bgarner[s|ed|ing]?\b',
    r'\bintricate\b', r'\bvibrant\b', r'\btestament\b',
    r'\bworth noting\b', r'\bimportant to note\b',
    r'\bexceptional\b',
]

COPULA_AVOIDANCE = [
    r'\brepresents\b', r'\bconstitutes\b', r'\boccupies\b', r'\bserves as\b',
    r'\breflecting\b',
]

US_TO_BRITISH = {
    'characterize': 'characterise',
    'characterized': 'characterised',
    'characterizing': 'characterising',
    'prioritization': 'prioritisation',
    'prioritize': 'prioritise',
    'program': 'programme', # when referring to surveillance program
    'programs': 'programmes',
    'catalog': 'catalogue',
    'catalogs': 'catalogues',
    'homolog': 'homologue',
    'homologs': 'homologues',
    'modeled': 'modelled',
    'modeling': 'modelling',
    'neighbor': 'neighbour',
    'neighborhood': 'neighbourhood',
    'neighbors': 'neighbours',
    'color': 'colour',
    'colored': 'coloured',
}

def analyze_paragraph(text, section_name=""):
    # Split text into sentences
    sentences = [s.strip() for s in re.split(r'(?<=[.!?])\s+', text) if s.strip()]
    if not sentences:
        return []
    
    issues = []
    
    # Sentence length statistics
    lengths = [len(s.split()) for s in sentences]
    short_count = sum(1 for l in lengths if l < 15)
    pct_short = (short_count / len(lengths)) * 100 if lengths else 0
    
    if len(sentences) >= 3 and pct_short < 30:
        issues.append(f"Low short-sentence ratio: {short_count}/{len(lengths)} ({pct_short:.1f}%) < 15 words.")
    
    # Check 3 consecutive sentences of similar length (within 3 words)
    for i in range(len(lengths) - 2):
        if abs(lengths[i] - lengths[i+1]) <= 2 and abs(lengths[i+1] - lengths[i+2]) <= 2 and lengths[i] > 18:
            issues.append(f"Monotonous rhythm at sentence {i+1}-{i+3}: lengths {lengths[i]}, {lengths[i+1]}, {lengths[i+2]}.")
            
    # Check long sentence followed by long sentence (> 28 words)
    for i in range(len(lengths) - 1):
        if lengths[i] > 28 and lengths[i+1] > 24:
            issues.append(f"Back-to-back long sentences at {i+1}-{i+2}: lengths {lengths[i]} and {lengths[i+1]}.")

    # Check banned words
    for pat in BANNED_WORDS:
        m = re.findall(pat, text, re.IGNORECASE)
        if m:
            # allow 'highlight' if in Section heading Highlights
            if 'highlight' in pat and 'Highlights' in section_name:
                continue
            issues.append(f"Banned word found: {m}")

    # Check copula avoidance
    for pat in COPULA_AVOIDANCE:
        m = re.findall(pat, text, re.IGNORECASE)
        if m:
            issues.append(f"Copula-avoidance verb found: {m}")

    # Check em-dashes / double-dashes (ignore CLI flags like --ultra-sensitive, --auto)
    dashes = re.findall(r'(?<![a-zA-Z0-9\-])--(?![\w\-])', text)
    if dashes:
        issues.append(f"Rhetorical double-dash '--' found: {len(dashes)} instances.")
    if '—' in text or '–' in text:
        # Check if en-dash used as rhetorical dash rather than number range
        rhetorical_en = re.findall(r'\s[–—]\s', text)
        if rhetorical_en:
            issues.append(f"Rhetorical em/en-dash found: {len(rhetorical_en)} instances.")

    # Check 'enriched'
    enriched_matches = re.finditer(r'\benriched\b', text, re.IGNORECASE)
    for em in enriched_matches:
        sub = text[max(0, em.start()-15):min(len(text), em.end()+15)]
        if 'culture' not in sub.lower():
            issues.append(f"Suspicious 'enriched' without 'culture': '{sub}'")

    # Check US spelling
    for us, br in US_TO_BRITISH.items():
        m = re.findall(r'\b' + us + r'\b', text, re.IGNORECASE)
        if m:
            issues.append(f"US spelling '{us}' should be British '{br}': {len(m)} instances.")

    # Check trailing -ing clauses (e.g. ', highlighting', ', reflecting', ', underscoring', ', ensuring')
    trailing_ing = re.findall(r',\s+(?:highlighting|reflecting|underscoring|emphasizing|emphasising|creating|ensuring)\b', text, re.IGNORECASE)
    if trailing_ing:
        issues.append(f"Superficial trailing -ing clause found: {trailing_ing}")

    return issues

if __name__ == "__main__":
    print("Validator module ready.")
