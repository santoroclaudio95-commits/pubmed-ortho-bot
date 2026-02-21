import os
from Bio import Entrez
import requests
import re
from datetime import datetime, timedelta

# ---------------- CONFIG ----------------
Entrez.email = "santoro.claudio95@gmail.com"   # la tua email reale
TELEGRAM_TOKEN = os.environ.get("TELEGRAM_TOKEN")
CHAT_ID = os.environ.get("CHAT_ID")

# ---------------- QUERY -----------------
QUERY = """
(
  anterior cruciate ligament OR ACL
  OR meniscus OR meniscal
  OR ankle
  OR foot
  OR achilles tendon OR achilles
  OR hip
)
AND
(
  orthopedic OR orthopaedic
  OR surgery OR surgical
  OR reconstruction
  OR repair
  OR rehabilitation
)
AND
(
  randomized OR randomised
  OR meta-analysis
  OR systematic review
)
AND english[Language]
OR
(
  sport OR sports OR athlete OR athletes OR "return to sport"
)
"""

# ---------- Funzioni utili --------------
def has_significant_p(text):
    patterns = [
        r"p\s*<\s*0\.05",
        r"p\s*=\s*0\.0[0-4]",
        r"statistically significant"
    ]
    return any(re.search(p, text, re.I) for p in patterns)

def score_quality(text):
    score = 0
    if "meta-analysis" in text.lower() or "systematic review" in text.lower():
        score += 10
    if "randomized" in text.lower() or "randomised" in text.lower():
        score += 8
    if "multicenter" in text.lower():
        score += 1
    if "double-blind" in text.lower():
        score += 1
    return score

def extract_p_value(text):
    match = re.search(r"p\s*[=≤<]\s*0\.\d+", text, re.I)
    if match:
        return match.group(0)
    return "p non riportato"

def extract_ci(text):
    match = re.search(r"95%?\s*CI[:\s]*[\(\[]?\s*\d+\.?\d*\s*[-–]\s*\d+\.?\d*", text, re.I)
    if match:
        return match.group(0)
    return "CI non riportato"

def summarize(text):
    sentences = re.split(r'(?<=[.!?]) +', text)
    return " ".join(sentences[:2])  # primi 2 periodi come riassunto clinico

def tag_article(text):
    tags = []
    text_lower = text.lower()
    if any(k in text_lower for k in ["acl", "anterior cruciate ligament", "meniscus", "meniscal"]):
        tags.append("🦵 Knee")
    if any(k in text_lower for k in ["ankle", "foot", "achilles"]):
        tags.append("🦶 Foot/Ankle")
    if "hip" in text_lower:
        tags.append("🦴 Hip")
    if any(k in text_lower for k in ["sport", "athlete", "return to sport"]):
        tags.append("🏃 Sport")
    return " ".join(tags) if tags else ""

# ------------- Cerca articoli -------------
today = datetime.today()
yesterday = today - timedelta(days=1)

handle = Entrez.esearch(
    db="pubmed",
    term=QUERY,
    mindate=yesterday.strftime("%Y/%m/%d"),
    maxdate=today.strftime("%Y/%m/%d"),
    datetype="pdat",
    retmax=20
)

ids = Entrez.read(handle)["IdList"]
messages = []

if ids:
    fetch = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    articles = Entrez.read(fetch)

    for art in articles["PubmedArticle"]:
        article = art["MedlineCitation"]["Article"]
        title = str(article["ArticleTitle"])
        abstract = " ".join(str(x) for x in article.get("Abstract", {}).get("AbstractText", []))

        if not abstract:
            continue
        if not has_significant_p(abstract):
            continue

        quality = score_quality(abstract)
        if quality < 8:
            continue

        p_value = extract_p_value(abstract)
        ci_value = extract_ci(abstract)

        pmid = art["MedlineCitation"]["PMID"]
        link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}"
        summary = summarize(abstract)
        tags = tag_article(title + " " + abstract)

        msg = (
            f"{tags} — High-Quality Evidence\n\n"
            f"{title}\n\n"
            f"P-value: {p_value}\n"
            f"CI: {ci_value}\n"
            f"Quality: {quality}/10\n\n"
            f"{summary}\n\n"
            f"{link}"
        )
        messages.append(msg)

# --------- Messaggio Telegram ----------
if messages:
    text = "\n\n".join(messages)
else:
    text = "🦵 Nessun articolo nuovo oggi su ACL, menisco, caviglia, piede, achilles o anca che soddisfi i criteri di qualità."

url = f"https://api.telegram.org/bot{TELEGRAM_TOKEN}/sendMessage"
r = requests.post(url, data={"chat_id": CHAT_ID, "text": text})

print("Invio messaggio Telegram...")
print("Risposta Telegram:", r.json())
