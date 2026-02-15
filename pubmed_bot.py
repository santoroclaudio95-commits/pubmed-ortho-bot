from Bio import Entrez
import requests
import re
from datetime import datetime, timedelta

# CONFIG
Entrez.email = "santoro.claudio95@gmail.com"
TELEGRAM_TOKEN = "{{TELEGRAM_TOKEN}}"
CHAT_ID = "{{CHAT_ID}}"

QUERY = """
(anterior cruciate ligament OR ACL OR ankle)
AND (orthopedic OR orthopaedic OR surgery OR reconstruction OR rehabilitation)
AND (randomized OR randomised OR meta-analysis)
AND english[Language]
"""

def has_significant_p(text):
    patterns = [
        r"p\s*<\s*0\.05",
        r"p\s*=\s*0\.0[0-4]",
        r"statistically significant"
    ]
    return any(re.search(p, text, re.I) for p in patterns)

def score_quality(text):
    score = 0
    if "meta-analysis" in text.lower():
        score += 10
    if "randomized" in text.lower():
        score += 8
    if "multicenter" in text.lower():
        score += 1
    if "double-blind" in text.lower():
        score += 1
    return score

def summarize(text):
    # riassunto semplice clinico
    sentences = re.split(r'(?<=[.!?]) +', text)
    return " ".join(sentences[:2])

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

        pmid = art["MedlineCitation"]["PMID"]
        link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}"
        summary = summarize(abstract)

        msg = (
            f"🦵 ACL / ANKLE — High-Quality Evidence\n\n"
            f"{title}\n"
            f"Quality: {quality}/10\n\n"
            f"{summary}\n\n"
            f"{link}"
        )
        messages.append(msg)

if messages:
    text = "\n\n".join(messages)
else:
    text = "🦵 Nessun articolo nuovo oggi su ACL o caviglia che soddisfi i criteri di qualità."
    url = f"https://api.telegram.org/bot{TELEGRAM_TOKEN}/sendMessage"
    requests.post(url, data={"chat_id": CHAT_ID, "text": text})
    print("Invio messaggio Telegram...")
    print(requests.json())  # mostrerà il JSON di risposta di Telegram
