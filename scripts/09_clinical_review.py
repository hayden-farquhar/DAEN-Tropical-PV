"""
07 — Clinical Face-Validity Review (Phase 7, Section 9.7)

Blinded pharmacological plausibility assessment for all 68 signals.
Rating is based solely on pharmacological and clinical knowledge,
BEFORE consulting literature or comparator databases.

3-point scale:
  - Established: well-known from pharmacology or clinical literature
  - Plausible novel: biologically plausible but not well-documented
  - Noise: lacks biological plausibility; likely artefact or confounding

Outputs:
  - outputs/tables/clinical_review_log.csv   Blinded ratings + rationale
  - outputs/tables/clinical_review_summary.txt

Usage:
    python scripts/09_clinical_review.py
"""

from pathlib import Path

import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "outputs" / "tables"

# Pre-specified blinded ratings based on pharmacological knowledge.
# These are rated BEFORE any literature search or comparator lookup.
BLINDED_RATINGS = {
    # ── Brown snake antivenom ──────────────────────────────────────────
    ("Brown snake antivenom", "Hypofibrinogenaemia"): (
        "Noise",
        "Indication confounding: hypofibrinogenaemia is a feature of VICC "
        "(venom-induced consumptive coagulopathy), the condition being treated. "
        "The PT reflects the envenomation, not an antivenom adverse effect.",
    ),
    ("Brown snake antivenom", "Underdose"): (
        "Noise",
        "Medication error PT, not a pharmacological adverse effect. Reflects "
        "dosing practice (antivenom dose titrated to coagulopathy reversal).",
    ),
    ("Brown snake antivenom", "Serum sickness"): (
        "Established",
        "Type III hypersensitivity to equine F(ab')2 immunoglobulin fragments. "
        "Well-characterised delayed reaction (5-14 days) to heterologous serum "
        "products. Expected class effect for all equine-derived antivenoms.",
    ),
    ("Brown snake antivenom", "Cerebral haemorrhage"): (
        "Plausible novel",
        "Could represent: (a) complication of brown snake VICC with delayed "
        "coagulopathy reversal despite antivenom; (b) secondary haemorrhagic "
        "event during the peri-treatment period. Biologically plausible via "
        "residual coagulopathy pathway, but not a direct antivenom effect.",
    ),
    # ── Box jellyfish antivenom ────────────────────────────────────────
    ("Box jellyfish antivenom", "Cardio-respiratory arrest"): (
        "Plausible novel",
        "Dual mechanism: (a) Chironex fleckeri envenomation causes direct "
        "cardiotoxicity (QT prolongation, cardiovascular collapse); (b) "
        "antivenom anaphylaxis could also produce cardio-respiratory arrest. "
        "Cannot distinguish drug effect from indication effect in PV data. "
        "Small n (3) in a rare product.",
    ),
    # ── Tiger snake antivenom ──────────────────────────────────────────
    ("Tiger snake antivenom", "Type I hypersensitivity"): (
        "Established",
        "IgE-mediated immediate hypersensitivity to equine-derived protein. "
        "Well-established class effect for all heterologous antivenoms.",
    ),
    ("Tiger snake antivenom", "Anaphylactic reaction"): (
        "Established",
        "Anaphylaxis to equine immunoglobulin. The most well-characterised "
        "adverse effect of antivenom therapy. Occurs in 5-25% of antivenom "
        "administrations per Australian prospective data.",
    ),
    ("Tiger snake antivenom", "Erythema"): (
        "Established",
        "Cutaneous manifestation of immediate hypersensitivity. Component of "
        "the anaphylaxis/urticaria/angioedema spectrum for equine antivenoms.",
    ),
    ("Tiger snake antivenom", "Hypotension"): (
        "Established",
        "Haemodynamic consequence of anaphylaxis (distributive shock from "
        "histamine-mediated vasodilation). Expected in the antivenom "
        "hypersensitivity reaction spectrum.",
    ),
    ("Tiger snake antivenom", "Urticaria"): (
        "Established",
        "IgE-mediated cutaneous reaction to equine protein. Component of "
        "the type I hypersensitivity spectrum for antivenoms.",
    ),
    # ── Polyvalent snake antivenom ─────────────────────────────────────
    ("Polyvalent snake antivenom", "Anaphylactic reaction"): (
        "Established",
        "Higher protein load (five-species polyvalent) increases anaphylaxis "
        "risk vs monovalent. Well-established; positive control for this study.",
    ),
    # ── Primaquine ─────────────────────────────────────────────────────
    ("Primaquine", "Methaemoglobinaemia"): (
        "Established",
        "Direct pharmacological effect via oxidant metabolites (primarily "
        "5-hydroxyprimaquine). Dose-dependent; markedly potentiated by G6PD "
        "deficiency. H5 positive control.",
    ),
    # ── JEV vaccines ──────────────────────────────────────────────────
    ("JEV vaccines", "Contraindication to vaccination"): (
        "Noise",
        "Administrative PT, not a pharmacological adverse effect. Reflects "
        "vaccination programme quality (vaccine given despite contraindication).",
    ),
    ("JEV vaccines", "Incorrect route of product administration"): (
        "Noise",
        "Medication error PT. Reflects mass vaccination programme logistics, "
        "not vaccine pharmacology. Surge in 2022 outbreak rollout.",
    ),
    ("JEV vaccines", "Oedema mouth"): (
        "Plausible novel",
        "Angioedema variant from vaccine hypersensitivity. Biologically "
        "plausible (IgE-mediated mucosal oedema) but uncommon for inactivated "
        "JEV vaccine. Low n (3).",
    ),
    ("JEV vaccines", "Product prescribing issue"): (
        "Noise",
        "Administrative/error PT, not a pharmacological effect.",
    ),
    ("JEV vaccines", "Injection site inflammation"): (
        "Established",
        "Local reactogenicity: innate immune activation at injection site. "
        "Expected for all injectable vaccines. Listed in product information.",
    ),
    ("JEV vaccines", "Exposure during pregnancy"): (
        "Noise",
        "Risk situation PT, not an adverse effect. Reflects inappropriate "
        "administration (JEV vaccine category B3 in pregnancy).",
    ),
    ("JEV vaccines", "Photophobia"): (
        "Plausible novel",
        "Neurological symptom. Could represent: (a) post-vaccination headache "
        "with photophobia (common); (b) rare aseptic meningitis. Biologically "
        "plausible given known neurotropism of JEV. Low n (6).",
    ),
    ("JEV vaccines", "Product administered to patient of inappropriate age"): (
        "Noise",
        "Medication error PT. Reflects vaccination programme logistics during "
        "emergency rollout (age-group eligibility confusion).",
    ),
    ("JEV vaccines", "Vaccination error"): (
        "Noise",
        "Generic medication error PT. Dominated by 2022 outbreak rollout "
        "administrative errors.",
    ),
    ("JEV vaccines", "Urticaria"): (
        "Established",
        "Type I hypersensitivity to vaccine components (stabilisers, "
        "protamine sulphate residual). Listed in product information.",
    ),
    # ── Ivermectin ─────────────────────────────────────────────────────
    ("Ivermectin", "Hypervitaminosis B6"): (
        "Noise",
        "COVID-era confounding: co-administration with high-dose vitamin B6 "
        "supplements during off-label COVID-19 use. Not a pharmacological "
        "effect of ivermectin. TGA ivermectin advisory September 2021.",
    ),
    ("Ivermectin", "Vitamin B6 increased"): (
        "Noise",
        "Same COVID-era confounding as hypervitaminosis B6. Concomitant "
        "supplement use, not ivermectin pharmacology.",
    ),
    ("Ivermectin", "Abdominal distension"): (
        "Established",
        "GI effect of ivermectin. Abdominal symptoms (distension, nausea, "
        "diarrhoea) are listed adverse effects, particularly at higher doses "
        "or during treatment of heavy parasite burden (Mazzotti-type reaction "
        "from parasite die-off).",
    ),
    # ── Doxycycline ───────────────────────────────────────────────────
    ("Doxycycline", "Photosensitivity reaction"): (
        "Established",
        "Tetracycline class phototoxicity via UV-A-activated oxygen radicals. "
        "Dose-dependent. Listed CMI warning. Study positive control.",
    ),
    ("Doxycycline", "Oesophagitis"): (
        "Established",
        "Direct mucosal chemical irritation from capsule dissolution in the "
        "oesophagus. Well-characterised 'pill oesophagitis'. Prevented by "
        "taking with a full glass of water and remaining upright.",
    ),
    ("Doxycycline", "Oesophagitis ulcerative"): (
        "Established",
        "Severe end of pill oesophagitis spectrum. Same mechanism as "
        "oesophagitis: direct chemical mucosal injury.",
    ),
    ("Doxycycline", "Oesophageal ulcer"): (
        "Established",
        "Focal mucosal ulceration from pill contact injury. Same spectrum as "
        "oesophagitis/oesophagitis ulcerative.",
    ),
    ("Doxycycline", "Dysphagia"): (
        "Established",
        "Symptom of oesophageal irritation/pill oesophagitis. Same mechanism.",
    ),
    ("Doxycycline", "Jarisch-Herxheimer reaction"): (
        "Established",
        "Cytokine storm from spirochetal endotoxin release during treatment "
        "of Lyme disease, syphilis, leptospirosis. Not a drug-specific effect "
        "but a pathogen-treatment interaction well-known for doxycycline.",
    ),
    ("Doxycycline", "Intracranial pressure increased"): (
        "Established",
        "Pseudotumour cerebri (benign intracranial hypertension). Tetracycline "
        "class effect. Risk increased with concomitant retinoids.",
    ),
    ("Doxycycline", "Papilloedema"): (
        "Established",
        "Sign of raised intracranial pressure (same mechanism as above). "
        "Downstream manifestation of pseudotumour cerebri.",
    ),
    ("Doxycycline", "Drug-induced liver injury"): (
        "Established",
        "Hepatotoxicity is a recognised tetracycline class effect. "
        "Doxycycline-specific DILI is uncommon but documented in LiverTox.",
    ),
    ("Doxycycline", "Hepatitis"): (
        "Established",
        "Hepatocellular injury pattern. Same class-effect mechanism.",
    ),
    ("Doxycycline", "Hepatitis cholestatic"): (
        "Established",
        "Cholestatic or mixed pattern DILI. Documented for doxycycline.",
    ),
    ("Doxycycline", "Jaundice"): (
        "Established",
        "Clinical sign of hepatotoxicity or cholestasis. Downstream "
        "manifestation of doxycycline DILI.",
    ),
    ("Doxycycline", "Fixed eruption"): (
        "Established",
        "Doxycycline is one of the most common causes of fixed drug eruption. "
        "Cell-mediated hypersensitivity with tissue-resident memory T cells.",
    ),
    ("Doxycycline", "Genital ulceration"): (
        "Established",
        "Fixed drug eruption with genital predilection. Doxycycline is "
        "specifically known for genital fixed drug eruption.",
    ),
    ("Doxycycline", "Acute generalised exanthematous pustulosis"): (
        "Plausible novel",
        "AGEP is a T-cell-mediated severe cutaneous adverse reaction. "
        "Doxycycline is not among the classic AGEP triggers (aminopenicillins, "
        "macrolides) but case reports exist. Biologically plausible via "
        "delayed-type hypersensitivity.",
    ),
    ("Doxycycline", "Toxic epidermal necrolysis"): (
        "Plausible novel",
        "TEN is HLA-mediated severe drug hypersensitivity. Tetracyclines "
        "are not classic TEN triggers but case reports exist. The signal "
        "may reflect confounding (concomitant drugs) or rare genuine "
        "hypersensitivity.",
    ),
    ("Doxycycline", "Dermatitis bullous"): (
        "Plausible novel",
        "Bullous skin reaction. Could represent phototoxic bullous dermatosis "
        "(extension of phototoxicity) or immunobullous reaction. Biologically "
        "plausible via phototoxic pathway.",
    ),
    ("Doxycycline", "Face oedema"): (
        "Established",
        "Angioedema from drug hypersensitivity. Documented for doxycycline; "
        "can occur as part of DRESS or isolated angioedema.",
    ),
    ("Doxycycline", "Nail disorder"): (
        "Established",
        "Photo-onycholysis: UV-A-mediated detachment of nail plate from nail "
        "bed. Well-characterised tetracycline class effect, same phototoxic "
        "mechanism as photosensitivity.",
    ),
    ("Doxycycline", "Pathogen resistance"): (
        "Established",
        "Antimicrobial resistance selection. Expected for any antibiotic with "
        "widespread use. Not a pharmacological ADR per se but a recognised "
        "consequence of tetracycline prescribing.",
    ),
    ("Doxycycline", "Candida infection"): (
        "Established",
        "Opportunistic fungal infection from antibiotic-induced disruption "
        "of commensal flora. Expected class effect for broad-spectrum "
        "antibiotics.",
    ),
    ("Doxycycline", "Aspergillus infection"): (
        "Plausible novel",
        "Opportunistic fungal infection. More typically associated with "
        "immunosuppression than antibiotics alone. May reflect confounding "
        "(doxycycline in immunocompromised patients) rather than direct "
        "causation. Low plausibility for direct effect.",
    ),
    ("Doxycycline", "Pyroglutamic acidosis"): (
        "Noise",
        "Pyroglutamic acidosis (5-oxoprolinuria) is associated with "
        "paracetamol/acetaminophen in malnourished patients, not doxycycline. "
        "Likely confounding from concomitant paracetamol use.",
    ),
    ("Doxycycline", "Hidradenitis"): (
        "Noise",
        "Doxycycline is a treatment for hidradenitis suppurativa, not a "
        "cause. This signal likely reflects indication bias (the condition "
        "being treated, not an adverse effect).",
    ),
    ("Doxycycline", "Completed suicide"): (
        "Noise",
        "No pharmacological mechanism for doxycycline-induced suicidality. "
        "Likely confounding: doxycycline prescribed for acne/skin conditions "
        "in demographics with higher baseline mental health burden, or "
        "concomitant isotretinoin use.",
    ),
    ("Doxycycline", "Drug abuse"): (
        "Noise",
        "Doxycycline has no abuse potential (no CNS reward pathway activity). "
        "Signal artefact from: (a) confounding with co-prescribed "
        "drugs in substance-use populations; (b) coding of off-label use.",
    ),
    ("Doxycycline", "Intentional overdose"): (
        "Noise",
        "Not a pharmacological ADR. Reflects the use of available medications "
        "in self-harm. No drug-specific signal.",
    ),
    ("Doxycycline", "Somnambulism"): (
        "Plausible novel",
        "Parasomnia. No established mechanism for doxycycline-induced "
        "somnambulism, but CNS penetration occurs. Could also reflect "
        "confounding with concomitant CNS-active drugs. Weak plausibility.",
    ),
    ("Doxycycline", "Abnormal sleep-related event"): (
        "Plausible novel",
        "Broader sleep disruption PT. Same assessment as somnambulism: "
        "weak but not impossible given CNS penetration. May cluster with "
        "somnambulism signal.",
    ),
    ("Doxycycline", "Amnesia"): (
        "Noise",
        "No established mechanism for tetracycline-induced amnesia. "
        "Likely confounding (concomitant benzodiazepines, zolpidem) or "
        "co-reported with the suicidality/overdose cluster.",
    ),
    ("Doxycycline", "Parosmia"): (
        "Noise",
        "Olfactory distortion. In the temporal context of this dataset, "
        "parosmia is likely COVID-19-related (prominent post-COVID symptom), "
        "confounding doxycycline reports from 2020-2022.",
    ),
    ("Doxycycline", "Cytopenia"): (
        "Plausible novel",
        "Haematological suppression. Rare case reports of doxycycline-associated "
        "thrombocytopenia and neutropenia exist. Mechanism unclear; possibly "
        "immune-mediated. Low frequency but biologically plausible.",
    ),
    ("Doxycycline", "Cauda equina syndrome"): (
        "Noise",
        "No pharmacological mechanism. Doxycycline does not affect spinal "
        "cord or nerve roots. Likely coincidental co-reporting or confounding.",
    ),
    ("Doxycycline", "Metabolic acidosis"): (
        "Plausible novel",
        "Tetracyclines (especially outdated formulations) can cause Fanconi "
        "syndrome with proximal renal tubular acidosis. Modern doxycycline "
        "is less nephrotoxic but the pathway exists. May also reflect "
        "pyroglutamic acidosis confounding (see separate PT).",
    ),
    ("Doxycycline", "Sputum discoloured"): (
        "Noise",
        "No established mechanism. Tetracyclines can stain tissues but "
        "sputum discolouration is not a known effect. May reflect infection "
        "being treated rather than drug effect.",
    ),
    ("Doxycycline", "Rales"): (
        "Noise",
        "Auscultation finding reflecting pulmonary pathology (infection, "
        "oedema). Likely the indication being treated, not an ADR.",
    ),
    ("Doxycycline", "Product ineffective for unapproved use"): (
        "Noise",
        "Administrative PT reflecting off-label prescribing outcomes. Not a "
        "pharmacological adverse effect.",
    ),
    ("Doxycycline", "Product use in unapproved indication"): (
        "Noise",
        "Administrative PT. Reflects prescribing practice, not pharmacology.",
    ),
    ("Doxycycline", "Treatment failure"): (
        "Noise",
        "Efficacy PT, not a safety signal. May co-occur with pathogen "
        "resistance signal.",
    ),
    ("Doxycycline", "Toxicity to various agents"): (
        "Noise",
        "Nonspecific PT. Cannot attribute to doxycycline specifically. "
        "Likely a catch-all coding artefact for multi-drug reports.",
    ),
    ("Doxycycline", "Ankle arthroplasty"): (
        "Noise",
        "Surgical procedure PT. No pharmacological mechanism. Likely coding "
        "artefact (procedure co-reported with perioperative drug) or "
        "prophylactic doxycycline use around surgery.",
    ),
    ("Doxycycline", "Gastric mucosal lesion"): (
        "Plausible novel",
        "GI mucosal injury. While oesophageal irritation is well-established, "
        "gastric mucosal damage is less characterised for doxycycline. "
        "Biologically plausible as extension of direct chemical irritation.",
    ),
    ("Doxycycline", "Product monitoring error"): (
        "Noise",
        "Administrative PT. Not a pharmacological effect.",
    ),
}


def main():
    print("=" * 70)
    print("  Clinical Face-Validity Review (Phase 7) — DAEN Tropical PV")
    print("=" * 70)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    signals = pd.read_csv(RESULTS_DIR / "signals_primary.csv")
    print(f"\nTotal signals to review: {len(signals)}")

    rows = []
    unrated = []
    for _, sig in signals.iterrows():
        key = (sig["canonical_name"], sig["reaction_pt"])
        if key in BLINDED_RATINGS:
            rating, rationale = BLINDED_RATINGS[key]
        else:
            rating = "UNRATED"
            rationale = "Not found in blinded rating dictionary"
            unrated.append(key)

        rows.append({
            "product_id": sig["product_id"],
            "product_set": sig["product_set"],
            "canonical_name": sig["canonical_name"],
            "reaction_pt": sig["reaction_pt"],
            "n": sig["n"],
            "prr": sig["prr"],
            "ic025": sig["ic025"],
            "eb05": sig["eb05"],
            "blinded_rating": rating,
            "clinical_rationale": rationale,
        })

    review_df = pd.DataFrame(rows)

    # Summary statistics
    rating_counts = review_df["blinded_rating"].value_counts()
    print("\nBlinded ratings:")
    for rating, count in rating_counts.items():
        print(f"  {rating:20s}  {count}")

    if unrated:
        print(f"\nWARNING: {len(unrated)} signals not in rating dictionary:")
        for drug, pt in unrated:
            print(f"  {drug} × {pt}")

    # By product set
    print("\nBy product set:")
    for ps, grp in review_df.groupby("product_set"):
        counts = grp["blinded_rating"].value_counts().to_dict()
        print(f"  Set {ps}: {dict(counts)}")

    # Key findings
    established = review_df[review_df["blinded_rating"] == "Established"]
    plausible = review_df[review_df["blinded_rating"] == "Plausible novel"]
    noise = review_df[review_df["blinded_rating"] == "Noise"]

    print(f"\n--- Established signals ({len(established)}) ---")
    for _, row in established.sort_values(["canonical_name", "prr"], ascending=[True, False]).iterrows():
        print(f"  {row['canonical_name']:35s} × {row['reaction_pt']:35s} PRR={row['prr']:.1f}")

    print(f"\n--- Plausible novel signals ({len(plausible)}) ---")
    for _, row in plausible.sort_values(["canonical_name", "prr"], ascending=[True, False]).iterrows():
        print(f"  {row['canonical_name']:35s} × {row['reaction_pt']:35s} PRR={row['prr']:.1f}")

    # Save
    review_df.to_csv(RESULTS_DIR / "clinical_review_log.csv", index=False)

    summary_lines = [
        "CLINICAL FACE-VALIDITY REVIEW SUMMARY (Phase 7, Section 9.7)",
        "=" * 60,
        "",
        f"Total signals reviewed: {len(review_df)}",
        "",
        "Blinded ratings:",
    ]
    for rating, count in rating_counts.items():
        pct = 100 * count / len(review_df)
        summary_lines.append(f"  {rating:20s}  {count:3d}  ({pct:.0f}%)")

    summary_lines.extend([
        "",
        f"Established: {len(established)} (validated known associations)",
        f"Plausible novel: {len(plausible)} (biologically plausible, requiring further investigation)",
        f"Noise: {len(noise)} (artefact, confounding, or non-pharmacological PTs)",
        "",
        "Novel signals by product (candidates for literature comparison):",
    ])
    for _, row in plausible.iterrows():
        summary_lines.append(f"  {row['canonical_name']} × {row['reaction_pt']} (PRR={row['prr']:.1f}, n={row['n']})")

    summary_text = "\n".join(summary_lines)
    (RESULTS_DIR / "clinical_review_summary.txt").write_text(summary_text)

    print(f"\n{summary_text}")
    print(f"\nSaved to: {RESULTS_DIR}/")
    print(f"\n{'=' * 70}")
    print("  Clinical face-validity review complete (Phase 7).")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
