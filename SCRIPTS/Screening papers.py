# Google Scholar search for papers on elevated CO2 effects on grain quality,
# ionome, biomass and photosynthesis in C3 crops.
# excludes simulations, models and meta-analyses.

import pandas as pd
from scholarly import scholarly
import time
import random
import re
from datetime import datetime


query = (
    # target C3 crops (grain crops; no model plants)
    "(wheat OR rice OR barley OR rye OR oat OR "
    "soybean OR chickpea OR lentil OR bean OR pea) "

    # treatment: elevated CO2
    "AND (\"elevated CO2\" OR eCO2 OR \"high CO2\" OR \"CO2 enrichment\" OR "
    "\"atmospheric CO2\" OR \"elevated atmospheric CO2\") "

    # measured responses
    "AND ("
    "\"grain quality\" OR \"grain composition\" OR ionomics OR ionomic OR ionome OR "
    "\"mineral content\" OR \"protein content\" OR \"starch content\" OR "
    "\"amino acid\" OR \"nitrogen content\" OR \"N assimilation\" OR "
    "\"carbon assimilation\" OR \"C:N ratio\" OR "
    "photosynthesis OR \"photosynthetic rate\" OR \"Rubisco\" OR "
    "\"Rubisco abundance\" OR Vcmax OR Jmax OR \"stomatal conductance\" OR "
    "\"water use efficiency\" OR WUE OR \"transpiration rate\" OR "
    "biomass OR \"shoot biomass\" OR \"root biomass\" OR \"dry weight\" OR "
    "\"nutrient uptake\" OR \"micronutrient\" OR iron OR zinc OR magnesium OR "
    "phosphorus OR potassium OR calcium"
    ") "

    # exclude: simulations, models, reviews, model plants, other stresses
    "-simulation -model -modeling -modelling -review -\"meta-analysis\" "
    "-\"meta analysis\" -\"meta?analysis\" -metaanalysis "
    "-C4 -maize -corn -\"C4 plant\" -ozone -salt -salinity -drought -flooding "
    "-\"low water\" -infection -pathogen -disease "
    "-arabidopsis -Arabidopsis -Brassica -tomato -poplar -tobacco"
)



EXCLUDE_KEYWORDS = [
    'simulation', 'model', 'modeling', 'modelling', 'meta-analysis', 'meta analysis',
    'metaanalysis', 'review', 'theoretical', 'experiment', 'virtual',
    'machine learning', 'artificial neural', 'dynamic systems', 'crop model',
    'arabidopsis', 'brassica', 'tomato', 'poplar', 'tobacco'
]

QUALITY_KEYWORDS = [
    'ionomic', 'ionomics', 'ionome', 'grain quality', 'grain composition',
    'mineral content', 'nutritional quality', 'protein', 'nitrogen assimilation',
    'photosynthesis', 'stomatal', 'water use efficiency'
]



def is_experimental_study(result):
    # Keep experimental studies; drop simulations/models
    title = result.get('bib', {}).get('title', '').lower()
    abstract = result.get('bib', {}).get('abstract', '').lower()

    full_text = title + ' ' + abstract
    for word in EXCLUDE_KEYWORDS:
        if word in full_text:
            return False

    for word in QUALITY_KEYWORDS:
        if word in title:
            return True

    return True

def has_useful_data(result):
    # Require title, year and a minimal abstract
    bib = result.get('bib', {})
    has_title = bool(bib.get('title', '').strip())
    has_year = bool(bib.get('pub_year', ''))
    has_abstract = len(bib.get('abstract', '').strip()) > 50
    return has_title and has_year and has_abstract

def extract_doi(result):
    # Extract DOI from several possible sources
    bib = result.get('bib', {})

    doi = bib.get('doi', '').strip()
    if doi:
        return doi

    url = result.get('pub_url', '')
    if url:
        match = re.search(r'10\.\d{4,9}/[-._;()/:a-zA-Z0-9]+', str(url))
        if match:
            return match.group(0)

    for field in [bib.get('abstract', ''), str(bib)]:
        match = re.search(r'10\.\d{4,9}/[-._;()/:a-zA-Z0-9]+', str(field))
        if match:
            return match.group(0)

    return "N/A"

def extract_authors(result):
    authors = result.get('bib', {}).get('author', '')
    if authors:
        if isinstance(authors, list):
            return '; '.join([str(a) for a in authors[:3]])
        elif isinstance(authors, str):
            return authors[:200]
    return "N/A"

def process_result(result, index):
    try:
        bib = result.get('bib', {})

        title = bib.get('title', 'N/A')
        year = bib.get('pub_year', 'N/A')
        authors = extract_authors(result)
        journal = bib.get('venue', 'N/A')
        doi = extract_doi(result)
        url = result.get('pub_url', '')
        citations = result.get('num_citations', 0)
        abstract = bib.get('abstract', '')

        if not has_useful_data(result):
            return None

        if not is_experimental_study(result):
            return None

        variables = detect_variables(abstract)

        data = {
            'Index': index,
            'Title': title,
            'Year': year,
            'Authors (first 3)': authors,
            'Journal/Venue': journal,
            'DOI': doi,
            'URL': url,
            'Citations': citations,
            'Measured variables': variables,
            'Study type': detect_study_type(title, abstract),
            'Plant species': detect_species(title, abstract),
            'CO2 treatment': detect_co2(title, abstract),
            'Abstract': abstract[:300] + '...' if len(abstract) > 300 else abstract
        }

        return data

    except Exception as e:
        print(f"  Error processing result {index}: {str(e)}")
        return None

def detect_variables(text):
    if not text:
        return "Not specified"

    found = []

    variable_map = {
        'ionomic': 'Ionome',
        'mineral content': 'Mineral content',
        'photosynthesis': 'Photosynthesis',
        'stomatal': 'Stomatal conductance',
        'WUE': 'Water use efficiency',
        'water use efficiency': 'Water use efficiency',
        'grain yield': 'Grain yield',
        'biomass': 'Biomass',
        'protein': 'Protein',
        'nitrogen': 'Nitrogen',
        'rubisco': 'Rubisco',
        'vcmax': 'Vcmax',
        'jmax': 'Jmax'
    }

    text_lower = text.lower()
    for var, label in variable_map.items():
        if var in text_lower and label not in found:
            found.append(label)

    return '; '.join(found) if found else "Not specified"

def detect_study_type(title, abstract):
    text = (title + ' ' + abstract).lower()

    if 'field' in text:
        return 'Field'
    elif 'greenhouse' in text or 'growth chamber' in text:
        return 'Greenhouse/Chamber'
    elif 'pot' in text:
        return 'Pot'
    elif 'hydroponic' in text:
        return 'Hydroponic'
    else:
        return 'Not specified'

def detect_species(title, abstract):
    species = {
        'wheat': 'Wheat',
        'rice': 'Rice',
        'barley': 'Barley',
        'rye': 'Rye',
        'oat': 'Oat',
        'soybean': 'Soybean',
        'chickpea': 'Chickpea',
        'bean': 'Bean',
        'lentil': 'Lentil',
        'pea': 'Pea'
    }

    text = (title + ' ' + abstract).lower()
    found = []

    for eng, label in species.items():
        if eng in text and label not in found:
            found.append(label)

    return '; '.join(found) if found else 'Not specified'

def detect_co2(title, abstract):
    text = (title + ' ' + abstract).lower()

    if 'eCO2' in title + abstract or 'e-CO2' in title + abstract:
        return 'Elevated CO2 (eCO2)'
    elif 'elevated' in text:
        return 'Elevated CO2'
    elif 'high CO2' in text:
        return 'High CO2'
    else:
        return 'Not specified'



def search_papers(query, limit=500):
    print("=" * 80)
    print("PAPER SEARCH FOR META-ANALYSIS")
    print("Elevated CO2 effects on grain quality, ionome, biomass and photosynthesis (C3 crops)")
    print("=" * 80)
    print(f"\nQuery: {query}\n")
    print(f"Starting search (limit: {limit} papers)...\n")

    search_query = scholarly.search_pubs(query)
    results = []
    attempts = 0
    consecutive_errors = 0

    while attempts < limit:
        try:
            res = next(search_query)
            attempts += 1

            data = process_result(res, attempts)

            if data:
                results.append(data)
                print(f"[{attempts}] {data['Title'][:60]}... ({data['Year']})")
                consecutive_errors = 0
            else:
                print(f"[{attempts}] Filtered out (does not meet criteria)")

            if attempts % 10 == 0:
                print(f"\n  -> Processed {attempts} results, {len(results)} valid")
                wait = random.uniform(8, 15)
                print(f"  -> Pausing {wait:.1f}s to avoid blocking...\n")
                time.sleep(wait)
            else:
                time.sleep(random.uniform(1, 3))

        except StopIteration:
            print(f"\nReached all available results after {attempts} attempts")
            break

        except Exception as e:
            consecutive_errors += 1

            if consecutive_errors >= 3:
                print(f"\nToo many consecutive errors. Stopping search.")
                print(f"  Error: {str(e)[:100]}")
                break

            print(f"Error on attempt {attempts}: {str(e)[:80]}")
            wait = 30 + (consecutive_errors * 20)
            print(f"  Waiting {wait}s before retrying...")
            time.sleep(wait)

    return results


def save_results(papers_found):
    if not papers_found:
        print("\nNo papers matched the criteria.")
        return

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    df = pd.DataFrame(papers_found)

    excel_file = f'papers_CO2_C3_metaanalysis_{timestamp}.xlsx'
    try:
        df.to_excel(excel_file, index=False, engine='openpyxl')
        print(f"\nExcel saved: {excel_file}")
    except Exception as e:
        print(f"\nError saving Excel: {e}")

    csv_file = f'papers_CO2_C3_metaanalysis_{timestamp}.csv'
    df.to_csv(csv_file, index=False, encoding='utf-8-sig')
    print(f"CSV saved: {csv_file}")

    print("\n" + "=" * 80)
    print("SEARCH STATISTICS")
    print("=" * 80)
    print(f"Total papers found: {len(df)}")
    print(f"\nBy year:")
    print(df['Year'].value_counts().sort_index(ascending=False).head(10))
    print(f"\nBy journal:")
    print(df['Journal/Venue'].value_counts().head(10))
    print(f"\nStudy type:")
    print(df['Study type'].value_counts())
    print(f"\nSpecies studied:")
    print(df['Plant species'].value_counts())

    return excel_file, csv_file


if __name__ == "__main__":
    try:
        papers_found = search_papers(query, limit=1500)

        if papers_found:
            excel_file, csv_file = save_results(papers_found)
            print(f"\nSEARCH COMPLETE")


    except KeyboardInterrupt:
        print("\n\nSearch interrupted by user")
    except Exception as e:
        print(f"\n\nFatal error: {e}")
