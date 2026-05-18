import streamlit as st
import pandas as pd
import requests
import time
import re
from bs4 import BeautifulSoup
from io import BytesIO

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="HRMS Annotator - Zingiber", page_icon="🌿", layout="wide")

# --- UTILITY & API FUNCTIONS ---
@st.cache_data(show_spinner=False, max_entries=500)
def query_pubchem(formula):
    """Query PubChem PUG REST API by formula."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/{formula}/JSON"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            compounds = data.get('PC_Compounds', [])
            if not compounds:
                return []
            
            # Extract basic data for top 3 isomers to avoid massive payloads
            results = []
            for cid_data in compounds[:3]:
                cid = cid_data.get('id', {}).get('id', {}).get('cid', 'Unknown')
                props = cid_data.get('props', [])
                name = next((p['value']['sval'] for p in props if p['urn']['label'] == 'IUPAC Name'), f"CID: {cid}")
                results.append({"cid": cid, "name": name})
            return results
        elif response.status_code == 503: # Rate limited
            time.sleep(2)
            return query_pubchem(formula) # Retry once
    except Exception as e:
        pass
    return []

@st.cache_data(show_spinner=False, max_entries=500)
def query_knapsack(formula):
    """Scrape KNApSAcK DB for formula and species occurrence."""
    url = f"http://www.knapsackfamily.com/knapsack_core/info.php?sname=formula&word={formula}"
    results = []
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            soup = BeautifulSoup(response.content, 'html.parser')
            # KNApSAcK results are typically in a table
            tables = soup.find_all('table')
            for table in tables:
                rows = table.find_all('tr')
                for row in rows[1:]: # Skip header
                    cols = row.find_all('td')
                    if len(cols) >= 6:
                        c_id = cols[0].text.strip()
                        name = cols[2].text.strip()
                        organism = cols[5].text.strip()
                        results.append({"c_id": c_id, "name": name, "organism": organism})
    except Exception:
        pass
    return results

@st.cache_data(show_spinner=False, max_entries=500)
def query_chemspider(formula, api_key):
    """Query ChemSpider API (requires RSC API Key)."""
    if not api_key:
        return []
    # Using the standard ChemSpider REST API endpoint structure
    url = "https://api.rsc.org/compounds/v1/filter/formula"
    headers = {"apikey": api_key, "Content-Type": "application/json", "Accept": "application/json"}
    try:
        # 1. Initiate search
        res = requests.post(url, headers=headers, json={"formula": formula}, timeout=10)
        if res.status_code == 200:
            query_id = res.json().get('queryId')
            # 2. Wait and get results
            time.sleep(1) 
            status_url = f"https://api.rsc.org/compounds/v1/filter/{query_id}/results"
            res_status = requests.get(status_url, headers=headers)
            if res_status.status_code == 200:
                results = res_status.json().get('results', [])
                return [f"CSID: {csid}" for csid in results[:3]] # Return top 3 IDs
    except Exception:
        pass
    return []

def calculate_zingiber_score(pubchem_res, knapsack_res):
    """
    Weighted scoring mechanism:
    +50 if 'Zingiber montanum' found in KNApSAcK organism list.
    +30 if 'Zingiber' (genus) found in KNApSAcK organism list or PubChem name.
    +10 if compound is generally found in KNApSAcK (confirms natural product origin).
    """
    score = 0
    tags = []
    
    # Check KNApSAcK
    is_natural_product = len(knapsack_res) > 0
    if is_natural_product:
        score += 10
        tags.append("Natural Product")
        
    for k_item in knapsack_res:
        organism = k_item['organism'].lower()
        if 'zingiber montanum' in organism:
            score += 50
            tags.append("Z. montanum (Exact)")
            break # Max score achieved for this tier
        elif 'zingiber' in organism:
            score += 30
            tags.append("Zingiber (Genus)")
            break
            
    # Check PubChem names for common ginger compound markers (e.g., zingiberene, zerumbone)
    for p_item in pubchem_res:
        name = p_item['name'].lower()
        if 'zingiber' in name or 'zerumb' in name or 'ginger' in name:
            if "Zingiber (Genus)" not in tags:
                score += 30
                tags.append("Zingiber (Name Match)")
                
    return score, list(set(tags))

def clean_formula(formula):
    """Remove spaces from formulas (e.g., 'C11 H12 O3' -> 'C11H12O3')."""
    if pd.isna(formula):
        return ""
    return re.sub(r'\s+', '', str(formula))

# --- UI ---
st.title("🌿 HRMS Annotator: *Zingiber montanum* Pipeline")
st.markdown("Upload your High-Resolution Mass Spectrometry data. The app will query PubChem, KNApSAcK, and ChemSpider, scoring candidate compounds based on their likelihood of originating from *Zingiber montanum*.")

with st.sidebar:
    st.header("⚙️ Configuration")
    chemspider_key = st.text_input("ChemSpider API Key (Optional)", type="password", help="Requires an RSC Developer account. If left blank, ChemSpider queries will be skipped.")
    st.markdown("---")
    st.markdown("**Required Columns:**\n* `Formula`\n* `Calc. MW` / `Calc. Molecular Weight`\n* `RT [min]`\n* `Sample Area`\n* `DeltaMass [ppm]`")

uploaded_file = st.file_uploader("Upload HRMS Dataset (.csv or .xlsx)", type=['csv', 'xlsx'])

if uploaded_file is not None:
    try:
        # Handle file types and arbitrary header row locations
        if uploaded_file.name.endswith('.csv'):
            df_raw = pd.read_csv(uploaded_file, header=None)
        else:
            df_raw = pd.read_excel(uploaded_file, header=None)
            
        # Find the actual header row by looking for 'Formula'
        header_idx = df_raw[df_raw.apply(lambda r: r.astype(str).str.contains('Formula', case=False, na=False).any(), axis=1)].index
        if len(header_idx) > 0:
            df = df_raw.iloc[header_idx[0]+1:].copy()
            df.columns = df_raw.iloc[header_idx[0]]
        else:
            df = df_raw.copy()
            
        df = df.reset_index(drop=True)
        
        # Standardize column names dynamically based on user prompt & provided file structure
        col_mapping = {}
        for col in df.columns:
            col_str = str(col).lower()
            if 'formula' in col_str: col_mapping[col] = 'Formula'
            elif 'delta' in col_str and 'mass' in col_str: col_mapping[col] = 'DeltaMass [ppm]'
            elif 'calc' in col_str and 'mw' in col_str or 'molecular weight' in col_str: col_mapping[col] = 'Calc. Molecular Weight'
            elif 'rt' in col_str and 'min' in col_str: col_mapping[col] = 'RT [min]'
            elif 'area' in col_str: col_mapping[col] = 'Sample Area'
            elif 'annot' in col_str and 'mass' not in col_str: col_mapping[col] = 'Annotation'
            
        df = df.rename(columns=col_mapping)
        
        # Ensure 'Formula' exists
        if 'Formula' not in df.columns:
            st.error("Could not find a 'Formula' column. Please check your dataset.")
            st.stop()
            
        st.success(f"File parsed successfully. Found {len(df)} rows.")
        st.dataframe(df.head(5))
        
        if st.button("🚀 Start Annotation Pipeline"):
            annotated_data = []
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            for idx, row in df.iterrows():
                # Progress UI
                progress = (idx + 1) / len(df)
                progress_bar.progress(progress)
                status_text.text(f"Processing row {idx+1}/{len(df)}: {row.get('Formula', 'Unknown')}")
                
                raw_formula = row.get('Formula', '')
                formula = clean_formula(raw_formula)
                
                if not formula:
                    continue
                    
                # 1. API Queries
                pubchem_hits = query_pubchem(formula)
                knapsack_hits = query_knapsack(formula)
                chemspider_hits = query_chemspider(formula, chemspider_key)
                
                # 2. Scoring & Filtering
                score, tags = calculate_zingiber_score(pubchem_hits, knapsack_hits)
                
                # 3. Compile Results
                pc_names = ", ".join([h['name'] for h in pubchem_hits]) if pubchem_hits else "None Found"
                ks_orgs = ", ".join(list(set([h['organism'] for h in knapsack_hits]))) if knapsack_hits else "None Found"
                cs_ids = ", ".join(chemspider_hits) if chemspider_hits else ("Skipped" if not chemspider_key else "None Found")
                
                row_dict = row.to_dict()
                row_dict['Cleaned Formula'] = formula
                row_dict['PubChem Candidates'] = pc_names
                row_dict['KNApSAcK Organisms'] = ks_orgs
                row_dict['ChemSpider IDs'] = cs_ids
                row_dict['Zingiber Score'] = score
                row_dict['Botanical Tags'] = ", ".join(tags)
                
                annotated_data.append(row_dict)
                
                # Respect rate limits between rows
                time.sleep(0.3)
                
            progress_bar.empty()
            status_text.success("Annotation Complete!")
            
            # Display output
            results_df = pd.DataFrame(annotated_data)
            # Sort by Zingiber Score descending
            results_df = results_df.sort_values(by='Zingiber Score', ascending=False).reset_index(drop=True)
            
            st.subheader("📊 Annotation Results")
            st.dataframe(results_df)
            
            # Download Button
            output = BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                results_df.to_excel(writer, index=False, sheet_name='Annotated_Results')
            processed_data = output.getvalue()
            
            st.download_button(
                label="📥 Download Annotated Dataset (.xlsx)",
                data=processed_data,
                file_name="Zingiber_Annotated_HRMS.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

    except Exception as e:
        st.error(f"An error occurred while processing the file: {e}")
