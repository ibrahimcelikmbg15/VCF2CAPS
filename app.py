# -*- coding: utf-8 -*-
import streamlit as st
import pandas as pd
import io
import tempfile
import os
import shutil
import gzip
import re
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio import Restriction
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import primer3

st.set_page_config(page_title="GBS2CAPS: TaqMan & CAPS Marker Designer", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
    <style>
    .metric-box {
        background-color: #f8fafc;
        border: 1px solid #cbd5e1;
        border-radius: 8px;
        padding: 15px;
        text-align: center;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05);
    }
    @media (prefers-color-scheme: dark) {
        .metric-box {
            background-color: #1e293b;
            border: 1px solid #334155;
        }
    }
    .metric-value { font-size: 2rem; font-weight: bold; color: #10b981; }
    .metric-label { font-size: 1rem; color: #64748b; }
    .credit-text {
        font-size: 1.1rem;
        color: #0ea5e9;
        font-weight: 700;
        border-left: 4px solid #0ea5e9;
        padding-left: 15px;
        margin-top: 5px;
    }
    </style>
""", unsafe_allow_html=True)

col_logo, col_title = st.columns([1, 4])
with col_logo:
    logo_url = "https://raw.githubusercontent.com/ibrahimcelikmbg15/-SNP-Spesifik-CAPS-Mark-r-Tasarlay-c-/main/genom.jpeg"
    try:
        st.image(logo_url, use_container_width=True)
    except:
        st.markdown("🧬")

with col_title:
    st.markdown("<h1>🧬 GBS2CAPS: TaqMan & CAPS Marker Designer</h1>", unsafe_allow_html=True)
    st.markdown("### **Genom Okulu** | Bioinformatics and Molecular Breeding Platform")
    st.markdown("<div class='credit-text'>System Architect & Lead Developer: Assoc. Prof. Dr. Ibrahim Celik</div>", unsafe_allow_html=True)

def parse_vcf(file_obj):
    lines = []
    for line in io.TextIOWrapper(file_obj, encoding='utf-8'):
        if line.startswith('##'): continue
        lines.append(line)
    df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')
    df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)
    if 'ID' in df.columns:
        df['ID'] = df.apply(lambda row: f"{row['CHROM']}_{row['POS']}" if row['ID'] == '.' or pd.isna(row['ID']) else row['ID'], axis=1)
    return df

def find_fasta_chrom(chrom, fasta_keys):
    if chrom in fasta_keys: return chrom
    chrom_base = chrom.split('.')[0] if '.' in chrom else chrom
    for k in fasta_keys:
        k_base = k.split(' ')[0]
        k_base_no_ver = k_base.split('.')[0] if '.' in k_base else k_base
        if chrom_base == k_base_no_ver: return k
        if chrom_base.lower().replace('chr', '') == k_base_no_ver.lower().replace('chr', ''): return k
    return None

def get_fragments(cut_positions, total_len):
    if not cut_positions: return f"{total_len}bp (No Cut)"
    frags, last_pos = [], 0
    for cut in sorted(cut_positions):
        frags.append(cut - last_pos)
        last_pos = cut
    frags.append(total_len - last_pos)
    return " + ".join([f"{f}bp" for f in frags if f > 0])

def get_frag_list(cuts, total):
    if not cuts: return [total]
    f, last = [], 0
    for c in sorted(cuts):
        f.append(c - last)
        last = c
    f.append(total - last)
    return [x for x in f if x > 0]

def analyze_amplicon(amplicon_ref, amplicon_alt, ref, alt):
    batch = Restriction.CommOnly
    seq_ref, seq_alt = Seq(str(amplicon_ref)), Seq(str(amplicon_alt))
    res_ref, res_alt = batch.search(seq_ref), batch.search(seq_alt)
    
    profiles = []
    enz_list = []
    for enz in batch:
        c_ref, c_alt = res_ref.get(enz, []), res_alt.get(enz, [])
        if c_ref != c_alt:
            if len(c_ref) > 2 or len(c_alt) > 2: continue
            f_ref_list, f_alt_list = get_frag_list(c_ref, len(amplicon_ref)), get_frag_list(c_alt, len(amplicon_alt))
            if any(x < 15 for x in f_ref_list) or any(x < 15 for x in f_alt_list): continue
            
            frag_ref = get_fragments(c_ref, len(amplicon_ref))
            frag_alt = get_fragments(c_alt, len(amplicon_alt))
            profiles.append(f"{enz.__name__} -> Allele {ref}: [{frag_ref}] | Allele {alt}: [{frag_alt}]")
            enz_list.append(enz.__name__)
            
    return " || ".join(profiles) if profiles else "No Suitable/Selective Enzyme Found", len(profiles), enz_list

def design_marker(left_flank, ref, alt, right_flank, snp_id, p_size, p_tm, tm_diff, p_gc, p_prod, f_len):
    half_f = f_len // 2
    c_left = left_flank[-half_f:] if len(left_flank) >= half_f else left_flank
    formatli_snp = f"{c_left}[{ref}/{alt}]{right_flank[:half_f]}"
    
    seq_ref = left_flank + ref + right_flank
    seq_alt = left_flank + alt + right_flank
    target_pos = len(left_flank)
    
    seq_dict = {
        'SEQUENCE_ID': str(snp_id),
        'SEQUENCE_TEMPLATE': seq_ref,
        'SEQUENCE_TARGET': [target_pos, len(ref)]
    }
    primer_opt = {
        'PRIMER_OPT_SIZE': p_size,
        'PRIMER_MIN_SIZE': max(15, p_size - 5),
        'PRIMER_MAX_SIZE': min(35, p_size + 5),
        'PRIMER_OPT_TM': p_tm,
        'PRIMER_MIN_TM': p_tm - tm_diff,
        'PRIMER_MAX_TM': p_tm + tm_diff,
        'PRIMER_PRODUCT_SIZE_RANGE': [list(p_prod)],
        'PRIMER_MIN_GC': p_gc[0],
        'PRIMER_MAX_GC': p_gc[1]
    }
    
    try:
        primers = primer3.bindings.design_primers(seq_dict, primer_opt)
        if primers['PRIMER_LEFT_NUM_RETURNED'] > 0:
            f_seq = primers['PRIMER_LEFT_0_SEQUENCE']
            r_seq = primers['PRIMER_RIGHT_0_SEQUENCE']
            f_tm_val = round(primers['PRIMER_LEFT_0_TM'], 1)
            r_tm_val = round(primers['PRIMER_RIGHT_0_TM'], 1)
            pcr_size = primers['PRIMER_PAIR_0_PRODUCT_SIZE']
            
            left_pos = primers['PRIMER_LEFT_0'][0]
            right_pos = primers['PRIMER_RIGHT_0'][0]
            
            amplicon_ref = seq_ref[left_pos:right_pos+1]
            amplicon_alt = seq_alt[left_pos:right_pos+1]
            
            amplicon_ref = re.sub(r'[^ACGTURYSWKMBDHVN]', 'N', amplicon_ref.upper())
            amplicon_alt = re.sub(r'[^ACGTURYSWKMBDHVN]', 'N', amplicon_alt.upper())
            
            snp_amp_pos = target_pos - left_pos
            probe_start, probe_end = max(0, snp_amp_pos - 10), min(len(amplicon_ref), snp_amp_pos + 11)
            p_ref = amplicon_ref[probe_start:probe_end]
            p_alt = amplicon_alt[probe_start:probe_end]
            
            c_prof, enz_c, enz_l = analyze_amplicon(amplicon_ref, amplicon_alt, ref, alt)
            return f_seq, r_seq, f_tm_val, r_tm_val, pcr_size, p_ref, p_alt, formatli_snp, c_prof, enz_c, enz_l
    except:
        pass
        
    return "Primer Not Found", "Primer Not Found", 0, 0, 0, "-", "-", formatli_snp, "No Suitable Amplicon", 0, []

def plot_physical_map(df, selected_enzyme):
    fig, ax = plt.subplots(figsize=(12, min(8, len(df['CHROM'].unique()) * 0.8)))
    
    if selected_enzyme == "All Successful Enzymes":
        plot_df = df[df['Compatible_Enzyme_Count'] > 0]
    else:
        plot_df = df[df['Used_Enzyme_List'].apply(lambda x: selected_enzyme in x)]

    if plot_df.empty:
        ax.text(0.5, 0.5, "No SNPs found matching this filter.", ha='center', va='center')
        ax.axis('off')
        return fig

    chroms = sorted(plot_df['CHROM'].unique())
    y_ticks, y_labels = [], []
    
    for i, chrom in enumerate(chroms):
        y_pos = i
        y_ticks.append(y_pos)
        y_labels.append(chrom)
        chrom_df = plot_df[plot_df['CHROM'] == chrom]
        
        max_pos = chrom_df['POS'].max() / 1_000_000
        ax.hlines(y_pos, 0, max_pos + 1, color='lightgray', linewidth=5, zorder=1) 
        ax.plot(chrom_df['POS'] / 1_000_000, [y_pos]*len(chrom_df), '|', color='#10b981', markersize=20, markeredgewidth=2, zorder=2)

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("Physical Position (Mb)")
    ax.set_title(f"Genome Distribution: {selected_enzyme}", fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.tight_layout()
    return fig

def parse_profile_for_gel(profile_str, enzyme, ref, alt, pcr_size):
    frag_data = {ref: [pcr_size], alt: [pcr_size]}
    if not profile_str or "No Suitable" in profile_str:
        return frag_data
    
    enzymes_profiles = profile_str.split(" || ")
    target_prof = None
    for prof in enzymes_profiles:
        if prof.startswith(enzyme + " ->"):
            target_prof = prof
            break
            
    if not target_prof: return frag_data
    
    alleles_split = target_prof.split(" | ")
    for part in alleles_split:
        m = re.search(r'Allele (.*?): \[(.*?)\]', part)
        if m:
            al = m.group(1).strip()
            sizes_str = m.group(2)
            if "No Cut" in sizes_str:
                frag_data[al] = [pcr_size]
            else:
                sizes = [int(s.replace("bp", "").strip()) for s in sizes_str.split("+")]
                frag_data[al] = sizes
    return frag_data

def render_gel_simulator(frag_data, ref, alt, pcr_size, enzyme):
    bg_color = '#020617' 
    fig, ax = plt.subplots(figsize=(6, 8), facecolor=bg_color)
    ax.set_facecolor(bg_color)
    
    max_bp = ((pcr_size // 100) + 1) * 100 if pcr_size % 100 != 0 else pcr_size
    marker_bands = sorted(list(range(50, max_bp + 100, 50)), reverse=True)
    
    frags1 = frag_data.get(ref, [pcr_size])
    frags2 = frag_data.get(alt, [pcr_size])
    
    lanes = {
        "DNA\nMarker": marker_bands,
        f"Allele 1\n({ref}{ref})": frags1,
        f"Allele 2\n({alt}{alt})": frags2,
        f"Hybrid\n({ref}{alt})": sorted(list(set(frags1 + frags2)), reverse=True)
    }
    
    well_y = max_bp + 50
    for i in range(1, len(lanes) + 1):
        well = patches.Rectangle((i - 0.25, well_y), 0.5, max_bp * 0.03, fill=False, edgecolor='#334155', lw=2)
        ax.add_patch(well)
        ax.text(i, well_y + (max_bp * 0.04), "WELL", color='#94a3b8', ha='center', fontsize=9, fontweight='bold')
        
    for i, (lane_name, bands) in enumerate(lanes.items()):
        x = i + 1
        for band in bands:
            alpha_core = 0.9 if lane_name != "DNA\nMarker" else 0.5
            neon_color = '#10b981' if lane_name != "DNA\nMarker" else '#38bdf8'
            ax.hlines(band, x - 0.25, x + 0.25, colors=neon_color, linewidth=18, alpha=alpha_core * 0.25)
            ax.hlines(band, x - 0.2, x + 0.2, colors='#ffffff', linewidth=4, alpha=alpha_core) 
            if lane_name == "DNA\nMarker":
                ax.text(x - 0.35, band, f"{band}", color='#cbd5e1', va='center', ha='right', fontsize=10)

    ax.set_xticks(list(range(1, len(lanes) + 1)))
    ax.set_xticklabels(lanes.keys(), color='#ffffff', fontsize=11, fontweight='bold')
    ax.set_ylim(0, max_bp + (max_bp * 0.15))
    ax.tick_params(colors='#ffffff', axis='y', labelsize=10)
    for spine in ax.spines.values():
        spine.set_color('#1e293b')
    ax.set_ylabel("DNA Size (bp)", color='#ffffff', fontsize=12, fontweight='bold')
    ax.set_title(f"{enzyme} Restriction Profile", color='#10b981', fontweight='bold', fontsize=14, pad=20)
    return fig

# --- SIDEBAR SETTINGS ---
st.sidebar.header("🌡️ Marker & Probe Settings")
with st.sidebar.form("settings_form"):
    vcf_file = st.file_uploader("1. VCF File (.vcf)", type=['vcf'])
    fasta_file = st.file_uploader("2. Reference Genome (.fasta, .fa, .gz)", type=['fasta', 'fa', 'fna', 'gz'])
    
    flank_len = st.number_input("FASTA Flank Window (bp):", min_value=150, max_value=1000, value=250, step=50)
    f_len = st.slider("Formatted SNP Sequence Length (bp):", 21, 201, 51, step=10)
    p_prod = st.slider("Target PCR Size (bp):", 50, 600, (100, 300), step=10)
    p_size = st.slider("Ideal Primer Length (nt):", 15, 30, 20)
    p_tm = st.slider("Ideal Melting Temp (Tm °C):", 50.0, 72.0, 60.0, step=0.5)
    tm_diff = st.number_input("Max Tm Difference (°C):", min_value=0.1, max_value=10.0, value=2.0, step=0.1)
    p_gc = st.slider("GC Content Range (%):", 20.0, 80.0, (40.0, 60.0))
    run_btn = st.form_submit_button("🚀 Run Pipeline")

if run_btn:
    if not vcf_file or not fasta_file:
        st.error("Please upload both VCF and FASTA files.")
    else:
        df = parse_vcf(vcf_file)
        valid_bases = ['A', 'C', 'G', 'T']
        df = df[~df['ALT'].str.contains(',', na=False)]
        df = df[df['REF'].isin(valid_bases) & df['ALT'].isin(valid_bases)]
        df = df.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']).reset_index(drop=True)
        
        if len(df) == 0: 
            st.warning("No suitable Biallelic SNP found in this VCF file!")
        else:
            file_extension = ".fa.gz" if fasta_file.name.endswith(".gz") else ".fasta"
            with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp_fasta:
                if fasta_file.name.endswith(".gz"):
                    with gzip.open(fasta_file, 'rb') as f_in:
                        shutil.copyfileobj(f_in, tmp_fasta)
                else:
                    shutil.copyfileobj(fasta_file, tmp_fasta)
                fasta_path = tmp_fasta.name
            
            try:
                fasta_obj = Fasta(fasta_path)
                fasta_keys = list(fasta_obj.keys())
                
                results, all_enzymes = [], []
                progress_bar = st.progress(0)
                
                for i, row in df.iterrows():
                    chrom, pos, ref, alt, snp_id = str(row['CHROM']), int(row['POS']), str(row['REF']), str(row['ALT']), str(row['ID'])
                    actual_chrom = find_fasta_chrom(chrom, fasta_keys)
                    
                    if actual_chrom:
                        pos_0 = pos - 1
                        actual_ref_base = str(fasta_obj[actual_chrom][pos_0 : pos_0 + len(ref)].seq).upper()
                        
                        if actual_ref_base != ref.upper():
                            results.append({
                                'SNP_ID': snp_id, 'CHROM': chrom, 'FASTA_Header': actual_chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
                                'Forward_Primer': "REF Mismatch", 'Reverse_Primer': "-", 'PCR_Size_bp': 0,
                                'TaqMan_Probe_REF': "-", 'TaqMan_Probe_ALT': "-", 'Formatted_SNP_[REF/ALT]': "-",
                                'Compatible_Enzyme_Count': 0, 'Used_Enzyme_List': [], 'CAPS_Restriction_Profile': f"ERROR: VCF REF ({ref}) != FASTA ({actual_ref_base})"
                            })
                            continue

                        left_flank = str(fasta_obj[actual_chrom][max(0, pos_0 - flank_len):pos_0].seq).upper()
                        right_flank = str(fasta_obj[actual_chrom][pos_0 + len(ref):pos_0 + len(ref) + flank_len].seq).upper()
                        
                        f_seq, r_seq, f_tm, r_tm, p_size_val, p_ref, p_alt, form_snp, c_prof, enz_c, enz_l = design_marker(
                            left_flank, ref, alt, right_flank, snp_id, p_size, p_tm, tm_diff, p_gc, p_prod, f_len
                        )
                        all_enzymes.extend(enz_l)
                        results.append({
                            'SNP_ID': snp_id, 'CHROM': chrom, 'FASTA_Header': actual_chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
                            'Forward_Primer': f_seq, 'Reverse_Primer': r_seq, 'PCR_Size_bp': p_size_val,
                            'TaqMan_Probe_REF': p_ref, 'TaqMan_Probe_ALT': p_alt, 'Formatted_SNP_[REF/ALT]': form_snp,
                            'Compatible_Enzyme_Count': enz_c, 'Used_Enzyme_List': enz_l, 'CAPS_Restriction_Profile': c_prof
                        })
                    else:
                        results.append({
                            'SNP_ID': snp_id, 'CHROM': chrom, 'FASTA_Header': "-", 'POS': pos, 'REF': ref, 'ALT': alt,
                            'Forward_Primer': "Not Found in FASTA", 'Reverse_Primer': "-", 'PCR_Size_bp': 0,
                            'TaqMan_Probe_REF': "-", 'TaqMan_Probe_ALT': "-", 'Formatted_SNP_[REF/ALT]': "-",
                            'Compatible_Enzyme_Count': 0, 'Used_Enzyme_List': [], 'CAPS_Restriction_Profile': f"VCF '{chrom}' did not match FASTA"
                        })
                    
                    if i % 10 == 0 or i == len(df) - 1: progress_bar.progress(min((i + 1) / len(df), 1.0))
                
                st.session_state['result_df'] = pd.DataFrame(results)
                st.session_state['all_enzymes'] = all_enzymes
                
            except Exception as e:
                st.error(f"An error occurred during processing: {str(e)}")
            finally:
                fasta_obj.close()
                os.remove(fasta_path)
                if os.path.exists(fasta_path + ".fai"): os.remove(fasta_path + ".fai")

# --- DASHBOARD UI ---
if 'result_df' in st.session_state:
    res_df = st.session_state['result_df']
    all_enz = st.session_state['all_enzymes']
    
    st.markdown("## 📊 1. Genome Statistics")
    total_snp = len(res_df)
    mapped_snp = len(res_df[res_df['FASTA_Header'] != "-"])
    primer_success = len(res_df[(res_df['Forward_Primer'] != "Primer Not Found") & (res_df['Forward_Primer'] != "Not Found in FASTA") & (res_df['Forward_Primer'] != "REF Mismatch")])
    caps_success = len(res_df[res_df['Compatible_Enzyme_Count'] > 0])
    
    m1, m2, m3, m4 = st.columns(4)
    m1.markdown(f"<div class='metric-box'><div class='metric-value'>{total_snp}</div><div class='metric-label'>Processed VCF SNPs</div></div>", unsafe_allow_html=True)
    m2.markdown(f"<div class='metric-box'><div class='metric-value'>{mapped_snp}</div><div class='metric-label'>Matched in Genome</div></div>", unsafe_allow_html=True)
    m3.markdown(f"<div class='metric-box'><div class='metric-value'>{primer_success}</div><div class='metric-label'>Successful Primers</div></div>", unsafe_allow_html=True)
    m4.markdown(f"<div class='metric-box'><div class='metric-value'>{caps_success}</div><div class='metric-label'>CAPS Conversions</div></div>", unsafe_allow_html=True)
    
    unique_enzymes = list(set(all_enz))
    st.info(f"🧬 **Enzyme Screening Summary:** A total of **{len(unique_enzymes)} different restriction enzymes** produced polymorphic cuts across the genome.")
    
    st.divider()
    
    st.markdown("## 🌍 2. Chromosome Physical Map")
    col_map_sel, col_map_plot = st.columns([1, 4])
    with col_map_sel:
        st.markdown("Filter enzyme distribution:")
        map_enz_list = ["All Successful Enzymes"] + sorted(unique_enzymes)
        selected_map_enz = st.selectbox("Map Filter:", map_enz_list)
    with col_map_plot:
        fig_map = plot_physical_map(res_df, selected_map_enz)
        st.pyplot(fig_map)
        
    st.divider()
    
    st.markdown("## 🖥️ 3. Gel Electrophoresis Simulation")
    gel_candidates = res_df[res_df['Compatible_Enzyme_Count'] > 0]
    
    if gel_candidates.empty:
        st.warning("No suitable polymorphic cut profile found for gel simulation.")
    else:
        col_gel1, col_gel2 = st.columns([1, 2])
        with col_gel1:
            selected_snp_id = st.selectbox("Select SNP to Visualize:", gel_candidates['SNP_ID'].tolist())
            snp_data = gel_candidates[gel_candidates['SNP_ID'] == selected_snp_id].iloc[0]
            
            ref_allele = snp_data['REF']
            alt_allele = snp_data['ALT']
            pcr_size = snp_data['PCR_Size_bp']
            enz_list = snp_data['Used_Enzyme_List']
            profile_str = snp_data['CAPS_Restriction_Profile']
            
            selected_snp_enz = st.selectbox("Select Enzyme to Use:", enz_list)
            
            frag_dict = parse_profile_for_gel(profile_str, selected_snp_enz, ref_allele, alt_allele, pcr_size)
            
            st.info(f"**SNP ID:** {selected_snp_id}\n\n"
                    f"**PCR Product:** {pcr_size} bp\n\n"
                    f"**Allele 1 ({ref_allele}):** {frag_dict.get(ref_allele, [pcr_size])} bp\n\n"
                    f"**Allele 2 ({alt_allele}):** {frag_dict.get(alt_allele, [pcr_size])} bp")
                    
        with col_gel2:
            fig_gel = render_gel_simulator(frag_dict, ref_allele, alt_allele, pcr_size, selected_snp_enz)
            st.pyplot(fig_gel)

    st.divider()
    st.markdown("## 💾 4. Data Export")
    st.dataframe(res_df.drop(columns=['Used_Enzyme_List']))
    
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        export_df = res_df.copy()
        export_df['Used_Enzyme_List'] = export_df['Used_Enzyme_List'].apply(lambda x: ", ".join(x))
        export_df.to_excel(writer, sheet_name='Marker_Designs', index=False)
    st.download_button("📂 Download Detailed Report as Excel", output.getvalue(), "GBS2CAPS_Analysis_Report.xlsx")
