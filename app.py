# -*- coding: utf-8 -*-
import streamlit as st
import pandas as pd
import io
import tempfile
import os
import shutil
import gzip
import re
import time
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio import Restriction
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import plotly.graph_objects as go
import primer3

st.set_page_config(page_title="AutoCAPS: Advanced Marker Designer", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
    <style>
    .metric-box {
        background-color: rgba(16, 185, 129, 0.1);
        border: 1px solid #10b981;
        border-radius: 8px;
        padding: 15px;
        text-align: center;
        margin-bottom: 10px;
    }
    .metric-value { font-size: 2rem; font-weight: bold; color: #10b981; }
    .metric-label { font-size: 1rem; font-weight: 600; }
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

st.markdown("<h1>🧬 AutoCAPS: Advanced CAPS Marker Designer</h1>", unsafe_allow_html=True)
st.markdown("### Standalone Bioinformatics Pipeline for Automated Marker Conversion")
st.markdown("<div class='credit-text'>System Architect & Lead Developer: Assoc. Prof. Dr. Ibrahim Celik</div>", unsafe_allow_html=True)
st.divider()

def parse_vcf(file_obj):
    lines = []
    for line in io.TextIOWrapper(file_obj, encoding='utf-8', errors='replace'):
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
    
    profiles, enz_list = [], []
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

def design_marker(left_flank, ref, alt, right_flank, snp_id, p_args, f_len):
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
        'PRIMER_OPT_SIZE': p_args['size_opt'],
        'PRIMER_MIN_SIZE': p_args['size_min'],
        'PRIMER_MAX_SIZE': p_args['size_max'],
        'PRIMER_OPT_TM': p_args['tm_opt'],
        'PRIMER_MIN_TM': p_args['tm_min'],
        'PRIMER_MAX_TM': p_args['tm_max'],
        'PRIMER_MIN_GC': p_args['gc_min'],
        'PRIMER_OPT_GC_PERCENT': p_args['gc_opt'],
        'PRIMER_MAX_GC': p_args['gc_max'],
        'PRIMER_PRODUCT_SIZE_RANGE': [p_args['prod_size']],
        
        'PRIMER_SALT_MONOVALENT': p_args['salt_mono'],
        'PRIMER_SALT_DIVALENT': p_args['salt_div'],
        'PRIMER_DNTP_CONC': p_args['salt_dntp'],
        'PRIMER_DNA_CONC': p_args['dna_conc'],
        'PRIMER_MAX_SELF_ANY_TH': p_args['max_self_any'],
        'PRIMER_MAX_SELF_END_TH': p_args['max_self_end'],
        'PRIMER_MAX_HAIRPIN_TH': p_args['max_hairpin'],
        'PRIMER_MAX_NS_ACCEPTED': p_args['max_ns'],
        'PRIMER_MAX_POLY_X': p_args['max_poly_x'],
        'PRIMER_CG_CLAMP': p_args['cg_clamp'],
        'PRIMER_MAX_GC_IN_3END': p_args['max_gc_3end'],
        
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_OPT_SIZE': p_args['io_s_opt'],
        'PRIMER_INTERNAL_MIN_SIZE': p_args['io_s_min'],
        'PRIMER_INTERNAL_MAX_SIZE': p_args['io_s_max'],
        'PRIMER_INTERNAL_OPT_TM': p_args['io_t_opt'],
        'PRIMER_INTERNAL_MIN_TM': p_args['io_t_min'],
        'PRIMER_INTERNAL_MAX_TM': p_args['io_t_max'],
        'PRIMER_INTERNAL_MIN_GC': p_args['io_g_min'],
        'PRIMER_INTERNAL_OPT_GC_PERCENT': p_args['io_g_opt'],
        'PRIMER_INTERNAL_MAX_GC': p_args['io_g_max']
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
    except Exception as e:
        pass
        
    return "Primer Not Found", "Primer Not Found", 0, 0, 0, "-", "-", formatli_snp, "No Suitable Amplicon", 0, []

def plot_vertical_map(df, chrom_lengths, selected_enzyme):
    fig = go.Figure()
    
    if selected_enzyme == "All Successful Enzymes":
        plot_df = df[df['Compatible_Enzyme_Count'] > 0]
    else:
        plot_df = df[df['Used_Enzyme_List'].apply(lambda x: selected_enzyme in x)]

    valid_chroms = sorted([ch for ch in chrom_lengths.keys() if chrom_lengths[ch] > 0])
    
    max_len_mb = 0
    for i, chrom in enumerate(valid_chroms):
        length_mb = chrom_lengths[chrom] / 1_000_000
        if length_mb > max_len_mb: max_len_mb = length_mb
        
        fig.add_trace(go.Scatter(
            x=[chrom, chrom], y=[0, length_mb], 
            mode='lines', line=dict(color='lightgray', width=8), 
            showlegend=False, hoverinfo='skip'
        ))
        
        chrom_df = plot_df[plot_df['FASTA_Header'] == chrom]
        if not chrom_df.empty:
            y_vals = chrom_df['POS'] / 1_000_000
            hover_texts = chrom_df.apply(lambda r: f"SNP: {r['SNP_ID']}<br>POS: {r['POS']}<br>Enzymes: {len(r['Used_Enzyme_List'])}", axis=1)
            
            fig.add_trace(go.Scatter(
                x=[chrom] * len(y_vals), y=y_vals,
                mode='markers',
                marker=dict(symbol='line-ew', color='#10b981', size=25, line=dict(width=3, color='#10b981')),
                text=hover_texts, hoverinfo='text', showlegend=False
            ))

    fig.update_layout(
        title=f"Interactive Physical Map: {selected_enzyme}",
        yaxis_title="Physical Position (Mb)",
        yaxis=dict(autorange='reversed', range=[max_len_mb * 1.05, - (max_len_mb * 0.05)]),
        xaxis_title="Chromosomes",
        plot_bgcolor='rgba(0,0,0,0)',
        height=700,
        margin=dict(l=40, r=40, t=60, b=40)
    )
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
    fig, ax = plt.subplots(figsize=(6, 8))
    
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
        well = patches.Rectangle((i - 0.25, well_y), 0.5, max_bp * 0.03, fill=False, edgecolor='black', lw=2)
        ax.add_patch(well)
        ax.text(i, well_y + (max_bp * 0.04), "WELL", color='gray', ha='center', fontsize=9, fontweight='bold')
        
    for i, (lane_name, bands) in enumerate(lanes.items()):
        x = i + 1
        for band in bands:
            alpha_core = 0.9 if lane_name != "DNA\nMarker" else 0.5
            b_color = '#10b981' if lane_name != "DNA\nMarker" else 'gray'
            ax.hlines(band, x - 0.25, x + 0.25, colors=b_color, linewidth=18, alpha=alpha_core * 0.3)
            ax.hlines(band, x - 0.2, x + 0.2, colors='black', linewidth=4, alpha=alpha_core) 
            if lane_name == "DNA\nMarker":
                ax.text(x - 0.35, band, f"{band}", color='black', va='center', ha='right', fontsize=10)

    ax.set_xticks(list(range(1, len(lanes) + 1)))
    ax.set_xticklabels(lanes.keys(), fontsize=11, fontweight='bold')
    ax.set_ylim(0, max_bp + (max_bp * 0.15))
    ax.set_ylabel("DNA Size (bp)", fontsize=12, fontweight='bold')
    ax.set_title(f"{enzyme} Restriction Profile", color='#10b981', fontweight='bold', fontsize=14, pad=20)
    return fig


# --- SIDEBAR SETTINGS ---
st.sidebar.header("📁 Data Input")
with st.sidebar.form("settings_form"):
    vcf_file = st.file_uploader("1. VCF File (.vcf)", type=['vcf'])
    fasta_file = st.file_uploader("2. Reference Genome (.fasta, .fa, .gz)", type=['fasta', 'fa', 'fna', 'gz'])
    
    st.markdown("### 🌡️ General PCR Settings")
    flank_len = st.number_input("FASTA Flank Window (bp):", min_value=150, max_value=1000, value=250, step=50)
    f_len = st.slider("Formatted SNP Length (bp):", 21, 201, 51, step=10)
    
    st.markdown("**Product Size Range (bp)**")
    col_p1, col_p2 = st.columns(2)
    prod_min = col_p1.number_input("Min", value=100, step=10, key='prmin')
    prod_max = col_p2.number_input("Max", value=300, step=10, key='prmax')

    st.markdown("**Primer Size (nt)**")
    col_s1, col_s2, col_s3 = st.columns(3)
    p_s_min = col_s1.number_input("Min", value=18, key='psmin')
    p_s_opt = col_s2.number_input("Opt", value=20, key='psopt')
    p_s_max = col_s3.number_input("Max", value=27, key='psmax')
    
    st.markdown("**Primer Tm (°C)**")
    col_t1, col_t2, col_t3 = st.columns(3)
    p_tm_min = col_t1.number_input("Min", value=57.0, step=0.5, key='ptmin')
    p_tm_opt = col_t2.number_input("Opt", value=60.0, step=0.5, key='ptopt')
    p_tm_max = col_t3.number_input("Max", value=63.0, step=0.5, key='ptmax')
    tm_diff = st.number_input("Max Tm Difference (°C):", value=2.0, step=0.1)
    
    st.markdown("**Primer GC (%)**")
    col_g1, col_g2, col_g3 = st.columns(3)
    p_gc_min = col_g1.number_input("Min", value=20.0, step=1.0, key='pgmin')
    p_gc_opt = col_g2.number_input("Opt", value=50.0, step=1.0, key='pgopt')
    p_gc_max = col_g3.number_input("Max", value=80.0, step=1.0, key='pgmax')

    with st.expander("⚙️ Advanced Primer3 Settings", expanded=False):
        adv_tabs = st.tabs(["Thermodynamics & Salt", "Sequence Penalties", "Internal Oligo (Probe)"])
        
        with adv_tabs[0]:
            st.markdown("**Salt Concentrations (mM)**")
            salt_mono = st.number_input("Monovalent Cations", value=50.0)
            salt_div = st.number_input("Divalent Cations", value=1.5)
            salt_dntp = st.number_input("dNTPs", value=0.6)
            dna_conc = st.number_input("Annealing Oligo (nM)", value=50.0)
            
            st.markdown("**Thermodynamic Alignments**")
            max_self_any = st.number_input("Max Self Complementarity", value=45.0)
            max_self_end = st.number_input("Max 3' Self Comp.", value=35.0)
            max_hairpin = st.number_input("Max Primer Hairpin", value=24.0)

        with adv_tabs[1]:
            max_ns = st.number_input("Max #N's accepted", value=0, step=1)
            max_poly_x = st.number_input("Max Poly-X", value=5, step=1)
            cg_clamp = st.number_input("CG Clamp", value=0, step=1)
            max_gc_3end = st.number_input("Max GC in 3' end", value=5, step=1)
            
        with adv_tabs[2]:
            st.markdown("**Internal Oligo Size (nt)**")
            c1, c2, c3 = st.columns(3)
            io_s_min = c1.number_input("Min", value=18, key='ismin')
            io_s_opt = c2.number_input("Opt", value=20, key='isopt')
            io_s_max = c3.number_input("Max", value=27, key='ismax')
            
            st.markdown("**Internal Oligo Tm (°C)**")
            c1, c2, c3 = st.columns(3)
            io_t_min = c1.number_input("Min", value=57.0, key='itmin')
            io_t_opt = c2.number_input("Opt", value=60.0, key='itopt')
            io_t_max = c3.number_input("Max", value=63.0, key='itmax')
            
            st.markdown("**Internal Oligo GC (%)**")
            c1, c2, c3 = st.columns(3)
            io_g_min = c1.number_input("Min", value=20.0, key='igmin')
            io_g_opt = c2.number_input("Opt", value=50.0, key='igopt')
            io_g_max = c3.number_input("Max", value=80.0, key='igmax')

    primer3_args = {
        'prod_size': [prod_min, prod_max],
        'size_min': p_s_min, 'size_opt': p_s_opt, 'size_max': p_s_max,
        'tm_min': p_tm_min, 'tm_opt': p_tm_opt, 'tm_max': p_tm_max,
        'gc_min': p_gc_min, 'gc_opt': p_gc_opt, 'gc_max': p_gc_max,
        'salt_mono': salt_mono, 'salt_div': salt_div, 'salt_dntp': salt_dntp, 'dna_conc': dna_conc,
        'max_self_any': max_self_any, 'max_self_end': max_self_end, 'max_hairpin': max_hairpin,
        'max_ns': max_ns, 'max_poly_x': max_poly_x, 'cg_clamp': cg_clamp, 'max_gc_3end': max_gc_3end,
        'io_s_min': io_s_min, 'io_s_opt': io_s_opt, 'io_s_max': io_s_max,
        'io_t_min': io_t_min, 'io_t_opt': io_t_opt, 'io_t_max': io_t_max,
        'io_g_min': io_g_min, 'io_g_opt': io_g_opt, 'io_g_max': io_g_max
    }

    run_btn = st.form_submit_button("🚀 Execute Pipeline")

if run_btn:
    if not vcf_file or not fasta_file:
        st.error("Please upload both VCF and FASTA files.")
    else:
        # Uyarı mesajı ve güvenli sayaç alanı
        st.warning("⚠️ **Pipeline is running!** Depending on your VCF size, this may take a few minutes. **Please DO NOT close or refresh this tab.**")
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        start_time = time.time()
        status_text.markdown("⚙️ **0% Completed** | ⏱️ **Elapsed Time:** 0 sec (Initializing...)")
        
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
                chrom_lengths = {}
                total_snps = len(df)
                batch_size = 100  # 100'erli güvenli paketler halinde işleme
                
                elapsed_secs = int(time.time() - start_time)
                status_text.markdown(f"⚙️ **0% Completed** (0/{total_snps} SNPs) | ⏱️ **Elapsed Time:** {elapsed_secs} sec")
                
                for start_idx in range(0, total_snps, batch_size):
                    end_idx = min(start_idx + batch_size, total_snps)
                    batch_df = df.iloc[start_idx:end_idx]
                    
                    for _, row in batch_df.iterrows():
                        chrom, pos, ref, alt, snp_id = str(row['CHROM']), int(row['POS']), str(row['REF']), str(row['ALT']), str(row['ID'])
                        actual_chrom = find_fasta_chrom(chrom, fasta_keys)
                        
                        if actual_chrom:
                            if actual_chrom not in chrom_lengths:
                                chrom_lengths[actual_chrom] = len(fasta_obj[actual_chrom])
                                
                            pos_0 = pos - 1
                            actual_ref_base = str(fasta_obj[actual_chrom][pos_0 : pos_0 + len(ref)].seq).upper()
                            
                            if actual_ref_base != ref.upper():
                                results.append({
                                    'SNP_ID': snp_id, 'CHROM': chrom, 'FASTA_Header': actual_chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
                                    'Forward_Primer': "REF Mismatch", 'Reverse_Primer': "-", 'PCR_Size_bp': 0,
                                    'Probe_REF': "-", 'Probe_ALT': "-", 'Formatted_SNP_[REF/ALT]': "-",
                                    'Compatible_Enzyme_Count': 0, 'Used_Enzyme_List': [], 
                                    'CAPS_Restriction_Profile': f"ERROR: VCF REF ({ref}) != FASTA ({actual_ref_base})"
                                })
                                continue

                            left_flank = str(fasta_obj[actual_chrom][max(0, pos_0 - flank_len):pos_0].seq).upper()
                            right_flank = str(fasta_obj[actual_chrom][pos_0 + len(ref):pos_0 + len(ref) + flank_len].seq).upper()
                            
                            f_seq, r_seq, f_tm, r_tm, p_size_val, p_ref, p_alt, form_snp, c_prof, enz_c, enz_l = design_marker(
                                left_flank, ref, alt, right_flank, snp_id, primer3_args, f_len
                            )
                            
                            all_enzymes.extend(enz_l)
                            results.append({
                                'SNP_ID': snp_id, 'CHROM': chrom, 'FASTA_Header': actual_chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
                                'Forward_Primer': f_seq, 'Reverse_Primer': r_seq, 'PCR_Size_bp': p_size_val,
                                'Probe_REF': p_ref, 'Probe_ALT': p_alt, 'Formatted_SNP_[REF/ALT]': form_snp,
                                'Compatible_Enzyme_Count': enz_c, 'Used_Enzyme_List': enz_l, 'CAPS_Restriction_Profile': c_prof
                            })
                        else:
                            results.append({
                                'SNP_ID': snp_id, 'CHROM': chrom, 'FASTA_Header': "-", 'POS': pos, 'REF': ref, 'ALT': alt,
                                'Forward_Primer': "Not Found in FASTA", 'Reverse_Primer': "-", 'PCR_Size_bp': 0,
                                'Probe_REF': "-", 'Probe_ALT': "-", 'Formatted_SNP_[REF/ALT]': "-",
                                'Compatible_Enzyme_Count': 0, 'Used_Enzyme_List': [], 
                                'CAPS_Restriction_Profile': f"VCF '{chrom}' did not match FASTA"
                            })
                    
                    # HER 100LÜK BATCH SONUNDA EKRANI VE SÜREYİ GÜNCELLE
                    processed_count = end_idx
                    elapsed_secs = int(time.time() - start_time)
                    mins, secs = divmod(elapsed_secs, 60)
                    time_str = f"{mins} min {secs} sec" if mins > 0 else f"{secs} sec"
                    
                    pct = int((processed_count / total_snps) * 100)
                    progress_bar.progress(pct / 100.0)
                    status_text.markdown(f"⚙️ **{pct}% Completed** ({processed_count}/{total_snps} SNPs) | ⏱️ **Elapsed Time:** {time_str}")
                
                elapsed_secs = int(time.time() - start_time)
                mins, secs = divmod(elapsed_secs, 60)
                time_str = f"{mins} min {secs} sec" if mins > 0 else f"{secs} sec"
                
                status_text.markdown(f"✅ **100% Completed** ({total_snps}/{total_snps} SNPs) in **{time_str}** - Rendering results...")
                time.sleep(0.5)
                
                st.session_state['result_df'] = pd.DataFrame(results)
                st.session_state['all_enzymes'] = all_enzymes
                st.session_state['chrom_lengths'] = chrom_lengths
                
            except Exception as e:
                st.error(f"An error occurred during processing: {str(e)}")
            finally:
                fasta_obj.close()
                os.remove(fasta_path)
                if os.path.exists(fasta_path + ".fai"): os.remove(fasta_path + ".fai")

# --- DASHBOARD UI (OTURUM HAFIZASINDAN BESLENİR, KAYBOLMAZ) ---
if 'result_df' in st.session_state:
    res_df = st.session_state['result_df']
    all_enz = st.session_state['all_enzymes']
    chrom_lengths = st.session_state['chrom_lengths']
    
    st.markdown("## 🔍 Targeted Analysis & Filtering")
    c1, c2, c3, c4 = st.columns(4)
    all_chroms = ["All"] + sorted(list(res_df['FASTA_Header'].unique()))
    if "-" in all_chroms: all_chroms.remove("-")
    
    sel_chrom = c1.selectbox("Target Chromosome", all_chroms)
    sel_start = c2.number_input("Region Start (bp)", value=0, min_value=0)
    sel_end = c3.number_input("Region End (bp) (0 = All)", value=0, min_value=0)
    
    success_df = res_df[res_df['Compatible_Enzyme_Count'] > 0]
    unique_enzymes = list(set(all_enz))
    all_enz_opts = ["All Successful Enzymes"] + sorted(unique_enzymes)
    sel_enz = c4.selectbox("Restriction Enzyme", all_enz_opts)
    
    filt_df = res_df.copy()
    if sel_chrom != "All":
        filt_df = filt_df[filt_df['FASTA_Header'] == sel_chrom]
        if sel_end > 0:
            filt_df = filt_df[(filt_df['POS'] >= sel_start) & (filt_df['POS'] <= sel_end)]
    
    if sel_enz != "All Successful Enzymes":
        filt_df = filt_df[filt_df['Used_Enzyme_List'].apply(lambda x: sel_enz in x)]
        
    st.divider()

    st.markdown("## 📊 1. Genome & Chromosome Statistics")
    
    total_snp = len(filt_df)
    mapped_snp = len(filt_df[filt_df['FASTA_Header'] != "-"])
    primer_success = len(filt_df[(filt_df['Forward_Primer'] != "Primer Not Found") & (filt_df['Forward_Primer'] != "Not Found in FASTA") & (filt_df['Forward_Primer'] != "REF Mismatch")])
    caps_success = len(filt_df[filt_df['Compatible_Enzyme_Count'] > 0])
    
    m1, m2, m3, m4 = st.columns(4)
    m1.markdown(f"<div class='metric-box'><div class='metric-value'>{total_snp}</div><div class='metric-label'>Filtered VCF SNPs</div></div>", unsafe_allow_html=True)
    m2.markdown(f"<div class='metric-box'><div class='metric-value'>{mapped_snp}</div><div class='metric-label'>Genome Matched</div></div>", unsafe_allow_html=True)
    m3.markdown(f"<div class='metric-box'><div class='metric-value'>{primer_success}</div><div class='metric-label'>Successful Primers</div></div>", unsafe_allow_html=True)
    m4.markdown(f"<div class='metric-box'><div class='metric-value'>{caps_success}</div><div class='metric-label'>CAPS Convertible</div></div>", unsafe_allow_html=True)
    
    stats_list = []
    total_genome_len = sum(chrom_lengths.values())
    
    for ch, length in chrom_lengths.items():
        ch_df = res_df[(res_df['FASTA_Header'] == ch) & (res_df['Compatible_Enzyme_Count'] > 0)]
        marker_count = len(ch_df)
        resolution = round((length / 1000) / marker_count, 2) if marker_count > 0 else "N/A"
        stats_list.append({"Chromosome": ch, "Length (Mb)": round(length/1_000_000, 2), "CAPS Markers": marker_count, "Resolution (kb/Marker)": resolution})
    
    total_markers = sum([s["CAPS Markers"] for s in stats_list])
    global_res = round((total_genome_len / 1000) / total_markers, 2) if total_markers > 0 else "N/A"
    stats_list.append({"Chromosome": "TOTAL GENOME", "Length (Mb)": round(total_genome_len/1_000_000, 2), "CAPS Markers": total_markers, "Resolution (kb/Marker)": global_res})
    
    st.table(pd.DataFrame(stats_list))
    
    st.divider()
    
    st.markdown("## 🌍 2. Physical Chromosome Map")
    if caps_success > 0:
        fig_map = plot_vertical_map(filt_df, chrom_lengths, sel_enz)
        st.plotly_chart(fig_map, use_container_width=True)
    else:
        st.info("No successful CAPS markers found to map with current filters.")
        
    st.divider()
    
    st.markdown("## 🖥️ 3. Gel Electrophoresis Simulation")
    gel_candidates = filt_df[filt_df['Compatible_Enzyme_Count'] > 0]
    
    if gel_candidates.empty:
        st.warning("No suitable polymorphic cut profile found for gel simulation under current filters.")
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
            
            selected_snp_enz = st.selectbox("Select Target Enzyme:", enz_list)
            frag_dict = parse_profile_for_gel(profile_str, selected_snp_enz, ref_allele, alt_allele, pcr_size)
            
            st.info(f"**Target SNP:** {selected_snp_id}\n\n"
                    f"**PCR Amplicon Size:** {pcr_size} bp\n\n"
                    f"**Allele 1 ({ref_allele}):** {frag_dict.get(ref_allele, [pcr_size])} bp\n\n"
                    f"**Allele 2 ({alt_allele}):** {frag_dict.get(alt_allele, [pcr_size])} bp")
                    
        with col_gel2:
            fig_gel = render_gel_simulator(frag_dict, ref_allele, alt_allele, pcr_size, selected_snp_enz)
            st.pyplot(fig_gel)

    st.divider()
    st.markdown("## 💾 4. Export Filtered Data")
    st.dataframe(filt_df.drop(columns=['Used_Enzyme_List']))
    
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        export_df = filt_df.copy()
        export_df['Used_Enzyme_List'] = export_df['Used_Enzyme_List'].apply(lambda x: ", ".join(x))
        export_df.to_excel(writer, sheet_name='Filtered_Markers', index=False)
    st.download_button("📂 Download Filtered Report as Excel", output.getvalue(), f"AutoCAPS_{sel_chrom}_{sel_enz}.xlsx")
