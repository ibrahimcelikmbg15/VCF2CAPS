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

st.set_page_config(page_title="GBS2CAPS: TaqMan & CAPS Tasarım", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
    <style>
    .stApp {background-color: #030712;}
    h1, h2, h3, p, label, li {color: #ffffff !important;}
    .stTextArea textarea {
        background-color: #020617 !important;
        color: #00ffcc !important; 
        font-family: 'Courier New', Consolas, monospace !important;
        font-size: 1.05rem !important;
        border: 2px solid #1e293b !important;
        border-radius: 8px;
    }
    .stTextArea textarea:focus {
        border: 2px solid #10b981 !important;
        box-shadow: 0 0 10px rgba(16, 185, 129, 0.5);
    }
    .stButton>button {
        background-color: transparent;
        color: #10b981 !important;
        font-weight: 800;
        font-size: 1.1rem;
        border-radius: 8px;
        border: 2px solid #10b981;
        box-shadow: 0 0 10px rgba(16, 185, 129, 0.3);
        transition: all 0.3s ease-in-out;
    }
    .stButton>button:hover {
        background-color: #10b981;
        color: #000000 !important;
        box-shadow: 0 0 20px rgba(16, 185, 129, 0.8);
        border: 2px solid #10b981;
    }
    .credit-text {
        font-size: 1.2rem;
        color: #00e5ff !important;
        font-weight: 800;
        border-left: 4px solid #00e5ff;
        padding-left: 15px;
        margin-top: 5px;
        letter-spacing: 0.5px;
    }
    .metric-box {
        background-color: #0f172a;
        border: 1px solid #1e293b;
        border-radius: 8px;
        padding: 15px;
        text-align: center;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.5);
    }
    .metric-value { font-size: 2rem; font-weight: bold; color: #10b981; }
    .metric-label { font-size: 1rem; color: #94a3b8; }
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
    st.markdown("<h1>🧬 SNP-Spesifik CAPS Markör Tasarlayıcı v2.0</h1>", unsafe_allow_html=True)
    st.markdown("### **Genom Okulu** | Biyoinformatik ve Moleküler Islah Platformu")
    st.markdown("<div class='credit-text'>Sistem Mimarı ve Baş Geliştirici: Doç. Dr. İbrahim Çelik</div>", unsafe_allow_html=True)

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

def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def calculate_tm(seq):
    return (4 * (seq.count('G') + seq.count('C'))) + (2 * (seq.count('A') + seq.count('T')))

def design_marker(left_flank, ref, alt, right_flank, snp_id, p_size, p_tm, tm_diff, p_gc, p_prod, f_len):
    p_size, min_prod, max_prod = int(p_size), int(p_prod[0]), int(p_prod[1])
    half_f = f_len // 2
    c_left = left_flank[-half_f:] if len(left_flank) >= half_f else left_flank
    formatli_snp = f"{c_left}[{ref}/{alt}]{right_flank[:half_f]}"
    
    valid_fwds, valid_revs = [], []
    for i in range(len(left_flank) - p_size + 1):
        cand = left_flank[i:i+p_size]
        if p_gc[0] <= ((cand.count('G') + cand.count('C')) / p_size) * 100 <= p_gc[1]:
            tm = calculate_tm(cand)
            if abs(tm - p_tm) <= tm_diff + 5: valid_fwds.append((i, cand, tm))
                
    for i in range(len(right_flank) - p_size + 1):
        cand, rev_cand = right_flank[i:i+p_size], get_reverse_complement(right_flank[i:i+p_size])
        if p_gc[0] <= ((rev_cand.count('G') + rev_cand.count('C')) / p_size) * 100 <= p_gc[1]:
            tm = calculate_tm(rev_cand)
            if abs(tm - p_tm) <= tm_diff + 5: valid_revs.append((i + p_size, rev_cand, tm))

    fwd_seq, rev_seq, fwd_tm, rev_tm, best_left, best_right, best_diff = "-", "-", 0, 0, 0, 0, 9999
    target_ideal = (min_prod + max_prod) / 2
    
    for l_pos, f_seq, f_tm in valid_fwds:
        for r_pos, r_seq, r_tm in valid_revs:
            if abs(f_tm - r_tm) <= tm_diff + 2:
                prod_size = (len(left_flank) - l_pos) + len(ref) + r_pos
                if min_prod <= prod_size <= max_prod:
                    diff = abs(prod_size - target_ideal)
                    if diff < best_diff:
                        best_diff, fwd_seq, fwd_tm, best_left, rev_seq, rev_tm, best_right = diff, f_seq, f_tm, l_pos, r_seq, r_tm, r_pos

    if fwd_seq != "-":
        prod_size = (len(left_flank) - best_left) + len(ref) + best_right
        amplicon_ref = left_flank[best_left:] + ref + right_flank[:best_right]
        amplicon_alt = left_flank[best_left:] + alt + right_flank[:best_right]
        
        amplicon_ref = re.sub(r'[^ACGTURYSWKMBDHVN]', 'N', amplicon_ref.upper())
        amplicon_alt = re.sub(r'[^ACGTURYSWKMBDHVN]', 'N', amplicon_alt.upper())
        
        snp_amp_pos = len(left_flank) - best_left
        probe_start, probe_end = max(0, snp_amp_pos - 10), min(len(amplicon_ref), snp_amp_pos + 11)
        
        cut_profile, enz_count, enz_list = analyze_amplicon(amplicon_ref, amplicon_alt, ref, alt)
        return (fwd_seq, rev_seq, round(fwd_tm, 1), round(rev_tm, 1), prod_size, amplicon_ref[probe_start:probe_end], amplicon_alt[probe_start:probe_end], formatli_snp, cut_profile, enz_count, enz_list)
        
    return ("Primer Bulunamadı", "Primer Bulunamadı", 0, 0, 0, "-", "-", formatli_snp, "Uygun Amplikon Yok", 0, [])

def plot_physical_map(df, selected_enzyme):
    fig, ax = plt.subplots(figsize=(12, min(8, len(df['CHROM'].unique()) * 0.8)))
    fig.patch.set_facecolor('#030712')
    ax.set_facecolor('#030712')
    
    if selected_enzyme == "Tüm Başarılı Enzimler":
        plot_df = df[df['Uyumlu_Enzim_Sayısı'] > 0]
    else:
        plot_df = df[df['Kullanılan_Enzim_Listesi'].apply(lambda x: selected_enzyme in x)]

    if plot_df.empty:
        ax.text(0.5, 0.5, "Bu filtreye uygun SNP bulunamadı.", color='white', ha='center', va='center')
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
        ax.hlines(y_pos, 0, max_pos + 1, color='#1e293b', linewidth=5, zorder=1) 
        ax.plot(chrom_df['POS'] / 1_000_000, [y_pos]*len(chrom_df), '|', color='#10b981', markersize=20, markeredgewidth=2, zorder=2)

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, color='white')
    ax.set_xlabel("Fiziksel Pozisyon (Mb)", color='white')
    ax.set_title(f"Genom Dağılımı: {selected_enzyme}", fontweight='bold', color='white')
    ax.tick_params(colors='white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_color('#1e293b')
    plt.tight_layout()
    return fig

def parse_profile_for_gel(profile_str, enzyme, ref, alt, pcr_size):
    # Örnek profile_str: "AanI -> Allele G: [200bp (No Cut)] | Allele A: [116bp + 84bp] || Baul -> ..."
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
        # Örnek parça: "AanI -> Allele G: [200bp (No Cut)]" veya "Allele A: [116bp + 84bp]"
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
        f"Alel 1\n({ref}{ref})": frags1,
        f"Alel 2\n({alt}{alt})": frags2,
        f"Melez\n({ref}{alt})": sorted(list(set(frags1 + frags2)), reverse=True)
    }
    
    well_y = max_bp + 50
    for i in range(1, len(lanes) + 1):
        well = patches.Rectangle((i - 0.25, well_y), 0.5, max_bp * 0.03, fill=False, edgecolor='#334155', lw=2)
        ax.add_patch(well)
        ax.text(i, well_y + (max_bp * 0.04), "KUYU", color='#94a3b8', ha='center', fontsize=9, fontweight='bold')
        
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
    ax.set_ylabel("DNA Büyüklüğü (bp)", color='#ffffff', fontsize=12, fontweight='bold')
    ax.set_title(f"{enzyme} Restriksiyon Profili", color='#10b981', fontweight='bold', fontsize=14, pad=20)
    return fig

# --- YAN PANEL AYARLARI ---
st.sidebar.header("🌡️ Markör & Prob Ayarları")
with st.sidebar.form("settings_form"):
    vcf_file = st.file_uploader("1. VCF Dosyası (.vcf)", type=['vcf'])
    fasta_file = st.file_uploader("2. Referans Genom (.fasta, .fa, .gz)", type=['fasta', 'fa', 'fna', 'gz'])
    
    flank_len = st.number_input("FASTA Okuma Penceresi (bp):", min_value=150, max_value=1000, value=250, step=50)
    f_len = st.slider("Formatlı SNP Dizisi Uzunluğu (bp):", 21, 201, 51, step=10)
    p_prod = st.slider("Hedef PCR Boyutu (bp):", 50, 600, (100, 300), step=10)
    p_size = st.slider("İdeal Primer Uzunluğu (nt):", 15, 30, 20)
    p_tm = st.slider("İdeal Erime Sıcaklığı (Tm °C):", 50.0, 72.0, 60.0, step=0.5)
    tm_diff = st.number_input("Maksimum Tm Farkı (°C):", min_value=0.1, max_value=10.0, value=2.0, step=0.1)
    p_gc = st.slider("GC Oranı Aralığı (%):", 20.0, 80.0, (40.0, 60.0))
    run_btn = st.form_submit_button("🚀 Boru Hattını Çalıştır")

if run_btn:
    if not vcf_file or not fasta_file:
        st.error("Lütfen VCF ve FASTA dosyalarını yükleyin.")
    else:
        df = parse_vcf(vcf_file)
        valid_bases = ['A', 'C', 'G', 'T']
        df = df[~df['ALT'].str.contains(',', na=False)]
        df = df[df['REF'].isin(valid_bases) & df['ALT'].isin(valid_bases)]
        df = df.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']).reset_index(drop=True)
        
        if len(df) == 0: 
            st.warning("Bu VCF dosyasında uygun Biallelik SNP bulunamadı!")
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
                                'SNP_ID': snp_id, 'CHROM': chrom, 'FASTA_Başlığı': actual_chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
                                'Forward_Primer': "REF Uyuşmazlığı", 'Reverse_Primer': "-", 'PCR_Boyutu_bp': 0,
                                'TaqMan_Prob_REF': "-", 'TaqMan_Prob_ALT': "-", 'Formatlı_SNP_[REF/ALT]': "-",
                                'Uyumlu_Enzim_Sayısı': 0, 'Kullanılan_Enzim_Listesi': [], 'CAPS_Kesim_Profili': f"HATA: VCF REF ({ref}) != FASTA ({actual_ref_base})"
                            })
                            continue

                        left_flank = str(fasta_obj[actual_chrom][max(0, pos_0 - flank_len):pos_0].seq).upper()
                        right_flank = str(fasta_obj[actual_chrom][pos_0 + len(ref):pos_0 + len(ref) + flank_len].seq).upper()
                        
                        f_seq, r_seq, f_tm, r_tm, p_size_val, p_ref, p_alt, form_snp, c_prof, enz_c, enz_l = design_marker(
                            left_flank, ref, alt, right_flank, snp_id, p_size, p_tm, tm_diff, p_gc, p_prod, f_len
                        )
                        all_enzymes.extend(enz_l)
                        results.append({
                            'SNP_ID': snp_id, 'CHROM': chrom, 'FASTA_Başlığı': actual_chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
                            'Forward_Primer': f_seq, 'Reverse_Primer': r_seq, 'PCR_Boyutu_bp': p_size_val,
                            'TaqMan_Prob_REF': p_ref, 'TaqMan_Prob_ALT': p_alt, 'Formatlı_SNP_[REF/ALT]': form_snp,
                            'Uyumlu_Enzim_Sayısı': enz_c, 'Kullanılan_Enzim_Listesi': enz_l, 'CAPS_Kesim_Profili': c_prof
                        })
                    else:
                        results.append({
                            'SNP_ID': snp_id, 'CHROM': chrom, 'FASTA_Başlığı': "-", 'POS': pos, 'REF': ref, 'ALT': alt,
                            'Forward_Primer': "FASTA'da Bulunamadı", 'Reverse_Primer': "-", 'PCR_Boyutu_bp': 0,
                            'TaqMan_Prob_REF': "-", 'TaqMan_Prob_ALT': "-", 'Formatlı_SNP_[REF/ALT]': "-",
                            'Uyumlu_Enzim_Sayısı': 0, 'Kullanılan_Enzim_Listesi': [], 'CAPS_Kesim_Profili': f"VCF '{chrom}' FASTA ile eşleşmedi"
                        })
                    
                    if i % 10 == 0 or i == len(df) - 1: progress_bar.progress(min((i + 1) / len(df), 1.0))
                
                st.session_state['result_df'] = pd.DataFrame(results)
                st.session_state['all_enzymes'] = all_enzymes
                
            except Exception as e:
                st.error(f"İşlem sırasında hata oluştu: {str(e)}")
            finally:
                fasta_obj.close()
                os.remove(fasta_path)
                if os.path.exists(fasta_path + ".fai"): os.remove(fasta_path + ".fai")

# --- DASHBOARD ARAYÜZÜ ---
if 'result_df' in st.session_state:
    res_df = st.session_state['result_df']
    all_enz = st.session_state['all_enzymes']
    
    st.markdown("## 📊 1. Genom İstatistikleri")
    total_snp = len(res_df)
    mapped_snp = len(res_df[res_df['FASTA_Başlığı'] != "-"])
    primer_success = len(res_df[(res_df['Forward_Primer'] != "Primer Bulunamadı") & (res_df['Forward_Primer'] != "FASTA'da Bulunamadı") & (res_df['Forward_Primer'] != "REF Uyuşmazlığı")])
    caps_success = len(res_df[res_df['Uyumlu_Enzim_Sayısı'] > 0])
    
    m1, m2, m3, m4 = st.columns(4)
    m1.markdown(f"<div class='metric-box'><div class='metric-value'>{total_snp}</div><div class='metric-label'>İşlenen VCF SNP</div></div>", unsafe_allow_html=True)
    m2.markdown(f"<div class='metric-box'><div class='metric-value'>{mapped_snp}</div><div class='metric-label'>Genomda Eşleşen</div></div>", unsafe_allow_html=True)
    m3.markdown(f"<div class='metric-box'><div class='metric-value'>{primer_success}</div><div class='metric-label'>Başarılı Primer</div></div>", unsafe_allow_html=True)
    m4.markdown(f"<div class='metric-box'><div class='metric-value'>{caps_success}</div><div class='metric-label'>CAPS Dönüşümü</div></div>", unsafe_allow_html=True)
    
    unique_enzymes = list(set(all_enz))
    st.info(f"🧬 **Enzim Tarama Özeti:** Toplam **{len(unique_enzymes)} farklı restriksiyon enzimi** genom genelinde polimorfik kesim noktası yakaladı.")
    
    st.divider()
    
    st.markdown("## 🌍 2. Kromozom Üzerinde Fiziksel Dağılım")
    col_map_sel, col_map_plot = st.columns([1, 4])
    with col_map_sel:
        st.markdown("Enzim dağılımını filtreleyin:")
        map_enz_list = ["Tüm Başarılı Enzimler"] + sorted(unique_enzymes)
        selected_map_enz = st.selectbox("Harita Filtresi:", map_enz_list)
    with col_map_plot:
        fig_map = plot_physical_map(res_df, selected_map_enz)
        st.pyplot(fig_map)
        
    st.divider()
    
    st.markdown("## 🖥️ 3. Jel Elektroforez Simülasyonu")
    gel_candidates = res_df[res_df['Uyumlu_Enzim_Sayısı'] > 0]
    
    if gel_candidates.empty:
        st.warning("Jel simülasyonu için uygun polimorfik kesim profili bulunamadı.")
    else:
        col_gel1, col_gel2 = st.columns([1, 2])
        with col_gel1:
            selected_snp_id = st.selectbox("Görüntülenecek SNP'yi Seçin:", gel_candidates['SNP_ID'].tolist())
            snp_data = gel_candidates[gel_candidates['SNP_ID'] == selected_snp_id].iloc[0]
            
            ref_allele = snp_data['REF']
            alt_allele = snp_data['ALT']
            pcr_size = snp_data['PCR_Boyutu_bp']
            enz_list = snp_data['Kullanılan_Enzim_Listesi']
            profile_str = snp_data['CAPS_Kesim_Profili']
            
            selected_snp_enz = st.selectbox("Kullanılacak Enzimi Seçin:", enz_list)
            
            frag_dict = parse_profile_for_gel(profile_str, selected_snp_enz, ref_allele, alt_allele, pcr_size)
            
            st.info(f"**SNP ID:** {selected_snp_id}\n\n"
                    f"**PCR Ürünü:** {pcr_size} bp\n\n"
                    f"**Allel 1 ({ref_allele}):** {frag_dict.get(ref_allele, [pcr_size])} bp\n\n"
                    f"**Allel 2 ({alt_allele}):** {frag_dict.get(alt_allele, [pcr_size])} bp")
                    
        with col_gel2:
            fig_gel = render_gel_simulator(frag_dict, ref_allele, alt_allele, pcr_size, selected_snp_enz)
            st.pyplot(fig_gel)

    st.divider()
    st.markdown("## 💾 4. Veri Dışa Aktarımı")
    st.dataframe(res_df.drop(columns=['Kullanılan_Enzim_Listesi']))
    
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        export_df = res_df.copy()
        export_df['Kullanılan_Enzim_Listesi'] = export_df['Kullanılan_Enzim_Listesi'].apply(lambda x: ", ".join(x))
        export_df.to_excel(writer, sheet_name='GBS2CAPS_Rapor', index=False)
    st.download_button("📂 Detaylı Raporu Excel Olarak İndir", output.getvalue(), "GBS2CAPS_Tasarimlar.xlsx")
