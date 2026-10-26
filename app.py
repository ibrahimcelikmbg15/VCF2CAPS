import streamlit as st
import pandas as pd
import io
import tempfile
import os
import shutil
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio import Restriction
from collections import Counter
import matplotlib.pyplot as plt

st.set_page_config(page_title="GBS2CAPS: TaqMan & CAPS Marker Designer", layout="wide")

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
            
            frag_ref, frag_alt = get_fragments(c_ref, len(amplicon_ref)), get_fragments(c_alt, len(amplicon_alt))
            profiles.append(f"{enz.__name__} -> Allele {ref}: [{frag_ref}] | Allele {alt}: [{frag_alt}]")
            enz_list.append(enz.__name__)
            
    return " || ".join(profiles) if profiles else "No Suitable/Selective Enzyme Found", len(profiles), enz_list

def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def calculate_tm(seq):
    return (4 * (seq.count('G') + seq.count('C'))) + (2 * (seq.count('A') + seq.count('T')))

def design_marker_pure_python(left_flank, ref, alt, right_flank, snp_id, p_size, p_tm, tm_diff, p_gc, p_prod, f_len):
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
        
        snp_amp_pos = len(left_flank) - best_left
        probe_start, probe_end = max(0, snp_amp_pos - 10), min(len(amplicon_ref), snp_amp_pos + 11)
        
        cut_profile, enz_count, enz_list = analyze_amplicon(amplicon_ref, amplicon_alt, ref, alt)
        return (fwd_seq, rev_seq, round(fwd_tm, 1), round(rev_tm, 1), prod_size, amplicon_ref[probe_start:probe_end], amplicon_alt[probe_start:probe_end], formatli_snp, cut_profile, enz_count, enz_list)
        
    return ("Primer Not Found", "Primer Not Found", 0, 0, 0, "-", "-", formatli_snp, "No Suitable Amplicon", 0, [])

def plot_physical_map(df, map_type, top_enzymes=None):
    fig, ax = plt.subplots(figsize=(12, min(8, len(df['CHROM'].unique()) * 0.8)))
    
    if map_type == "Successful PCR and Enzymes Only":
        plot_df = df[(df['Forward_Primer'] != "Primer Not Found") & (df['Forward_Primer'] != "Not Found in FASTA") & (df['Compatible_Enzyme_Count'] > 0)]
    elif map_type == "Top 3 Most Frequent Enzymes" and top_enzymes:
        top_3 = [e[0] for e in top_enzymes[:3]]
        plot_df = df[df['Used_Enzyme_List'].apply(lambda x: any(e in x for e in top_3))]
    else:
        plot_df = df[df['FASTA_Header'] != "-"]

    if plot_df.empty:
        st.warning("No SNPs found matching this filter.")
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
        ax.plot(chrom_df['POS'] / 1_000_000, [y_pos]*len(chrom_df), '|', color='darkred', markersize=20, markeredgewidth=2, zorder=2)

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("Physical Position (Mb)")
    ax.set_title(f"Chromosome Physical Map: {map_type}", fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.tight_layout()
    return fig

def main():
    st.title("✂️ GBS2CAPS & TaqMan: Advanced PCR Marker Design")
    st.markdown("Design specific PCR primers directly from VCF files. Supports `.fa.gz` and `.fasta.gz` compressed reference genomes for memory-efficient processing.")
    
    col1, col2 = st.columns(2)
    with col1: vcf_file = st.file_uploader("1. VCF File (.vcf)", type=['vcf'])
    with col2: fasta_file = st.file_uploader("2. Reference Genome (.fasta, .fa, .gz)", type=['fasta', 'fa', 'fna', 'gz'])

    st.sidebar.header("🌡️ Marker & Probe Settings")
    flank_len = st.sidebar.number_input("FASTA Flank Window (bp):", min_value=150, max_value=1000, value=250, step=50)
    f_len = st.sidebar.slider("Formatted SNP Sequence Length (bp):", 21, 201, 51, step=10)
    p_prod = st.sidebar.slider("Target PCR Size (bp):", 50, 600, (100, 300), step=10)
    p_size = st.sidebar.slider("Ideal Primer Length (nt):", 15, 30, 20)
    p_tm = st.sidebar.slider("Ideal Melting Temp (Tm °C):", 50.0, 72.0, 60.0, step=0.5)
    tm_diff = st.sidebar.number_input("Maximum Tm Difference (°C):", min_value=0.1, max_value=10.0, value=2.0, step=0.1)
    p_gc = st.sidebar.slider("GC Content Range (%):", 20.0, 80.0, (40.0, 60.0))
    
    if vcf_file and fasta_file:
        if st.button("🚀 Design Markers and Generate Report"):
            df = parse_vcf(vcf_file)
            df = df[~df['ALT'].str.contains(',', na=False)].drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']).reset_index(drop=True)
            
            if len(df) == 0: return st.warning("No suitable Biallelic SNP found in this VCF file!")
            
            file_extension = ".fa.gz" if fasta_file.name.endswith(".gz") else ".fasta"
            with tempfile.NamedTemporaryFile(delete=False, suffix=file_extension) as tmp_fasta:
                shutil.copyfileobj(fasta_file, tmp_fasta)
                fasta_path = tmp_fasta.name
            
            try:
                fasta_obj = Fasta(fasta_path)
                fasta_keys = list(fasta_obj.keys())
            except Exception as e:
                return st.error(f"Error indexing FASTA: {str(e)}")

            results, all_enzymes = [], []
            progress_bar = st.progress(0)
            
            for i, row in df.iterrows():
                chrom, pos, ref, alt, snp_id = str(row['CHROM']), int(row['POS']), str(row['REF']), str(row['ALT']), str(row['ID'])
                actual_chrom = find_fasta_chrom(chrom, fasta_keys)
                
                if actual_chrom:
                    pos_0 = pos - 1
                    left_flank = str(fasta_obj[actual_chrom][max(0, pos_0 - flank_len):pos_0].seq).upper()
                    right_flank = str(fasta_obj[actual_chrom][pos_0 + len(ref):pos_0 + len(ref) + flank_len].seq).upper()
                    
                    f_seq, r_seq, f_tm, r_tm, p_size_val, p_ref, p_alt, form_snp, c_prof, enz_c, enz_l = design_marker_pure_python(
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
                    
            fasta_obj.close()
            os.remove(fasta_path)
            if os.path.exists(fasta_path + ".fai"): os.remove(fasta_path + ".fai")
            if os.path.exists(fasta_path + ".gzi"): os.remove(fasta_path + ".gzi")
            
            result_df = pd.DataFrame(results)
            
            st.divider()
            st.header("📊 GBS2CAPS Analysis Summary")
            
            total_snp = len(result_df)
            mapped_snp = len(result_df[result_df['FASTA_Header'] != "-"])
            primer_success = len(result_df[(result_df['Forward_Primer'] != "Primer Not Found") & (result_df['Forward_Primer'] != "Not Found in FASTA")])
            caps_success = len(result_df[result_df['Compatible_Enzyme_Count'] > 0])
            multi_enz = len(result_df[result_df['Compatible_Enzyme_Count'] > 1])
            
            col_a, col_b, col_c, col_d = st.columns(4)
            col_a.metric("Total VCF SNPs (Biallelic)", total_snp)
            col_b.metric("FASTA Matched & Probe Designed", mapped_snp)
            col_c.metric("Successful PCR Primers", primer_success)
            col_d.metric("SNPs with CAPS Enzyme", caps_success)
            
            top_enzymes = Counter(all_enzymes).most_common(10)
            
            st.markdown(f"**🔬 Enzyme & Restriction Details:** A total of **{len(set(all_enzymes))} different restriction enzymes** produced polymorphic cuts. **{multi_enz} SNPs** can be cut by more than one enzyme type.")
            if top_enzymes:
                st.markdown("**🏆 Top 10 Most Frequent Enzymes:** " + ", ".join([f"{e[0]} ({e[1]} times)" for e in top_enzymes]))
            
            st.divider()
            st.header("🧬 Chromosome Physical Map (MapChart View)")
            map_option = st.radio("Select Data to Visualize:", ["All SNPs", "Successful PCR and Enzymes Only", "Top 3 Most Frequent Enzymes"], horizontal=True)
            
            fig = plot_physical_map(result_df, map_option, top_enzymes)
            st.pyplot(fig)
            
            st.divider()
            st.success("All data processed. You can review the detailed table below and download it.")
            st.dataframe(result_df.drop(columns=['Used_Enzyme_List']))
            
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                export_df = result_df.copy()
                export_df['Used_Enzyme_List'] = export_df['Used_Enzyme_List'].apply(lambda x: ", ".join(x))
                export_df.to_excel(writer, sheet_name='Marker_Designs', index=False)
            st.download_button("📂 Download Full Report as Excel", output.getvalue(), "GBS2CAPS_Analysis_Report.xlsx")

if __name__ == "__main__":
    main()