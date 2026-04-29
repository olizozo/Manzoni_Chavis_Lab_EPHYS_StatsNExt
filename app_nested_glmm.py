import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats import ttest_1samp, f_oneway
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# --- 1. CONFIGURATION ET STYLE ---
st.set_page_config(page_title="Chavis & Manzoni Lab - Master Stats", layout="wide")

col_l, col_r = st.columns([2, 5]) 
with col_l:
    try: st.image("logo_chavis_final.png", width=360) 
    except: st.info("Logo Chavis/Manzoni Lab")
with col_r:
    st.markdown("# 🧬 Master Stats Suite : Plasticité Synaptique")
    st.markdown("### Analyse de Population & Validation d'Induction")
    st.markdown("#### *Protocoles standardisés : GLMM, ANOVA & Corrections*")

st.divider()

# --- 2. GESTION DES GROUPES (SIDEBAR) ---
st.sidebar.header("⚙️ CONFIGURATION")
bin_size = st.sidebar.number_input("Bin d'alignement (min)", value=1.0, step=0.5)

if 'groups' not in st.session_state:
    st.session_state.groups = [{"name": "Group 1"}, {"name": "Group 2"}]

def add_group():
    new_id = len(st.session_state.groups) + 1
    st.session_state.groups.append({"name": f"Group {new_id}"})

st.sidebar.button("➕ Ajouter un Groupe", on_click=add_group)

all_data = []
colors = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e', '#8c564b', '#e377c2']

for idx, group in enumerate(st.session_state.groups):
    with st.sidebar.expander(f"📁 {group['name']}", expanded=True):
        group['name'] = st.text_input(f"Nom du groupe", value=group['name'], key=f"name_{idx}")
        files = st.file_uploader(f"Charger CSV ({group['name']})", type="csv", accept_multiple_files=True, key=f"files_{idx}")
        
        if files:
            st.markdown("**📝 Identifiants des Portées (Litters)**")
            for f in files:
                try:
                    animal_name = f.name.replace('.csv', '')
                    c1, c2 = st.columns([1, 1])
                    with c1: st.write(f"🐁 `{animal_name}`")
                    with c2: user_litter = st.text_input("Portée", value=animal_name, key=f"lit_{idx}_{animal_name}", label_visibility="collapsed")
                    
                    df = pd.read_csv(f)
                    df['Group'] = group['name']
                    df['Animal_ID'] = animal_name
                    df['Portee'] = user_litter.strip()
                    
                    t_col = next((c for c in ['Time_avg', 'Time_min'] if c in df.columns), df.columns[0])
                    df['Grid_Time'] = (np.round(df[t_col] / bin_size) * bin_size).round(2)
                    all_data.append(df)
                except Exception as e:
                    st.error(f"Erreur sur {f.name}: {e}")

# --- 3. ONGLETS D'ANALYSE ---
if all_data:
    master_df = pd.concat(all_data, ignore_index=True)
    amp_col = next((c for c in ['Amp_m', 'Amp_norm', 'Amplitude'] if c in master_df.columns), None)
    area_col = next((c for c in ['Area_m', 'Area_norm', 'Ch_m', 'Charge_norm'] if c in master_df.columns), None)
    
    for col in [amp_col, area_col]:
        if col is not None: master_df[col] = pd.to_numeric(master_df[col], errors='coerce')
    
    tab_plot, tab_glmm, tab_simple = st.tabs(["📉 Visualisation & Export", "🔬 Nested GLMM", "📊 Stats Simples (Époques)"])

    # --- TAB 1 : VISUALISATION & PRISM EXPORT ---
    with tab_plot:
        st.subheader("Visualisation de la Plasticité")
        err_mode = st.radio("Afficher la dispersion :", ["SEM", "95% CI"], horizontal=True)
        metrics_to_plot = [m for m in [amp_col, area_col] if m is not None]
        
        if metrics_to_plot:
            fig, axes = plt.subplots(len(metrics_to_plot), 1, figsize=(12, 5 * len(metrics_to_plot)), sharex=True)
            if len(metrics_to_plot) == 1: axes = [axes]
            for i, col in enumerate(metrics_to_plot):
                for g_idx, g_name in enumerate(master_df['Group'].unique()):
                    sub = master_df[(master_df['Group'] == g_name) & (master_df[col].notna())]
                    if sub.empty: continue
                    stats_df = sub.groupby('Grid_Time')[col].agg(['mean', 'std', 'count']).reset_index()
                    stats_df['sem'] = stats_df['std'] / np.sqrt(stats_df['count'].replace(0, 1))
                    err = stats_df['sem'] if err_mode == "SEM" else stats_df['sem'] * 1.96
                    axes[i].fill_between(stats_df['Grid_Time'], stats_df['mean']-err.fillna(0), stats_df['mean']+err.fillna(0), color=colors[g_idx % len(colors)], alpha=0.15)
                    axes[i].plot(stats_df['Grid_Time'], stats_df['mean'], 'o-', color=colors[g_idx % len(colors)], label=f"{g_name} (N={sub['Animal_ID'].nunique()})", markersize=4)
                axes[i].axhline(100, color='black', ls='--', alpha=0.4)
                axes[i].set_ylabel(f"{col} (%)")
                axes[i].legend(loc='upper right')
            st.pyplot(fig)

            # Export Wide pour Prism
            st.divider()
            export_df_list = []
            for g_name in master_df['Group'].unique():
                sub_df = master_df[master_df['Group'] == g_name]
                g_export = pd.DataFrame({'Grid_Time': sub_df['Grid_Time'].unique()}).sort_values('Grid_Time')
                g_export.insert(0, 'Group', g_name)
                for col_name, col_id in [("Amplitude", amp_col), ("Charge", area_col)]:
                    if col_id and col_id in sub_df.columns:
                        c_stats = sub_df.groupby('Grid_Time')[col_id].agg(['mean', 'std', 'count']).reset_index()
                        g_export = pd.merge(g_export, c_stats, on='Grid_Time', how='left')
                        g_export[f'{col_name}_Mean'] = g_export['mean']
                        g_export[f'{col_name}_SEM'] = g_export['std'] / np.sqrt(g_export['count'].replace(0, 1))
                        g_export[f'{col_name}_CI95'] = g_export[f'{col_name}_SEM'] * 1.96
                        g_export.drop(columns=['mean', 'std', 'count'], inplace=True)
                export_df_list.append(g_export)
            if export_df_list:
                full_export_df = pd.concat(export_df_list, ignore_index=True).round(2)
                st.download_button("📥 Télécharger Données Wide (Prism-Ready)", full_export_df.to_csv(index=False).encode('utf-8'), "plasticity_prism_ready.csv", "text/csv")

    # --- TAB 2 : NESTED GLMM AVEC INTERPRÉTATION ---
    with tab_glmm:
        st.subheader("🔬 Modèle Mixte Emboîté (Maintenance)")
        with st.expander("📖 Documentation Mathématique du GLMM", expanded=False):
            st.markdown(r"$$Y_{ijkl} = \beta_0 + \beta_1(Temps) + \beta_2(Groupe) + \beta_3(Temps \times Groupe) + u_k + v_{l(k)} + \epsilon_{ijkl}$$")
            st.write("Le modèle contrôle les effets de portée (Litter) et d'animal (Pseudo-réplication).")

        target = st.selectbox("Métrique GLMM", metrics_to_plot)
        t_start = st.number_input("Début analyse maintenance (min)", value=10.0)
        ref_group = st.selectbox("Référence Contrôle", master_df['Group'].unique())
        
        if st.button("🚀 EXÉCUTER LE NESTED GLMM"):
            df_stats = master_df[(master_df['Grid_Time'] >= t_start) & (master_df[target].notna())].copy()
            df_stats['Time_C'] = df_stats['Grid_Time'] - t_start
            df_stats['Group'] = pd.Categorical(df_stats['Group'], categories=[ref_group] + [g for g in df_stats['Group'].unique() if g != ref_group])
            try:
                model = smf.mixedlm(f"{target} ~ Time_C * C(Group)", df_stats, groups=df_stats["Portee"], vc_formula={"Animal": "0 + C(Animal_ID)"})
                result = model.fit()
                
                # Aide à l'interprétation
                st.markdown("### 📊 Interprétation Scientifique")
                p_values = result.pvalues
                st.info(f"**Comparaisons vs : {ref_group}**")
                
                interactions = [c for c in p_values.index if ":" in c]
                for inter in interactions:
                    other_g = inter.split('[T.')[1].split(']')[0]
                    p_inter = p_values[inter]
                    if p_inter < 0.05:
                        st.error(f"⚠️ **DIVERGENCE : {other_g} vs {ref_group} (p={p_inter:.3f})** : Les cinétiques diffèrent.")
                    else:
                        st.success(f"⚪ **STABILITÉ : {other_g} vs {ref_group} (p={p_inter:.3f})** : Les courbes sont parallèles.")
                
                st.table(result.summary().tables[1])
                st.download_button("📥 Télécharger Résultats GLMM", result.summary().tables[1].to_csv().encode('utf-8'), "glmm_results.csv", "text/csv")
            except Exception as e: st.error(f"Erreur GLMM : {e}")

    # --- TAB 3 : STATS SIMPLES (ÉPOQUES) AVEC ANOVA & CORRECTIONS ---
    with tab_simple:
        st.subheader("📊 Analyses par Époques (Avant/Après & Inter-Groupes)")
        
        with st.expander("📖 Documentation Mathématique : Époques & Corrections", expanded=False):
            st.markdown(r"""
            ### 1. Test de Plasticité (Benjamini-Hochberg FDR)
            Comparaison de chaque groupe à sa baseline ($100\%$) par un t-test : $t = \frac{\bar{x} - 100}{SEM}$.
            
            ### 2. Comparaison Inter-Groupes (ANOVA & Tukey HSD)
            L'ANOVA vérifie si une différence globale existe. Si $p < 0.05$ (Trigger), le test de **Tukey HSD** identifie les paires divergentes.
            """)

        c1, c2, c3 = st.columns(3)
        with c1: epoch1 = st.slider("Époque 1 (Baseline)", -30, 60, (-10, 0))
        with c2: epoch2 = st.slider("Époque 2 (Early)", 0, 60, (5, 15))
        with c3: epoch3 = st.slider("Époque 3 (Late)", 0, 60, (40, 50))
        target_s = st.selectbox("Métrique Stats Simples", metrics_to_plot)
        
        if st.button("📊 CALCULER & EXPORTER LE RAPPORT COMPLET"):
            simple_results_wide = []
            
            for name, (start, end) in [("Baseline", epoch1), ("Époque Early", epoch2), ("Époque Late", epoch3)]:
                st.markdown(f"### ⏱️ {name} [{start} à {end} min]")
                mask = (master_df['Grid_Time'] >= start) & (master_df['Grid_Time'] <= end)
                epoch_data = master_df[mask].groupby(['Group', 'Animal_ID'])[target_s].mean().reset_index()
                
                epoch_row = {"Époque": name, "Fenêtre": f"{start} à {end}"}
                p_raw = []; temp_stats = []; groups_vals = []; g_names = []
                
                for g in epoch_data['Group'].unique():
                    vals = epoch_data[epoch_data['Group'] == g][target_s].dropna()
                    if len(vals) > 1:
                        _, p = ttest_1samp(vals, 100)
                        p_raw.append(p); temp_stats.append((g, vals)); groups_vals.append(vals); g_names.append(g)

                # --- 1. TEST DE PLASTICITÉ VS BASELINE ---
                st.markdown("**1. Induction vs Baseline (FDR Benjamini-Hochberg)**")
                if p_raw:
                    _, p_corr, _, _ = multipletests(p_raw, method='fdr_bh')
                    for i, (g, v) in enumerate(temp_stats):
                        pc = p_corr[i]; mv = np.mean(v); sv = np.std(v, ddof=1)/np.sqrt(len(v))
                        res_txt = f"**{g}** : Moy={mv:.1f}%, p_BH={pc:.4f}"
                        if pc < 0.05: st.success(f"✅ {res_txt} (Significatif)")
                        else: st.info(f"⚪ {res_txt} (Non-significatif)")
                        epoch_row.update({f"{g}_Moy": mv, f"{g}_SEM": sv, f"{g}_p_BH": pc})

                # --- 2. ANOVA & TUKEY (INTER-GROUPES) ---
                st.markdown("**2. Comparaison Inter-Groupes (ANOVA & Tukey HSD)**")
                if len(groups_vals) > 1:
                    _, pa = f_oneway(*groups_vals)
                    epoch_row["ANOVA_p"] = pa; epoch_row["Groupes_Comp"] = " vs ".join(g_names)
                    if pa < 0.05:
                        st.error(f"⚠️ **ANOVA SIGNIFICATIVE (p={pa:.4f})**")
                        st.markdown("👉 **ACTION :** Consultez Tukey HSD ci-dessous.")
                        # Tukey
                        all_v = np.concatenate(groups_vals)
                        all_l = np.concatenate([[gn]*len(gv) for gn,gv in zip(g_names, groups_vals)])
                        tukey = pairwise_tukeyhsd(all_v, all_l)
                        st.text(tukey.summary())
                    else: st.success(f"⚪ **ANOVA NON-SIGNIFICATIVE (p={pa:.4f})**")
                
                simple_results_wide.append(epoch_row)
                st.divider()

            if simple_results_wide:
                st.download_button("📥 Télécharger Rapport Époques", pd.DataFrame(simple_results_wide).round(4).to_csv(index=False).encode('utf-8'), "rapport_stats_master.csv", "text/csv")

else:
    st.info("👈 Chargez vos fichiers CSV pour commencer.")