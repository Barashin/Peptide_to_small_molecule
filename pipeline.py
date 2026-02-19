"""
pipeline.py
===========
Peptide → Small Molecule 変換パイプライン メインスクリプト

使い方:
    venv/bin/python pipeline.py [PDB_FILE] [OPTIONS]

オプション:
    --protein-chain    A    タンパク質チェーン ID (default: A)
    --peptide-chain    B    ペプチドチェーン ID (default: B)
    --top-residues     3    設計に使う上位残基数 (default: 3)
    --cutoff           4.5  接触距離カットオフ Å (default: 4.5)
    --output-dir       .    出力ディレクトリ (default: results)
    --exhaustiveness   8    smina exhaustiveness (default: 8)
    --num-modes        9    smina num_modes (default: 9)
    --skip-docking          ドッキングをスキップ
    --skip-rescore          smina rescoring をスキップ
    --skip-sasa             ΔSASA解析をスキップ
    --skip-prodigy          PRODIGY予測をスキップ
"""

import argparse
import os
import json
import sys

from analyze_interactions  import load_structure, find_contacts, score_peptide_residues, summarize
from extract_pharmacophore import extract_pharmacophore, print_pharmacophore, \
                                   save_pharmacophore_csv, save_pharmacophore_pml
from design_small_molecule import select_key_residues, build_peptidomimetic, \
                                   print_candidates, save_candidates_sdf, calculate_drug_likeness
from visualize             import plot_residue_scores, plot_interaction_map, plot_pharmacophore_3d
from dock_with_smina       import (extract_chain, get_peptide_bbox,
                                   run_smina_per_molecule, print_docking_summary,
                                   save_docking_results, plot_docking_scores,
                                   rescore_all_docked, plot_score_decomposition,
                                   save_rescore_csv)
from analyze_sasa          import (calc_delta_sasa, print_sasa_summary,
                                   save_sasa_csv, plot_sasa_comparison,
                                   score_residues_by_sasa)
from analyze_prodigy       import (run_prodigy, print_prodigy_summary,
                                   save_prodigy_json, plot_prodigy_contacts)


def parse_args():
    p = argparse.ArgumentParser(description="Peptide to Small Molecule Pipeline")
    p.add_argument("pdb", nargs="?", default="Protein_Peptide.pdb")
    p.add_argument("--protein-chain",    default="A")
    p.add_argument("--peptide-chain",    default="B")
    p.add_argument("--top-residues",     type=int,   default=3)
    p.add_argument("--cutoff",           type=float, default=4.5)
    p.add_argument("--output-dir",       default="results")
    p.add_argument("--exhaustiveness",   type=int,   default=8)
    p.add_argument("--num-modes",        type=int,   default=9)
    p.add_argument("--skip-docking",     action="store_true")
    p.add_argument("--skip-rescore",     action="store_true")
    p.add_argument("--skip-sasa",        action="store_true")
    p.add_argument("--skip-prodigy",     action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print("  Peptide → Small Molecule Pipeline")
    print(f"{'='*60}")
    print(f"  PDB          : {args.pdb}")
    print(f"  タンパク質鎖 : {args.protein_chain}")
    print(f"  ペプチド鎖   : {args.peptide_chain}")
    print(f"  上位残基数   : {args.top_residues}")
    print(f"  接触カットオフ: {args.cutoff} Å")
    print(f"  出力先       : {args.output_dir}/\n")

    # ──────────────────────────────────────────
    # STEP 1: 構造読み込み & 相互作用解析
    # ──────────────────────────────────────────
    print("[STEP 1] 相互作用解析...")
    model    = load_structure(args.pdb)
    contacts = find_contacts(model, args.protein_chain, args.peptide_chain, args.cutoff)
    result   = summarize(contacts)
    scores   = result["residue_scores"]

    # 接触データをJSONに保存 (numpy型をPython nativeに変換)
    def to_native(v):
        if hasattr(v, "item"):   # numpy scalar
            return v.item()
        if hasattr(v, "tolist"): # numpy array
            return v.tolist()
        return v

    contacts_json = []
    for c in contacts:
        cj = {k: to_native(v) for k, v in c.items()
              if k not in ("peptide_coord", "protein_coord")}
        cj["peptide_coord"] = c["peptide_coord"].tolist()
        cj["protein_coord"] = c["protein_coord"].tolist()
        contacts_json.append(cj)
    with open(os.path.join(args.output_dir, "contacts.json"), "w") as f:
        json.dump(contacts_json, f, indent=2)
    print(f"  → {args.output_dir}/contacts.json 保存完了")

    # ──────────────────────────────────────────
    # STEP 1b: ΔSASA 解析
    # ──────────────────────────────────────────
    sasa_result = None
    sasa_scores = {}
    if not args.skip_sasa:
        print("\n[STEP 1b] ΔSASA 解析...")
        sasa_result = calc_delta_sasa(args.pdb,
                                      chain_protein=args.protein_chain,
                                      chain_peptide=args.peptide_chain)
        print_sasa_summary(sasa_result)
        save_sasa_csv(sasa_result,
                      os.path.join(args.output_dir, "sasa_analysis.csv"))
        plot_sasa_comparison(sasa_result,
                             os.path.join(args.output_dir, "sasa_comparison.png"))
        sasa_scores = score_residues_by_sasa(sasa_result["delta_peptide"])
    else:
        print("\n[STEP 1b] ΔSASA をスキップ (--skip-sasa)")

    # ──────────────────────────────────────────
    # STEP 1c: PRODIGY 結合親和性予測
    # ──────────────────────────────────────────
    prodigy_result = None
    if not args.skip_prodigy:
        print("\n[STEP 1c] PRODIGY 結合親和性予測...")
        prodigy_result = run_prodigy(args.pdb,
                                     chain_protein=args.protein_chain,
                                     chain_peptide=args.peptide_chain,
                                     sasa_result=sasa_result)
        print_prodigy_summary(prodigy_result)
        save_prodigy_json(prodigy_result,
                          os.path.join(args.output_dir, "prodigy_result.json"))
        plot_prodigy_contacts(prodigy_result,
                              os.path.join(args.output_dir, "prodigy_contacts.png"))
    else:
        print("\n[STEP 1c] PRODIGY をスキップ (--skip-prodigy)")

    # ──────────────────────────────────────────
    # スコア統合: 距離ベース + ΔSASA の加重平均
    # ──────────────────────────────────────────
    if sasa_scores:
        print("\n  [スコア統合] 距離ベース × ΔSASA を統合...")
        # 両スコアをそれぞれ最大値で正規化して合算 (重み 0.6 : 0.4)
        max_dist  = max(scores.values()) if scores else 1
        max_sasa  = max(sasa_scores.values()) if sasa_scores else 1
        combined  = {}
        all_keys  = set(scores.keys()) | set(sasa_scores.keys())
        for key in all_keys:
            d = scores.get(key, 0) / max_dist
            s = sasa_scores.get(key, 0) / max_sasa
            combined[key] = round(0.6 * d + 0.4 * s, 4)
        combined_scores = dict(sorted(combined.items(), key=lambda x: -x[1]))
        print("  統合スコア (上位5残基):")
        for (num, name), sc in list(combined_scores.items())[:5]:
            print(f"    {name}{num}: {sc:.3f}")
        # 以降の設計では統合スコアを使用
        scores = combined_scores

    # ──────────────────────────────────────────
    # STEP 2: ファーマコフォア抽出
    # ──────────────────────────────────────────
    print("\n[STEP 2] ファーマコフォア抽出...")
    key_resnums = [num for (num, _), score in scores.items() if score > 0]
    features    = extract_pharmacophore(model, args.peptide_chain, key_residues=key_resnums)
    print_pharmacophore(features)
    save_pharmacophore_csv(features, os.path.join(args.output_dir, "pharmacophore.csv"))
    save_pharmacophore_pml(features, os.path.join(args.output_dir, "pharmacophore.pml"))

    # ──────────────────────────────────────────
    # STEP 3: 低分子設計
    # ──────────────────────────────────────────
    print("\n[STEP 3] 低分子設計...")
    selected   = select_key_residues(scores, top_n=args.top_residues)
    print(f"  選択残基: {[(n, name, f'{sc:.1f}pt') for n, name, sc in selected]}")
    candidates = build_peptidomimetic(selected)
    print_candidates(candidates)

    sdf_path = os.path.join(args.output_dir, "candidate_ligands.sdf")
    save_candidates_sdf(candidates, sdf_path)

    # ──────────────────────────────────────────
    # STEP 4: 可視化
    # ──────────────────────────────────────────
    print("\n[STEP 4] 可視化...")
    try:
        plot_residue_scores(scores,
                            os.path.join(args.output_dir, "residue_scores.png"))
        plot_interaction_map(contacts,
                             os.path.join(args.output_dir, "interaction_map.png"))
        plot_pharmacophore_3d(features,
                              os.path.join(args.output_dir, "pharmacophore_3d.png"))
    except Exception as e:
        print(f"  可視化エラー (スキップ): {e}")

    # ──────────────────────────────────────────
    # STEP 5: サマリーレポート
    # ──────────────────────────────────────────
    print("\n[STEP 5] サマリーレポート生成...")
    report = {
        "pdb"              : args.pdb,
        "protein_chain"    : args.protein_chain,
        "peptide_chain"    : args.peptide_chain,
        "total_contacts"   : len(contacts),
        "contact_types"    : dict(result["type_counts"]),
        "key_residues"     : [{"resnum": n, "resname": name, "score": sc}
                               for n, name, sc in selected],
        "pharmacophore_count": len(features),
        "candidates": [
            {"name": c["name"], "smiles": c["smiles"],
             **calculate_drug_likeness(c["mol"])}
            for c in candidates
        ],
        "sasa": {
            "total_buried_sasa": sasa_result["total_buried_sasa"],
            "interface_peptide_residues": sasa_result["interface_residues_peptide"],
        } if sasa_result else None,
        "prodigy": {
            "dG_kcal_mol" : prodigy_result["dG_kcal_mol"],
            "Kd"          : prodigy_result["Kd_str"],
            "total_ic"    : prodigy_result["total_ic"],
        } if prodigy_result else None,
    }
    report_path = os.path.join(args.output_dir, "report.json")
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    # ──────────────────────────────────────────
    # STEP 6: smina ドッキング
    # ──────────────────────────────────────────
    dock_results = []
    if not args.skip_docking:
        print("\n[STEP 6] smina ドッキング...")
        dock_dir     = os.path.join(args.output_dir, "docking")
        os.makedirs(dock_dir, exist_ok=True)

        receptor_pdb = os.path.join(dock_dir, "receptor.pdb")
        peptide_ref  = os.path.join(dock_dir, "peptide_ref.pdb")
        extract_chain(args.pdb, args.protein_chain, receptor_pdb)
        extract_chain(args.pdb, args.peptide_chain,  peptide_ref)

        box = get_peptide_bbox(args.pdb, args.peptide_chain)
        print(f"  ドッキングボックス: center({box['center_x']:.1f}, "
              f"{box['center_y']:.1f}, {box['center_z']:.1f})  "
              f"size({box['size_x']:.0f} x {box['size_y']:.0f} x {box['size_z']:.0f} Å)")

        dock_results = run_smina_per_molecule(
            receptor_pdb, sdf_path, peptide_ref,
            dock_dir,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
        )
        print_docking_summary(dock_results)
        save_docking_results(dock_results,
                             os.path.join(dock_dir, "docking_results.csv"))
        plot_docking_scores(dock_results,
                            os.path.join(dock_dir, "docking_scores.png"))
        with open(os.path.join(dock_dir, "docking_results.json"), "w") as f:
            json.dump(dock_results, f, indent=2, ensure_ascii=False)

        # ── rescoring (score_only) ──
        if not args.skip_rescore:
            print("\n  [STEP 6b] smina score_only rescoring...")
            enriched = rescore_all_docked(dock_results, receptor_pdb, dock_dir)
            save_rescore_csv(enriched,
                             os.path.join(dock_dir, "rescore_terms.csv"))
            plot_score_decomposition(enriched,
                                     os.path.join(dock_dir, "score_decomposition.png"))
            dock_results = enriched
        else:
            print("  [STEP 6b] rescoring をスキップ (--skip-rescore)")

        # report にドッキング結果を追加
        report["docking"] = [
            {"name": r["name"], "best_score": r["best_score"]}
            for r in dock_results
        ]
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
    else:
        print("\n[STEP 6] ドッキングをスキップ (--skip-docking)")

    print(f"\n{'='*60}")
    print("  パイプライン完了!")
    print(f"  出力ファイル ({args.output_dir}/):")
    for fname in sorted(os.listdir(args.output_dir)):
        fpath = os.path.join(args.output_dir, fname)
        if os.path.isfile(fpath):
            size = os.path.getsize(fpath)
            print(f"    {fname:<35} ({size:>8,} bytes)")
        elif os.path.isdir(fpath):
            print(f"    {fname}/")
            for sub in sorted(os.listdir(fpath)):
                spath = os.path.join(fpath, sub)
                if os.path.isfile(spath):
                    size = os.path.getsize(spath)
                    print(f"      {sub:<33} ({size:>8,} bytes)")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
