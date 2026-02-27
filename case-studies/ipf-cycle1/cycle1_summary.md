# IPF / ALK5 (TGFBR1) Cycle 1 — summary

目标：先把 ALK5/TGFBR1 的 **DMTA baseline**（数据→QSAR→生成→规则过滤→docking→重排）跑通；优先级为 **安全/可成药性（慢性病） > 活性**。

---

## (1) 关键数据指标

### Data (ChEMBL)
- Target: **CHEMBL260 (TGFBR1/ALK5)**
- Raw activity records iterated: **16324**
- Unique molecules (median pChEMBL aggregated): **5067**
- File: `data/chembl_alk5_activity.csv`

### QSAR / GNN (scaffold split)
- Split: N=5067 (train=4053, valid=507, test=507)
- RF (Morgan, 500 trees): RMSE **0.586**, R² **0.567**
- GCN (3-layer, hidden=128): RMSE **0.821**, R² **0.150**（弱于 RF）
- File: `reports/model_metrics.md`

### Generation (SELFIES VAE)
- Trained 8 epochs；**KL collapse 明显**（AU≈0/32）
- Generated valid unique SMILES: **4936**
- File: `generated/generated_raw.csv`

### Candidate pool（rules-first）
- Rules-first ADMET（Lipinski+Veber + QED>0.3）+ novelty 后 pool：**N=2999**
- With single-step template routes + forward validation：**1689**
- Files: `reports/candidate_pool.csv`, `reports/routes_5.md`

### Docking QC（chem-docking）
- PDB: **1VJY**；co-crystal ligand resname: **460**
- Redocking: **passed**（best heavy-atom RMSD ≈ **0.05 Å**）
- Co-crystal redock mode1 score ≈ **-10.23 kcal/mol**
- Docking 扩展：Top300（实际成功打分 **297** 个）
  - `reports/top300_dock_scores.csv`
  - 重排评分表：`reports/top300_scored.csv`

---

## (2) 发现的问题 & 解决方案（本轮落地）

### 问题 A：生成分布偏离 kinase 抑制剂化学空间（比“筛选”更致命）
- 证据：对候选池做 hinge-binder 常见骨架覆盖率（heuristic SMARTS）统计，candidate_pool 里典型 kinase-like motif 基本接近 0，而训练集覆盖显著。
- 输出：`reports/scaffold_coverage.md` / `reports/scaffold_coverage.json`

### 问题 B：慢性病（IPF）场景下不希望出现 N–N / hydrazone / hydrazine 类片段
- 处理：将 N–N / hydrazone-like 作为 **hard denylist**（直接剔除，不降权）。

### 问题 C：docking 扩展规模时易因环境依赖导致超时/中断
- 处理：Top300 docking 拆分为 chunk 运行并合并；同时对 docking 脚本增加“缺少 meeko CLI 时的 OpenBabel pdbqt fallback”，提高可复现性。

---

## (3) Cycle 2 建议（重点：生成器改进方向）

本轮最关键结论：**Cycle 2 应优先改生成策略，把分子分布拉回 kinase/ALK5 常见化学空间**，否则再怎么扩 docking/筛选都只能在“错的空间里精挑细选”。

建议按优先级：
1. **加先验的生成**：用训练集/文献先验的 hinge-binder scaffold 做条件生成（或 scaffold-constraint），至少保证候选含合理的 hinge binder。
2. **改 VAE anti-collapse**：引入 KL anneal / cyclical annealing / free-bits / word dropout 等，先把 AU 提起来（避免 collapse）。
3. **改分子表示或模型家族**：考虑 graph-based / JT-VAE / fragment-based / transformer + constraint（SELFIES 虽保证有效性，但不保证“药物化学合理性”）。
4. **训练数据与目标对齐**：将 ALK5 活性数据中“高频 scaffold”做显式统计并反馈到生成端（prompt/conditioning/采样过滤）。

建议 Cycle 2 的最小可验证里程碑：
- hinge-binder motif 覆盖率（candidate vs train）不再接近 0；
- Top300 docking 的 best score 能接近 co-crystal（Δ 更小），而不是普遍 Δ>2–3。
