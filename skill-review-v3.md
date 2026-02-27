# skill-review-v3：12 个 chem skill 全量总结 + DRD2 闭环（含 docking / MLFF / experiment）复盘

时间：2026-02-26（UTC）  
覆盖 skills（12）：
- chem-literature
- chem-qsar
- chem-gnn
- chem-molgen
- chem-admet
- chem-retrosynthesis
- chem-rxn-conditions
- chem-llm
- chem-docking
- chem-3dgen
- chem-mlff
- chem-experiment

---

## (1) 12 个 skill：一句话核心规则（我需要“背下来”的那种）

1) **chem-literature**：所有“常见/罕见、SAR、先验”结论尽量给可核查证据与复现检索路径；API 429 时要有替代通道。

2) **chem-qsar**：默认 scaffold split；先 RF+Morgan baseline；报告 RMSE/MAE/R² 与 RMSE/std(y)；预处理只在 train 拟合。

3) **chem-gnn**：必须对照 RF baseline、报告参数量；GNN 只有在 scaffold 泛化显著提升才算赢，否则仅作一致性/不确定性参考。

4) **chem-molgen**：KL collapse 是默认；必须做 collapse 诊断；评估必须含 validity/uniqueness/novelty/diversity 与分布对照。

5) **chem-admet**：Rules-first（Lipinski/Veber/QED）+ traffic light 多端点 profile；用于过滤而不是事后解释。

6) **chem-retrosynthesis**：先官能团/断键逻辑再模板；每步 forward validate；输出必须带 SA Score，并按可合成性过滤。

7) **chem-rxn-conditions**：把 retro 的“纸面路线”变成“可操作条件”，并给产率区间与总体产率（区间乘法传播）。

8) **chem-llm**：LLM 做“解释与整合”，不做计算；所有定量声明必须 RDKit/脚本验证，并标注 LLM 产出为待验证。

9) **chem-docking**：从共晶结构起步；先 receptor/ligand/box 准备并 redocking（RMSD 过阈值）再 dock 新分子；只做相对比较并记录版本与输入。

10) **chem-3dgen**：只有当 3D 形状/口袋条件真正重要时才用；点云/3D 生成后的键指派必须 RDKit 验证；同样要报告 validity/uniqueness/novelty。

11) **chem-mlff**：MACE-OFF 用于有机分子几何/能量；必须对比 MMFF baseline；能量单位要明确（eV↔kcal/mol），能量差异常要触发排查。

12) **chem-experiment**：先计划、后计算；用 hard gates（SA/QED/一致性/motif 等）+ multi-objective scoring（活性/QED/SA/dock/一致性）做可执行决策，并记录每步用到的 skill 与原因。

---

## (2) DRD2 闭环 + docking/MLFF/experiment：做对的、犯错的、纠正后的改进

### 做对的（可迁移的正确做法）
- Scaffold split + RF baseline 先行，并把模型泛化好坏用于后续“信任策略”。
- ADMET rules-first 在筛选阶段执行，不让明显不合格分子进入最终候选。
- Retro 需要 forward validation，避免伪路线；并且 SA 必须进入筛选门槛。
- Docking 必须从共晶出发，先做 redocking 验证再 dock 候选（并落盘版本/输入）。
- 对“数值口径”保持可验证：RDKit 性质、Vina score、RMSD 计算方法必须明确。

### 犯错 → 被纠正 → 修正策略
- **SA>7 仍入 Top5**：我一开始把“有模板+forward”误当“容易合成”。修正：SA 设为硬门槛（默认 SA<=6）。
- **可疑 motif 未过滤**：RDKit sanitize/模板匹配不等价于化学可行。修正：增加 SMARTS denylist/红旗层。
- **RF vs GCN 分歧未当作不确定性**：修正：加入一致性阈值（|RF−GNN|<=1）并惩罚分歧。

### Docking 的额外教训
- RMSD 口径很关键：在 6CM4 上，`obrms` vs RDKit heavy-atom RMSD 会导致 pass/fail 翻转。
- “最小调参”应系统化：margin/num_modes/exhaustiveness/pH 做小网格，并用 20 modes 取 best RMSD。

### MLFF 的额外教训
- MACE 与 MMFF 的最优几何可能非常接近（RMSD 很小），但对方势能面下的能量会出现几 kcal/mol 级别偏好差，这本身就是“模型差异信号”。

---

## (3) 我当前可以独立完成的端到端工作流（文字流程图）

```text
[chem-experiment: plan + gates + score]
   │
   ├─▶ [chem-literature: 先验/骨架/SAR/风险]
   ├─▶ [chem-qsar: RF+Morgan baseline, scaffold split]
   ├─▶ [chem-gnn: GCN 对照; 不赢则只做一致性]
   ├─▶ [chem-molgen: SELFIES VAE 生成 + collapse诊断]
   ├─▶ [chem-admet: rules-first 过滤]
   ├─▶ [chem-retrosynthesis: 路线 + forward validate + SA]
   ├─▶ [chem-rxn-conditions: 每步条件 + 产率区间 + 总产率]
   ├─▶ [chem-docking: 共晶→redocking→dock候选（可选）]
   ├─▶ [chem-mlff: MACE refine/strain（可选，尤其用于 docking pose/构象能量）]
   ├─▶ [chem-3dgen: 3D 条件生成（可选，只有当口袋形状约束真需要）]
   └─▶ [chem-llm: 生成medchem叙述/报告，但所有定量声明必须验证]
```

---

## (4) 当前最薄弱的 2–3 个环节 + 改进计划（具体到可执行）

1) **真实合成可行性 beyond template/SA**
- 问题：模板+SA 仍可能给出实验上不靠谱的路线（保护基、兼容性、可购性）。
- 计划：加入“前体可得性/常见试剂启发式”与 2–3 步短路线搜索；对风险点打标签。

2) **docking/receptor prep 的鲁棒性与解释性**
- 问题：不同准备与 RMSD 口径会影响结论；目前仍偏“能跑通”而非“足够稳”。
- 计划：固化 RMSD 口径（heavy-atom + atom mapping 假设），并加可视化与额外 sanity checks。

3) **chem-literature 证据链（在 API 429/无 key 情况下）**
- 问题：只能给定性判断时，会削弱结论硬度。
- 计划：建立浏览器检索 + PubChem/ChEMBL 结构检索备份；对 actives 做子结构频率统计支撑“常见/不常见”。

---

## (5) 给陛下的建议：下一步加什么能力/做什么实验来进一步提升我

1) **把 DMTA 闭环固化成“cycle 目录规范”**：每次一个 cycle，自动生成 `plan.md / candidates.csv / top5.md / decisions.md`，便于迭代不重来。

2) **做一次跨靶点迁移测试**：换一个共晶更干净、配体更明确的 GPCR/kinase，验证 12-skill 工作流是否仍稳定。

3) **引入 3D generation 的对照实验**：在有高质量 pocket 的靶点上，对比 SELFIES VAE（2D）vs 3D gen（chem-3dgen）在 docking/redocking 对齐与 novelty 上是否真的增益。

4) **把“可疑 motif / 反应性风险”规则库产品化**：维护一个小型 SMARTS 规则集 + 红旗解释，作为所有生成/筛选的统一 gate。
