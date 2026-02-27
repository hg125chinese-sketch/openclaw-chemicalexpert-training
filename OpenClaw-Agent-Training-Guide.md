### 12.9 Cross-Target Migration Test: IPF/ALK5 DMTA Cycle 1

> **目标**：验证 ChemicalExpert 12-skill 训练成果能否迁移到全新靶点——不给提示，让 CE 自主完成从选靶到出候选的完整闭环。

#### 场景设定

用户给出开放式需求："找抗纤维化(IPF)药物候选物，吡非尼酮/尼达尼布疗效有限。"

CE 需要自主完成：选靶点 → 数据获取 → 模型训练 → 分子生成 → ADMET/docking/retro → 综合评分 → 质量诊断。

#### CE 的自主决策

| 决策点 | CE 的选择 | 理由 |
|--------|-----------|------|
| 靶点 | ALK5/TGFBR1 (CHEMBL260) | TGF-β核心通路，ChEMBL数据充足，有共晶结构 |
| 权重策略 | Safety-first (QED 0.25, SA 0.25, activity 0.20) | 慢性病用药，安全性优先于活性 |
| Docking QC | Redocking RMSD < 2Å | 标准 docking 验证阈值 |

**关键观察**：CE 在请求用户确认前自行列出了 3 个需要决策的点（靶点选择、权重优先级、docking 阈值），展示了 chem-experiment skill 的"plan before compute"原则内化。

#### Cycle 1 执行结果

| 指标 | 数值 | 备注 |
|------|------|------|
| ChEMBL actives | ~1200 | ALK5 IC50 数据 |
| 生成候选池 | ~3000 | SELFIES VAE |
| ADMET gates 通过 | 2999 | Lipinski + Veber + QED |
| Redocking RMSD | 0.052 Å | PDB 1VJY，极优 |
| Top5 Vina range | -2.78 ~ -7.58 | 初始 Top5（仅 docked 5） |
| Co-crystal Vina | -10.23 | |

#### CE 自主诊断的问题（最有价值的部分）

CE 在完成初始 Top5 后**主动指出 4 个问题**，不是被用户追问：

1. **Docking 竞争力不足**：Top5 最好 -7.58 vs 共晶 -10.23，Δ > 2.5 kcal/mol
2. **Motif 审查遗漏**：Top1 含 hydrazone-like N-N 键 → 慢性病代谢不稳定
3. **Docking 覆盖率过低**：2999 池中只 docked 5 → 统计不足
4. **根因诊断**：VAE KL collapse → 生成分布偏离 kinase 化学空间 → 生成策略问题 > 筛选问题

> **训练意义**：这 4 条诊断直接映射训练中建立的核心原则：
> - 问题 1 → chem-docking "always compare to reference"
> - 问题 2 → chem-admet "chronic disease → flag metabolically unstable motifs"
> - 问题 3 → chem-experiment "统计显著性需要覆盖"
> - 问题 4 → chem-molgen "KL collapse 是默认失败模式"

#### Cycle 1 补完

用户根据 CE 诊断下达 3 条指令，CE 全部执行完成：

1. **N-N/hydrazone 加入 hard denylist** → 2 条 SMARTS，旧 Top1/4/5 均被拦截 ✅
2. **Docking 扩到 Top 300** → 分 6 chunk，297/300 成功 ✅
3. **ALK5 骨架覆盖分析** → 证实诊断：pyrazole train 23.4% vs pool 0.6%

补完后 Top5 最好 Vina score: -9.591（vs 初始 -7.58）。

#### 工程经验

| 问题 | 解决方案 |
|------|----------|
| OAuth token 过期 → OpenClaw cooldown 循环 | onboard 重新认证 + 重启 gateway |
| meeko receptor prep 不稳定 | fallback 到 OpenBabel 直接写 pdbqt |
| 长任务 gateway timeout | 分 chunk 执行（Top300 → 6×50） |

#### 迁移测试结论

| 维度 | 评价 |
|------|------|
| 自主规划 | ✅ 独立选靶、设权重、定 QC 标准 |
| 技能调用 | ✅ 正确编排 12 个 skill 的子集 |
| 自我诊断 | ✅ 主动发现 4 个质量问题，根因准确 |
| 执行韧性 | ✅ timeout/中断后自主恢复 |
| 待改进 | ⚠️ 生成器化学空间覆盖（需 Cycle 2 验证） |

**结论**：12-skill 训练成功迁移到新靶点。最有价值的是自主质量诊断能力——CE 自己发现并分析了候选质量不足的根因。
