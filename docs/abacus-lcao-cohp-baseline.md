# ABACUS LCAO-COHP 项目基础认知

## 1. 方法定位

COHP, Crystal Orbital Hamilton Population, 是一种能量分辨的化学键分析方法。它将能带能量中的局域轨道对贡献分解到原子对或轨道对上，用于判断某一能量区间的电子态对特定键是成键、反键还是近似非键。

在局域轨道基下，Kohn-Sham 态可写为

```text
|psi_nk> = sum_mu c_mu,n(k) |phi_mu(k)>
```

能带能量满足

```text
epsilon_nk = sum_mu,nu c*_mu,n(k) H_mu,nu(k) c_nu,n(k)
```

因此，原子 I 与 J 之间的 LCAO-COHP 可形式化写为

```text
COHP_IJ(E) =
  sum_n,k w_k f_nk delta(E - epsilon_nk)
  sum_{i in I, j in J} Re[c*_{Ii,n}(k) H_{Ii,Jj}(k) c_{Jj,n}(k)]
```

常用绘图习惯是画 `-COHP`，使成键贡献显示在正方向。积分到费米能级的量通常称为 ICOHP，可作为键强弱的半定量指标。

## 2. 可解释性边界

COHP 分解的是 band-energy 相关项，而不是完整总能。它主要反映共价杂化和轨道相互作用，对离子静电、双计数修正、交换相关能、核-核相互作用等总能项没有直接分解意义。

COHP 与 Mulliken 布居分析类似，强依赖局域轨道基的选择。只要基函数的径向形状、截断半径、多 zeta 设置、极化轨道、正交化方式或投影空间发生变化，原子对/轨道对分解结果都可能明显改变。整体本征值由广义本征问题约束，但将 `c^dagger H c` 拆到具体轨道对后并不唯一。

因此，本项目中的直接 ABACUS LCAO-COHP 应理解为一种 ABACUS NAO 表示下的 COHP-like 分析量。它适合同一 ABACUS 版本、同一赝势、同一 NAO 基组和同一计算设置下做趋势比较，不应默认与 LOBSTER 的 pCOHP/ICOHP 数值逐点一致。

## 3. 常规 COHP 实践

当前社区中最成熟的 COHP 工作流通常由 LOBSTER 完成。典型路径是：

1. 使用平面波或 PAW DFT 程序完成 SCF，保留波函数。
2. LOBSTER 将平面波/PAW 波函数投影到一套化学分析用局域原子轨道基。
3. 检查投影质量，尤其是 charge spilling / absolute charge spilling。
4. 指定原子对、距离范围和是否 orbitalwise。
5. 读取 `COHPCAR`、`ICOHPLIST` 等输出，绘制 `-COHP` 并分析 ICOHP。

LOBSTER 的稳定性来自一整套投影定义和质量控制，而不仅是 COHP 公式本身。常见 basis set 包括 `bunge`、`koga`、`pbeVaspFit2015`。这些基组服务于化学解释，与 ABACUS LCAO 中为高效展开 Kohn-Sham 波函数而优化的数值原子轨道不是同一概念。

## 4. ABACUS LCAO 实现路线

本项目当前 Python 后处理脚本位于 `refs/cohp.py`，ABACUS 输出解析工具位于 `refs/read_abacus_out.py`。

直接 LCAO-COHP 路线需要 ABACUS 输出：

- `out_mat_hs 1`: `OUT.${suffix}/data-${ik}-H` 与 `OUT.${suffix}/data-${ik}-S`
- `out_wfc_lcao 1`: `OUT.${suffix}/WFC_NAO_K${ik}.txt` 或 `WFC_NAO_GAMMA${ik}.txt`
- k 点权重：优先读取 `OUT.${suffix}/kpoints`，缺失时可退化为等权
- 费米能级：优先读取 `OUT.${suffix}/running_scf.log`

计算流程为：

1. 读取每个 k 点的 `H(k)`、`S(k)`、本征矢 `C(k)`、本征值 `epsilon_nk`、占据数和 k 点权重。
2. 对目标原子轨道集合 `atomI_orbs` 和 `atomJ_orbs` 计算 `Re[c_i^* H_ij c_j]`。
3. 对轨道对求和，得到每个 `(n,k)` 上的原子对贡献。
4. 按 k 点权重汇总，并按本征能量聚合。
5. 对能量轴做 zero-padding 和可选 Gaussian smoothing。
6. 输出曲线数据和图像。

脚本同时保留 COOP 和 pCOHP-like 实验函数。pCOHP-like 路线用 `A^dagger C` 构造投影表示，但当前 `A` 仍主要来自 ABACUS NAO 重叠矩阵子块，因此不能等同于 LOBSTER 的标准外部投影基。

## 5. 已有尝试和结论

issue #3718 与本地 `refs/LCAO-COHP-implement-notes.md` 记录了两类尝试：

- 直接 ABACUS LCAO-COHP：使用 ABACUS NAO 基下的 `H(k)` 和 `C(k)` 直接计算 `c^* H c`。
- pCOHP-like 尝试：将 ABACUS LCAO 波函数再投影到一组选定轨道上，构造投影哈密顿量后计算类似 pCOHP 的量。

测试体系包括 diamond、GaAs、CsCl，并尝试 full basis 与 minimal basis。已有结论是：这些结果与 LOBSTER/pCOHP 的常规输出不一致。最合理解释不是单一代码错误，而是 COHP 表示本身对局域基/投影基强依赖。

## 6. 后续工程判断

若目标是实现 ABACUS 内部自洽的 LCAO-COHP 分析工具，当前路线是可行的，但文档和输出名称必须明确标记为 ABACUS-NAO-dependent COHP。

若目标是复现 LOBSTER 语义，需要额外建立：

- 标准化化学投影基，例如 STO/GTO-like 基；
- `A_{alpha,mu}=<phi_alpha|A_mu>` 的明确构造；
- spilling 或 band-overlap 类投影质量指标；
- 与 diamond、GaAs、CsCl、Na 等基准体系的逐步验证；
- 先验证 band/DOS/pDOS 投影重构质量，再比较 COHP/ICOHP。

项目当前应避免把“公式可计算”误解为“结果可与 LOBSTER 直接等价”。真正需要固化的是输入输出格式、轨道索引定义、基组依赖声明、以及可重复的端到端验证案例。

