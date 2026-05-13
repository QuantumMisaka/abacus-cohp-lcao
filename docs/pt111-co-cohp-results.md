# Pt(111)-CO ABACUS LCAO-COHP 结果记录

## 环境与作业

- SAI skill/environment 检查确认：`module avail abacus` 包含 `abacus/LTSv3.10.1-sm70-auto`。
- Slurm 模板参考：`/opt/sbatch_examples/gpu_abacus.sbatch`。
- 实际作业设置：4V100 分区，单卡 GPU，`OMP_NUM_THREADS=8`，不使用 `mpirun`，`ks_solver cusolver`。
- ASE 环境：`/home/pku-jianghong/liuzhaoqing/.conda/envs/ase/bin/python`，ASE `3.26.0b1`。
- 赝势轨道库：`~/PP_ORB/PP`、`~/PP_ORB/ORB`、`~/PP_ORB/apns-orbitals-precision-v1`。

作业 ID：

- `399760`: nspin=1 compact-basis relax
- `399761`: nspin=2 compact-basis relax
- `399774`: nspin=1 precision-basis final SCF
- `399783`: nspin=2 precision-basis final SCF

## 结构与输入

- 结构：Pt(111) 2x2 slab，4 层，共 16 个 Pt；top 位吸附 CO，共 18 原子。
- 真空方向：ABACUS 晶胞 B 方向。
- 约束：底两层 Pt 固定，其余 Pt、C、O 放开。
- 优化基组：`~/PP_ORB/ORB`。
- final SCF/COHP 基组：`~/PP_ORB/apns-orbitals-precision-v1`。
- 关键参数：`force_thr_ev 0.03`、`relax_method bfgs_trad`、`smearing_method gaussian`、`smearing_sigma 0.004`、`mixing_beta 0.4`、`kspacing 0.14`。

## 输出位置

- 工作流脚本：`scripts/pt111_co_workflow.py`
- nspin=1:
  - relax: `runs/pt111_co_top/nspin1/relax`
  - final SCF: `runs/pt111_co_top/nspin1/final_scf`
  - COHP: `runs/pt111_co_top/nspin1/cohp`
- nspin=2:
  - relax: `runs/pt111_co_top/nspin2/relax`
  - final SCF: `runs/pt111_co_top/nspin2/final_scf`
  - COHP: `runs/pt111_co_top/nspin2/cohp`

## COHP 数值摘要

COHP 后处理使用 `src/cohp.py` 完成。脚本从 final SCF 的 `OUT.ABACUS` 中读取 `data-*-H/S`、`WFC_NAO_K*.txt`、`kpoints` 和 `running_scf.log`，并根据 `mapping.json` 中记录的 ABACUS 全局 NAO 轨道范围选择 top Pt 与 C 的轨道集合。总 Pt-C 曲线使用 top Pt 的全部 NAO 与 C 的全部 NAO；分通道曲线进一步选择 Pt-d/C-p 或 Pt-d/C-s 轨道子集。输出包括能量分辨 `.dat` 文件和采用 `-COHP` 习惯绘制的 `.png` 文件。

| Setting | E_Fermi/eV | Pt-C/A | C-O/A | Channel | Spin | -ICOHP |
|---|---:|---:|---:|---|---|---:|
| nspin=1 | 1.768901 | 1.8479 | 1.1500 | top Pt-C total | sum | 0.058247 |
| nspin=1 | 1.768901 | 1.8479 | 1.1500 | top Pt-d / C-p | sum | 0.038627 |
| nspin=1 | 1.768901 | 1.8479 | 1.1500 | top Pt-d / C-s | sum | 0.013701 |
| nspin=2 | 1.768847 | 1.8479 | 1.1500 | top Pt-C total | up | 0.058249 |
| nspin=2 | 1.768847 | 1.8479 | 1.1500 | top Pt-C total | down | 0.058249 |
| nspin=2 | 1.768847 | 1.8479 | 1.1500 | top Pt-C total | sum | 0.116498 |
| nspin=2 | 1.768847 | 1.8479 | 1.1500 | top Pt-d / C-p | up | 0.038625 |
| nspin=2 | 1.768847 | 1.8479 | 1.1500 | top Pt-d / C-p | down | 0.038625 |
| nspin=2 | 1.768847 | 1.8479 | 1.1500 | top Pt-d / C-p | sum | 0.077250 |
| nspin=2 | 1.768847 | 1.8479 | 1.1500 | top Pt-d / C-s | up | 0.013700 |
| nspin=2 | 1.768847 | 1.8479 | 1.1500 | top Pt-d / C-s | down | 0.013700 |
| nspin=2 | 1.768847 | 1.8479 | 1.1500 | top Pt-d / C-s | sum | 0.027400 |

说明：绘图输出为 `-COHP`；数值表中的 `-ICOHP` 是 occupied energy range 上的 `-int COHP(E)dE`。当前 COHP 脚本遵循项目既有定义，不显式乘占据数。因此 nspin=2 的 up/down 单通道值与 nspin=1 更适合直接比较；nspin=2 的 `sum` 是两个 spin 通道相加。

## 科学结论

- nspin=1 和 nspin=2 的 top Pt-C 几何几乎相同：Pt-C 约 1.848 A，C-O 约 1.150 A。
- nspin=2 final SCF 的总磁矩最终降至约 `-7e-5 Bohr mag/cell`，说明该 Pt(111)-CO top 构型在当前设置下没有稳定自旋极化。
- nspin=2 的 up/down COHP 几乎完全相同；与 nspin=1 的 top Pt-C `-ICOHP` 也一致到 1e-5 量级。因此自旋设置没有改变 Pt-C 成键图像。
- top Pt-C 的主要 occupied bonding 贡献来自 Pt-d/C-p 通道：nspin=1 中 `0.038627 / 0.058247 ~= 66%`。Pt-d/C-s 贡献约 `24%`。这与 CO 在 Pt top 位吸附时 Pt d 态与 CO 2pi*/5sigma 相关前线轨道相互作用主导的图像相符，但这里应表述为 ABACUS-NAO COHP 下的定性结论。
- 除 top Pt-C 外，最近的其它 Pt-C 距离约 3.46 A，明显大于成键距离，因此本构型中可解释的 Pt-C COHP 主要集中在正下方 top Pt-C 对。

## 验证

- `python -m py_compile src/cohp.py src/read_abacus_out.py scripts/pt111_co_workflow.py`
- `python src/read_abacus_out.py`
- ABACUS final SCF 均生成 `data-*-H`、`data-*-S`、`WFC_NAO_K*.txt`、`kpoints`、`running_scf.log`。
