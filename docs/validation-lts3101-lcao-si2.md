# ABACUS LTS 3.10.x LCAO-COHP 端到端验证记录

## 环境

- 项目 submodule: `abacus-develop` at `v3.10.1-7-g1dccf29ca`
- 实际运行二进制: `/home/pku-jianghong/liuzhaoqing/apps/abacus/3.10-new-sai/bin/abacus`
- 二进制报告版本: `ABACUS v3.10.0`, commit `9300a8657`
- 测试目录: `runs/lts3101_lcao_si2`
- 测试体系: `abacus-develop/examples/scf/lcao_Si2`

说明：当前工作区没有由 submodule 源码直接构建出的 `abacus` 可执行文件，因此实际 SCF 使用本机已有 3.10.x LTS 二进制完成。submodule 源码中 `out_mat_hs` 和 `out_wfc_lcao` 的输出命名与本次实际输出一致：`data-${ik}-H/S` 与 `WFC_NAO_K${ik}.txt`。

## 输入改动

在示例 `INPUT` 基础上增加：

```text
out_mat_hs 1 8
out_wfc_lcao 1
out_app_flag 1
```

当前仓库复现实例已将 `pseudo_dir`、`orbital_dir` 改为 `../data/legacy-si` 相对路径，并在 `examples/data/legacy-si` 中保存该验证实际使用的 Si 赝势和轨道文件。

## ABACUS 输出

SCF 正常收敛，最终输出目录为 `runs/lts3101_lcao_si2/OUT.ABACUS`。关键文件包括：

```text
data-0-H ... data-7-H
data-0-S ... data-7-S
WFC_NAO_K1.txt ... WFC_NAO_K8.txt
kpoints
running_scf.log
istate.info
```

解析检查结果：

```text
nks = 8
H0/S0 shape = (26, 26)
C0 shape = (26, 14)
bands = 14
sum(k weights) ~= 1.0001
E_Fermi = 7.111283804 eV
```

## COHP 后处理命令

本项目的后处理入口是 `src/cohp.py`。它通过 `src/read_abacus_out.py` 读取 ABACUS LCAO SCF 输出目录中的 `data-*-H/S`、`WFC_NAO_K*.txt`、`kpoints` 和 `running_scf.log`，因此前置 SCF 必须打开 `out_mat_hs 1 8`、`out_wfc_lcao 1` 和 `out_app_flag 1`。后处理阶段需要手动给出两个原子或轨道组对应的 ABACUS 全局 NAO 编号。

Si2 示例中每个 Si 原子有 13 个 NAO，因此使用 0-12 与 13-25 作为 Si-Si 原子对轨道集合：

```bash
env MPLBACKEND=Agg python src/cohp.py \
  --out-dir runs/lts3101_lcao_si2/OUT.ABACUS \
  --atom-i-orbs 0,1,2,3,4,5,6,7,8,9,10,11,12 \
  --atom-j-orbs 13,14,15,16,17,18,19,20,21,22,23,24,25 \
  --method COHP \
  --de 0.1 \
  --smooth-nstddev 3 \
  --invert \
  --output-prefix runs/lts3101_lcao_si2/si_si_COHP
```

生成文件：

```text
runs/lts3101_lcao_si2/si_si_COHP.dat
runs/lts3101_lcao_si2/si_si_COHP.png
```

数据检查：

```text
points = 376
energy range = [-5.8703505520, 31.6296494480] eV
COHP range = [-0.0585821615, 0.1533233220]
integral up to E_Fermi = -0.0619846133
```

## 结论

当前 Python 脚本在补充 ABACUS 3.10.x 文本波函数格式兼容后，可以从实际 LCAO SCF 输出完成端到端 COHP 后处理，并输出曲线数据和图片。

需要注意：该验证确认的是 ABACUS LCAO/NAO 表示下的 COHP-like 后处理链路可运行；它不构成与 LOBSTER pCOHP 数值一致性的证明。
