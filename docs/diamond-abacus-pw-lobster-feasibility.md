# ABACUS PW 输出用于 LOBSTER COHP 的可行性测试

## 测试对象

- 体系：金刚石 primitive cell，2 个 C 原子。
- ABACUS：LTS v3.10.1，PW 基组，`ecutwfc 100` Ry，`nbands 16`，`KPT 7 7 7`，`symmetry -1`。
- ABACUS 关键输出参数：`out_wfc_pw 2`、`out_chg 1`、`out_band 1`。
- ABACUS 赝势：APNS `C.upf`，该文件为 norm-conserving UPF，`is_paw="F"`。
- LOBSTER：`/home/pku-jianghong/liuzhaoqing/work/sidereus/abacus-cohp/tools/lobster-5.1.1/lobster-5.1.1`。

## 依据文件

- ABACUS 波函数输出文档：`abacus-develop/docs/advanced/elec_properties/wfc.md`。
- ABACUS INPUT 参数文档：`abacus-develop/docs/advanced/input_files/input-main.md`。
- ABACUS PW 波函数二进制写出源码：`abacus-develop/source/module_io/write_wfc_pw.cpp`。
- LOBSTER Generic 格式说明：`tools/lobster-5.1.1/Generic/README.md`。
- LOBSTER 用户手册：`tools/lobster-5.1.1/Lobster_Users_Guide_5.1.1.pdf`。
- LOBSTER 方法论文页面：<https://www.oqi.ox.ac.uk/publication/1114472/europe-pubmed-central>，其摘要说明 LOBSTER 的投影方法基于 PAW DFT 计算。
- Quantum ESPRESSO 官方赝势页：<https://www.quantum-espresso.org/pseudopotentials/>，确认 QE 本身支持 NC、US 和 PAW 多类赝势。
- Quantum ESPRESSO UPF 说明：<https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format>，确认 UPF 容器可存储 NC、US、PAW 等不同类型赝势。

## 已完成输出检查

- ABACUS SCF 收敛：`True`。
- `WAVEFUNC*.dat` 数量：`343`。
- `istate.info`：`True`。
- `SPIN1_CHG.cube`：`True`。
- Generic manifest：`{"fft_grid": [32, 32, 32], "n_kpoints": 343, "status": "generated", "warning": "PAW data were taken from the LOBSTER QE diamond example, not from the APNS NC pseudopotential used by ABACUS. This is a format diagnostic only."}`。

## LOBSTER 对 PW 赝势类型的限制

需要区分两个问题：

1. QE 作为 PW-DFT 程序支持多类赝势。QE 官方赝势页明确列出 PAW、USPP 和 NC PP；UPF 官方说明也表明 UPF 是一个可容纳 NC、US、PAW 等类型的统一格式。因此，“LOBSTER 支持 QE 输出”并不等价于“LOBSTER 支持任意 QE 赝势类型输出”。
2. LOBSTER 对其可读取的 PW 结果有更强限制。LOBSTER 用户手册在 VASP 准备部分明确要求使用 PAW potential，而不是 ultrasoft pseudopotential；ABINIT 部分说明 VASP 的要求同样适用，并进一步限定为 PAW-XML 数据集；QE 部分说明 VASP 和 ABINIT 部分的要求同样适用于 QE，并在 FAQ 中要求读取本次计算使用的 PAW 数据 UPF 文件。

因此，对 LOBSTER 的标准 VASP/ABINIT/QE 接口而言，PW 计算结果需要来自 PAW 赝势/PAW 数据集。QE 虽然能用 NC、US、PAW 计算，但 LOBSTER 的 QE 接口要求的是 QE+PAW 数据，而不是 QE+NC 或 QE+US 的一般 PW 输出。

这一点也和 LOBSTER 的方法论文一致：其局域轨道投影方法建立在从 PAW DFT 计算重构化学信息的框架上。对本项目而言，APNS `C.upf` 的 `PP_HEADER` 为 `pseudo_type="NC"` 且 `is_paw="F"`，不含 LOBSTER Generic/PAW 重构所需的 augmentation charges、projectors、all-electron partial waves 等 PAW 数据。因此，当前 ABACUS PW+APNS NC 输出即使具备 PW coefficients，也缺少 LOBSTER 所需的 PAW 赝势侧信息。

## 直接读取 ABACUS PW 输出

LOBSTER 原生示例和二进制字符串显示其可自动识别 VASP、Quantum ESPRESSO、ABINIT 和 Generic 输入，但没有 ABACUS 原生检测路径。本测试将 ABACUS 的 `INPUT/STRU/KPT/WAVEFUNC*.dat/istate.info` 放入 LOBSTER 目录后运行。

结果文件存在：`True`。

错误摘要：

```text
LOBSTER v5.1.1 (g++ 9.3.0)
Copyright (C) 2024 by Chair of Solid-State and Quantum Chemistry, RWTH Aachen.
All rights reserved. Contributions by S. Maintz, V. L. Deringer, M. Esser, R. Nelson, C. Ertural, P. C. Mueller, M. Pauls, L. Sann, D. Schnieders, A. L. Tchougreeff, and R. Dronskowski
starting on host login-01.mr-sai.ai on 2026-05-14 at 18:28:06 CST using 16 threads
detecting used PAW program...
ERROR: could not determine which program was used for the quantum-chemical calculation.
ERROR: Make sure all required files are in the working directory and are not empty.
```

结论：ABACUS PW 输出不能像 VASP 的 `WAVECAR/vasprun.xml` 那样被 LOBSTER 直接消费。

## Generic 转换测试

转换脚本从 ABACUS `WAVEFUNC*.dat` 解析 k 点、PW 数量、G 网格索引、band 数、ecut 和 complex128 波函数系数，并写入 LOBSTER Generic 的 `LOBSTER_Kpoints/kPoint*.hdf5`：

- `Miller`: `(nPW, 3) int32`
- `PWCoeffs`: `(nBands, nPW) complex128`

能级和占据来自 `istate.info`，FFT grid 来自 `SPIN1_CHG.cube`，cell/atomic structure 来自 `STRU`。

重要限制：APNS `C.upf` 是 NC 赝势，不包含 LOBSTER Generic 文档要求的 PAW augmentation、projectors、all-electron wavefunctions。当前 Generic 诊断包为了测试格式通路，使用了 LOBSTER 自带 QE diamond 示例中的 C PAW UPF 填充 `paw` 字段。这与 ABACUS SCF 的赝势不一致，因此即使 LOBSTER 成功输出 COHP，也只能证明接口格式可能打通，不能作为严格科学结果。

Generic LOBSTER 结果文件存在：`True`。

错误摘要：

```text
LOBSTER v5.1.1 (g++ 9.3.0)
Copyright (C) 2024 by Chair of Solid-State and Quantum Chemistry, RWTH Aachen.
All rights reserved. Contributions by S. Maintz, V. L. Deringer, M. Esser, R. Nelson, C. Ertural, P. C. Mueller, M. Pauls, L. Sann, D. Schnieders, A. L. Tchougreeff, and R. Dronskowski
starting on host login-01.mr-sai.ai on 2026-05-14 at 18:33:33 CST using 16 threads
detecting used PAW program...
ERROR: could not determine which program was used for the quantum-chemical calculation.
ERROR: Make sure all required files are in the working directory and are not empty.
```

额外校验：将 LOBSTER 发行包自带的 Generic 示例原样复制到本测试目录并运行，本机 `lobster-5.1.1` 同样停在程序识别阶段，未打开 `LobsterInput.json`。这说明当前二进制/运行方式下 Generic detector 没有被触发；该现象独立于 ABACUS 转换数据本身。

发行包 Generic 示例摘要：

```text
LOBSTER v5.1.1 (g++ 9.3.0)
Copyright (C) 2024 by Chair of Solid-State and Quantum Chemistry, RWTH Aachen.
All rights reserved. Contributions by S. Maintz, V. L. Deringer, M. Esser, R. Nelson, C. Ertural, P. C. Mueller, M. Pauls, L. Sann, D. Schnieders, A. L. Tchougreeff, and R. Dronskowski
starting on host login-01.mr-sai.ai on 2026-05-14 at 18:30:50 CST using 24 threads
detecting used PAW program...
ERROR: could not determine which program was used for the quantum-chemical calculation.
ERROR: Make sure all required files are in the working directory and are not empty.
```

## 结论

1. **直接答案：不能直接做。** LOBSTER 不原生识别 ABACUS PW 输出。
2. **工程路径：理论上可以通过 LOBSTER Generic 做。** ABACUS `out_wfc_pw 2` 给出了构造 `PWCoeffs`/`Miller` 的核心数据，`out_band 1` 给出 `istate.info`，`out_chg 1` 可辅助获得 FFT grid。
3. **科学阻塞：PAW 数据一致性。** 严格 LOBSTER COHP 需要与 PW 波函数同源的 PAW projectors、augmentation charges、AE/PS partial waves。APNS 当前用于 ABACUS PW 的 C.upf 是 NC，不满足这一要求。
4. **推荐实践：** 当前生产级 ABACUS COHP 仍应使用本项目的 LCAO-COHP 后处理；ABACUS PW 到 LOBSTER 需要新增稳定的 Generic 导出器，并配套 ABACUS 可用且与 LOBSTER 数据结构一致的 PAW 势库。
