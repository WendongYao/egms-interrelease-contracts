# 运行记录：Confirm01_JAG

| 时间（UTC） | 步骤 | 文件与哈希 |
|---|---|---|
| 2026-11-17 10:04 | 用户同意冻结；合同冻结 | CONTRACT_CN.md `61cd70e914f947d22864d5cbc8a49ef9d4196fcd2c99ba70efa64dd1a9fa672a`；草案移入 `_superseded/` |
| 2026-11-17 10:05 | N 集抽签（云端） | code/s80_draw.py `7b1f8199…14b2c02`；DataHoldings.txt（当日下载）`6b2c2dc7…7961d3`；tiles_145.txt `7080583c…b2c02b7`；输出 data/N_draw_order.csv、N_draw_counts.json（候选数 Al 85、S 200、C 206、Nt 140） |
| 2026-11-17 10:06 | 冻结参数（O 集，在读取任何 N 集数据之前） | code/s81_frozen_params.py `57097c7b…bafbcfe0`；data/FROZEN_PARAMS.json `f2a066e6…ef016548`。B0：密集 8.81、中间 17.18、稀疏 45.15 |
| 2026-11-17 10:10 | 下载清单 | DOWNLOAD_LIST_CN.md；本机已有文件清单 data/local_egms_zips.txt |

说明：抽签候选规则未按国家过滤（合同 §2.2 原文），不少入选瓦片位于 EGMS 覆盖范围外（土耳其、乌克兰、白俄罗斯、俄罗斯、格陵兰等），按 §2.2 的跳过规则处理，以 Explorer 实际是否存在为准。这是按合同执行，不是偏差。
0310fba07ab7a7ca366d12873b0d4502c5e565f2251c688dec98b9511206b919  data/N_draw_order.csv
| 2026-11-17 10:33 | 云端：GHSL R2023A BUILT-S 与 POP（2020 年，100 m，Mollweide）共 24 个 1000 km 瓦片开始下载；NGL steps.txt 与 N 集候选瓦片（名次 ≤ 45）内 856 个站点的 IGS20 tenv3 已缓存 | code/s83_ghsl_download.sh、code/s84_ngl_download.sh（云端）；GHSL 瓦片编号公式已用 R3_C19 的地理变换核对 |
| 2026-11-17 10:35 | 全瓦片处理脚本 | code/s82_fulltile.py `f95a0f5b8ed01ec294a0cfe878134daca94167e5f298bb2bf687b4c8202963da`；code/rc_common.py（自 RelCons01 原样复制）`d4fe15eba92999dc5a33fc82e37e11012fb08aa1400a29cc18d569aae56411cc`；code/s85_manifest_O.py `0adf3825c301ca6a8cea47a19a4c992abe348b99c096489a607d74c31f4a1ccd`；data/manifest_O.csv（136 个 O 瓦片，文件选取规则同 RelCons01 s52）`3d99463cce2014e988a1b1df12c659306bc32f1ac60f566a530347e08d848d37` |
| 2026-11-17 10:36 | 冒烟（合同 §7：≤ 2 个 O 集瓦片，只检查运行与输出形状）：E15N10、E50N47 | 运行正常，峰值内存约 1.2 GB，每瓦片 3–7 秒；输出形状与字段齐全；40 个随机点上向量化拟合与 rc_common.fit_primary 的最大差 < 2e-6 mm/yr（float32 舍入）。未计算、未打印任何统计量。冒烟输出 data/fulltile_smoke/ 不进入结果（正式运行写入 data/fulltile/） |
| 2026-11-17 10:32 | 下载程序探测（Windows，run.bat） | 存档版本只剩 2020-2024 与 2020-2025；2019-2023 已不可下载。按方案未下载任何文件 |
| 2026-11-17 10:43 | 修订前的结构检查：只下载 E15N10 的 2020-2025 版（tools/egms_download/check/，sha256 3a2c9510…adeb42），与本机 2020-2024_1 比较结构 | 2020-2025 的日期为 20201209–20251125（303 期），2020-2024 为 20200108–20241218（302 期），共有日期 246 个；坐标共有 20,792 个（21,723 / 21,274）；pid 不同；共有日期上逐格字符串相同的比例 1.2%，mean_velocity 相同的比例 3.1%。结论：2020-2025 是窗口后移约一年的独立重处理，不是续接。未计算任何速度差统计 |
| 2026-11-17 10:47–11:39 | Windows 下载（run.bat，修订 1 的版本对） | tools/egms_download/download_log.csv：跳过（存档中无此瓦片）27 个，全部在覆盖范围外或边境；下载失败 11 个文件，待重跑续传 |
| 2026-11-17 13:15 | 合同修订 1 冻结；合同偏差记录追加偏差 1 | 见 README_HASH.txt |
| 2026-11-17 13:12 | 说明 | 修订 1 与偏差 1 正文写的冻结时间 13:15 UTC 为预估；实际以 README_HASH.txt 的时间戳 13:11:38 UTC 为准 |
| 2026-11-17 13:20 | s82 改为分块计算（E29N48 有 86.8 万个匹配点，整块读入超出 3 GB 内存被系统终止） | code/s82_fulltile.py `58b50de4757d37adadcc39dd2e7cec3dd94a1baf8aff264b1a547ac378749951`；在两个冒烟瓦片上与旧版输出逐位相同；最大瓦片峰值内存 1.25 GB。fulltile_O 中已有的 9 个瓦片输出由旧版生成，与新版逐位相同，保留 |
| 2026-11-17 14:0x | O 集全瓦片提取完成（136 个）；P2 清单 | code/s91_manifest_P2.py `3f5eb0a34e041960ad946352e29930b108961d6d9adc6720bb6effc2a8d78b30`；N′ 齐全 65 个、待补 7 个；O′ 齐全 133 个、待补 3 个 |
| 2026-11-17 14:4x | 分析脚本写定并记录哈希（读取任何 N′、O′ 抽取结果之前；O 集的 A、B 部分结果也尚未计算） | code/CODE_HASHES_frozen_before_analysis.txt（s86 抽样提取、s87 GNSS、s88 第一部分、s89 第二部分、s90 第三部分、s92 站点表、s31 与 rc_common 原样复用）。聚合方式等合同未规定的实现细节写在各脚本开头的说明里（例如 A1、A2 以瓦片中位数为主、A3 以加总比为主） |
| 2026-11-17 14:4x | 冒烟（合同 §7）：s88、s89、s90 在 2 个 O 集瓦片（E15N10、E50N47）上运行 | 运行正常、输出字段齐全；s88 的测点标准差按坐标匹配率 100%；GHSL 权重可读出。未查看、未打印任何统计量，冒烟输出已删除 |
| 2026-11-17 14:5x | 补下载完成（用户 13:51 重跑）；s91 修正：重跑时已下载文件记为 use_local、路径在 raw_new | code/s91_manifest_P2.py `add1a6c15926828703216d8fc8951e10fec66cabd29daaa3773b40f4554f162f`；N′ 72 个（Al 12、S 20、C 20、Nt 20），O′ 136 个；与第一版清单相比原有条目路径无变化 |
| 2026-11-17 15:0x | s86 种子修正：抽样种子按最终清单的排序位置定义；补下载后清单多出 7 个瓦片，51 个已提取瓦片的排序位置改变，其抽样结果移到 data/_to_delete/sample_Nprime_oldseed/ 并按最终清单重提（尚未被任何分析读取） | — |
| 2026-11-17 15:0x–15:46 | 云端运行：s87（N′ 705 站中 328 站入选；O′ 576 站中 269 站入选）、s88（N′、O′）、s90（N′、O′）、s89（O、O′、N′） | 输出见 results/ 与 results/RESULT_FILE_HASHES.txt；s89 改动见合同偏差 2 |
