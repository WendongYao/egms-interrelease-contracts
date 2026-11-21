# 运行记录：Confirm05_Example

| 时间（UTC） | 步骤 | 文件与哈希 |
|---|---|---|
| 2026-11-21 08:28 | 用户同意冻结；合同冻结 | CONTRACT_CN.md `fa8f0f9cfc1955297920fce2e1734946a3a6e793a5fcae26c3f38158695f0446`（见 README_HASH.txt） |
| 2026-11-21 08:30 | 写脚本；冒烟（两瓦片字段与形状；analyze 流程在合成 OSM 与临时目录上走通，输出已删除，未读数）；记录代码哈希 | code/CODE_HASHES_frozen_before_analysis.txt：s100_osm_fetch.py `7c3482da…c498`，s101_example.py `385b3b3b…bcbc` |
| 2026-11-21 08:31 | 选址（§1，云端 s101 select） | results/selection.json `b111c5bd…39e8`：瓦片 E39N31（N′），匹配点 353,140，GHSL 人口和 4,975,259（第二名 E55N17 3,669,171）；208 个瓦片无缺失栅格；窗口中心 E 3,942,750 / N 3,137,850（约 4.58°E, 51.23°N） |
| 2026-11-21 08:33 | OSM 下载（虚拟机 s100，Overpass，查询固定 2026-11-21T00:00:00Z） | data/osm/E39N31_overpass.json，19,113,885 字节（比合同预估的"几 MB"大），SHA-256 `f7154618…4c39`，16,190 个要素；查询文本 SHA-256 `c3b22a0f…76d2`；见 E39N31_meta.json |
| 2026-11-21 08:35 | 云端 s101 analyze（代码哈希复核通过） | results/example_counts.csv、summary.json、Fig_example_T5.png、Fig_example_T10.png；RESULT_FILE_HASHES.txt。实现自检：整瓦片 T5 C0 的移动点数与流失数与 Confirm03 A_tile.csv 中该瓦片逐位相同 |
| 2026-11-21 08:40 | 结果文件回传本机并复核哈希 | 直接回传的两张 PNG 字节被改写（像素逐一相同、哈希不同），改用 tar 回传后与 RESULT_FILE_HASHES.txt 全部一致；tar 留在 data/_xfer/ |

## 偏差

- 2026-11-21 10:00 UTC（呈现修订，不改任何数字）：论文版 `Fig06_antwerp_T5` 与补充图 `FigS_antwerp_T10`（脚本 `manuscript/scripts/build_all_figs_20261121.py`）重绘时调整了图层顺序（灰色无告警底图点 `zorder=1`、`s=0.35–0.45` 置于底图之下；OSM 主干铁路与高速公路线 `zorder=2`、`lw=0.65` 置于灰点之上；双版告警与翻转点 `zorder=3–4` 置于最上层），并将图幅改为 `7.1 × 4.95 in`、字号 `>= 7 pt`。两张图绘制前已对 `example_counts.csv` 全部区域的点数及 C0/C1 翻转计数做逐一断言校验（一致），不改变任何统计读数或结论。

## 实现说明

1. 城市名取 OSM 中离窗口中心最近的 place=city 节点：Antwerp（Antwerpen），距窗口中心约 12.6 km。
2. 主干铁路要素 4,851 条、高速公路要素 4,235 条（OSM way，含分段）。
3. 图中铁路、高速两列的缓冲区点（灰色）覆盖了线要素，线要素只在点稀疏处可见；论文版图可另画（只改呈现，不改数字），届时记为偏差。
4. 百分位计算用 Confirm03 A_tile.csv 中 O′ ∪ N′ 的 205 个瓦片（另 3 个瓦片 T5 无移动点，按定义不计）。
