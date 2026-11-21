# Confirm04_Robustness 运行记录

| 时间 (UTC) | 事件 | 文件 / 哈希 |
|---|---|---|
| 2026-11-21 06:45 | 用户回复"同意冻结"；合同冻结 | CONTRACT_CN.md `c57fbc7ed6e4aa01f11cbbc65ddbb5cde83af61fc573a712b01a876ebb640a58` |

| 2026-11-21 06:50–06:54 | 冒烟（每集前 2 个瓦片）：只检查输出字段与形状；R3-0（T5、T10 处与 Confirm03 A_tile 逐位相同）在冒烟瓦片上为 True；需判读的条目共 23 条；未打印、未查看统计量；冒烟输出已删除 | — |
| 2026-11-21 06:54 | 代码冻结，云端全量运行 | code/s99_confirm04.py `59b0bbba5d33c341b83c7d543d2fffe3fbeb5e012953e58eeae4e999e0f0b376` |
| 2026-11-21 06:54–07:05 | 全量运行完成；结果写入 RESULT_CN.md | results/，哈希见 results/RESULT_FILE_HASHES.txt |

## 实现说明（运行前写入代码文档字符串）

1. 见 s99 文档字符串：Harrell–Davis 权重与抖动随机数的种子规则；块自助中每次抽取的块数等于该统计量所用数据的块数，用原点估计加块区间代入原判据；R3 的瓦片中位流失只用有移动点的瓦片；R4 中 mean_velocity_std 缺失的移动点单独计数、不放入缓冲区。

## 偏差

无。
