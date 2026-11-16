# RQAOA 應用研究（2026 全國智慧運算與量子資訊創新應用競賽）

> 最新結果報告（網頁版）：https://claude.ai/artifact/4Vbfze9FsDE7vfMcmBF7C1 ，離線版在 `docs/rqaoa_report.html`。
>
> 量子方法適用性檢查的子專案在另一個 repo：https://github.com/johnny-0120/quantum-fit-check

## 2026-11 新增：二元噴漆換色與補救實驗

- `research/rqaoa_opportunity/`：噴漆換色（Binary Paint Shop）RQAOA、參數轉移、RSG 等古典基準。
- `research/rqaoa_opportunity/fixes.py`、`adaptive.py`：RQAOA 輸給禁忌搜尋的原因診斷，以及多路徑、QRR、回溯、半路交棒等補救實測（結果在 `results/fixes.csv`、`results/adaptive.csv`）。
- `research/rqaoa_maxkcut/`：三色 RQAOA 與 Wi-Fi 頻道規劃實測。
- `research/rqaoa_opportunity/sources/`：引用 V. Vijendran 的 MIT 授權程式（XQAOA-BPSP），保留原作者授權標頭。

全部是 p=1 的 CPU 精確模擬，沒有上量子真機。

---

# Recursive QAOA 實測專案

2026-11-07，Johnny 的第 9 題研究。這是本機 CPU 模擬研究，沒有呼叫付費 API 或 IBM 真機。

## 目前主交付：第 9 題圖結構研究

- `output/第9題研究報告_圖結構與Recursive_QAOA.md` 與 `output/pdf/第9題研究報告_圖結構與Recursive_QAOA.pdf`：以第 9 題為主線的新報告。
- 本輪 208 個主圖案例：八類圖 88 個、固定 16 點／32 邊／四正則的配對改接 120 個。另有 72 筆深度設定、45 筆小圖 Aer 雜訊設定，分別涉及 8 張、3 張圖，不能當作額外獨立圖加總。
- `research/q9_results/main.csv`、`graphs_and_traces.json`：全部指標、邊清單與三個 R-QAOA 版本軌跡。
- `research/q9_results/analysis.json`：配對 bootstrap、失敗步驟及參數網格檢查。
- `research/q9_protocol.json`：本輪跑前範圍；研究參考過先前 pilot，並非正式預先註冊。

重跑主研究（先依下方安裝鎖定套件；新增 CVXPY／SCS 作數值 GW 基準）：

```powershell
$env:PYTHONIOENCODING='utf-8'
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
.\.venv\Scripts\python.exe research\q9_study.py
.\.venv\Scripts\python.exe research\q9_depth_noise.py
.\.venv\Scripts\python.exe research\q9_depth_noise.py noise
.\.venv\Scripts\python.exe research\q9_analyze.py
.\.venv\Scripts\python.exe research\q9_report.py
```

會覆寫 `research/q9_results/` 與本輪報告，不會改寫舊 `results/benchmark.csv`。初次 Powell 深度試跑另存為 `depth_powell_pilot.csv`，後續暖啟動 L-BFGS-B 的方法修訂有獨立紀錄。先前的損失預測提案只保留作歷史候選，沒有成為本輪主方法。

## 先前 pilot 報告與原始實驗

- `output/pdf/Recursive_QAOA_實測報告.pdf`：白話解釋、實測圖表、失敗案例及參賽定位。
- `output/Recursive_QAOA_實測報告.md`：同內容可編輯版。
- `results/benchmark.csv`：98 個主實驗案例。
- `results/traces.json`：每一步縮圖、參數、相關性及最後分組。
- `results/shot_sensitivity.csv`：180 次有限量測／讀出擾動補測。
- `results/scale.csv`：30、50 點的 4 個二分圖案例。
- `results/cutoff_sensitivity.csv`：45 次終止大小補測。
- `results/failure_audit.json`：4 個失敗案例首次排除所有最佳解的步驟。
- `results/validation.json`、`results/audit.json`：公式、線路、縮圖與輸出驗證。

## 安裝與重跑（PowerShell）

```powershell
uv venv --python 3.11 .venv
uv pip install --python .venv\Scripts\python.exe -r requirements-lock.txt
$env:PYTHONIOENCODING='utf-8'
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
.\.venv\Scripts\python.exe rqaoa_experiment.py verify
.\.venv\Scripts\python.exe rqaoa_experiment.py run
.\.venv\Scripts\python.exe rqaoa_followup.py
.\.venv\Scripts\python.exe audit_results.py
.\.venv\Scripts\python.exe make_report.py
```

重跑會覆寫本專案 `results/` 與報告輸出，原始簡報、照片及雲端／Obsidian 筆記不會修改。報告製作使用 Windows 微軟正黑體。

## 方法與限制

QAOA 固定 p=1。以解析的淺層電路相關性優化參數，Qiskit Statevector 核對分布；Aer 另驗證 4,096 shots。R-QAOA 每次取全部點對中絕對相關性最大的點對，固定同／異組關係，消去一個變數，到 4 點後窮舉，再還原分組。縮圖後保留負權重，所有最終分數回原圖計算。

主實驗的 R-QAOA 用精確相關性。有限量測補測的參數仍以理想期望值訓練，只在消去決策時加入抽樣誤差和可選的獨立 3% 讀出翻轉；不是完整硬體噪音或端到端含 shot 訓練。標準 QAOA 的 256 shots 以精確 statevector 分布抽樣 20 批，不是 IBM job。原始圖皆為自建無權重圖，不是 Team 9 原始圖／程式的逐項重現。

R-QAOA 對照含：隨機縮圖、最重邊縮圖、32 起點古典單點翻轉。最後抽樣數加倍／對齊的欄位只對齊最後量測數，不能視為整體量子資源公平比較。未宣稱速度優勢、創新演算法或 p>=2 的結果。

文獻：Bravyi 等 arXiv:1910.08980；Ozaeta 等 arXiv:2012.03421；Bae 等 arXiv:2408.13207。本程式自行實作，沒有直接複製論文作者程式。`sources/provenance.json` 保存核讀來源及版本雜湊；Graphify 的專案索引過期，已回讀最新原文。
