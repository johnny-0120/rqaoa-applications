# 跑「量子 vs 古典消去訊號」實驗（signals.py）

同一套 RQAOA 流程，只換「每一步綁哪一對」的判斷訊號，比較誰綁錯得少。
測資是開放邊界的二維方格自旋玻璃（平面圖），每一步用窮舉精確檢查最佳解是否還在。

## 1. 安裝（另一台電腦第一次跑）

需要 Python 3.11。在 repo 根目錄（PowerShell）：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
```

## 2. 先跑快速檢查（約 5–10 分鐘）

```powershell
cd research\rqaoa_opportunity
$env:PYTHONIOENCODING='utf-8'
..\..\.venv\Scripts\python.exe signals.py quick gauss
```

每跑完一個格子會印一行：`gauss 4 4 2000 {'random': (1, 3.33), ...}`，括號裡是（有沒有綁錯, 和最佳解差多少）。

## 3. 正式實驗

```powershell
..\..\.venv\Scripts\python.exe signals.py full gauss
```

- 4×4 方格 30 個（QAOA p=1、2、3），4×5 方格 20 個（QAOA p=1、2）。
- **預估 6–8 小時**，大部分時間花在 4×5 的 QAOA p=2。電腦不能睡眠或關機。
- 可以中斷後續跑：結果逐列寫入 `results/signals_gauss_full.csv`，重新執行同一行指令會跳過已完成的部分。
- 時間夠的話再跑 ±1 耦合版本：`signals.py full pm`（±1 的最佳解大量重複，綁錯較少，差異可能較難看出）。

## 4. 跑完之後

把 `results/signals_gauss_full.csv`（以及 `signals_pm_full.csv`，如果有跑）commit 並 push：

```powershell
git add research/rqaoa_opportunity/results/signals_*.csv
git commit -m "Add signal comparison results"
git push
```

**只加這些 CSV，不要用 `git add -A`**，避免把本機的其他檔案一起推上去。

## 欄位說明

| 欄位 | 意思 |
|---|---|
| signal | 判斷訊號：random、local_r1～r3（古典，同半徑）、global_thermal、sdp（古典，全局參考）、qaoa_p1～p3（量子） |
| locked | 1＝過程中把所有最佳解都排除了（綁錯） |
| first_lock | 第幾步第一次綁錯 |
| gap | 最後答案比最佳解差多少（0＝找到最佳） |
| signal_seconds | 算訊號花的時間（古典模擬成本） |

主比較：qaoa_p{r} 對 local_r{r}（看得一樣遠）。global_thermal 與 sdp 只當參考上限，random 是下限。
