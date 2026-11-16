"""Append actual benchmark evidence, retain historical proposals and original words."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,os,shutil,re

ROOT=Path(__file__).resolve().parents[2]
V=Path(r'<Obsidian 筆記庫>')
P=V/'10 專案/全國智慧運算與量子資訊創新應用競賽.md'
O=V/'10 專案/專案總覽.md'
A=ROOT/'output/Recursive_QAOA_應用方向研究_2026-11-12.md'
REPORT=ROOT/'output/點雲法向量_RQAOA與古典方法比較_2026-11-12.md'
RESULTS=ROOT/'research/normal_orientation/results'
TZ=timezone(timedelta(hours=8))
def sha(b):return hashlib.sha256(b).hexdigest()

old={p:p.read_bytes() for p in [P,O,A]}
report=REPORT.read_bytes();summary=(RESULTS/'summary.json').read_bytes();audit=(RESULTS/'audit.json').read_bytes()
receipt=json.loads((ROOT/'research/application_research_save_receipt.json').read_text(encoding='utf-8'))
for p in [P,O]:
    if sha(old[p])!=receipt['hashes'][str(p)]:
        raise RuntimeError(f'Vault changed since reviewed research save; re-read and merge: {p}')
title='### 2026-11-12｜點雲法向量比較已執行：未見RQAOA優勢（Codex AI交接）'
project=old[P].decode('utf-8-sig').replace('\r\n','\n')
assert title not in project
next_old='- 11/12新增應用研究（AI建議，尚未定題）：優先評估點雲面片法向量定向、其次感測器極性校正；先做同圖古典強對照與端到端指標，再決定是否實作。完整模型、來源與檢索紀錄見下方「Recursive QAOA 應用方向研究」。本機11/7團隊報告已讀，本輪未重跑。'
next_new='- 11/12已依Johnny授權執行點雲面片方向比較：36組合成點雲、144個主圖設定及36個幾何補測圖。RQAOA未勝強古典對照；改善來自幾何模型。AI建議先保留負結果與基準，若續做須找真實掃描中的古典困難案例，不先主打量子優勢。詳見本頁「點雲法向量比較實測」，最終選題仍未決。'
assert project.count(next_old)==1
project=project.replace(next_old,next_new,1)
link='[完整實測報告](<./output/點雲法向量_RQAOA與古典方法比較_2026-11-12.md>)'
log=f'''{title}

Johnny要求直接實施比較，本輪已新增並執行`research/normal_orientation/`，舊圖結構研究未覆寫。36組合成點雲（球／圓環／薄板，12、16、20面片、4seeds）產生144個主圖設定，後續在同資料上加36個Dipole交互作用圖；180圖不是180組獨立點雲。108控制圖RQAOA最佳解105例、局部搜尋32起點107、SDP106；以每例RQAOA CPU耗時給古典搜尋則108/108。簡單鄰接模型在圓環／薄板即使窮舉也50%方向錯誤；換Dipole幾何圖後，RQAOA、生成樹、譜、SDP與局部搜尋32起點皆36/36方向全對，未見量子方法獨有收益。

已補144次有限shots決策、參數加密、3個失敗縮圖逐步窮舉及1,656筆原圖重計分（最大差約2.14e-14），核對真值未用於幾何解算。Dipole面片傳播／WNNC dense CPU為有來源與MIT授權的NumPy重實作，未跑官方完整神經網路／treecode環境；公開稀疏寬度l5的WNNC版本在36例也全對。給定局部法線軸、合成資料、理想p=1 CPU解析模擬，不是原始點雲完整重建、官方程式重現或真機優勢。{link}；原始CSV、協定、軌跡與稽核在`./research/normal_orientation/results/`。

本輪實測修正先前「優先研究此應用」的AI推論：能自然建模不代表有比較優勢。下一步建議先找真實掃描及強古典仍難解的條件；沒有替Johnny決定退題或最終題目。已查官方來源版本、備份後依共同鎖更新專案／總覽，讀回核對；Graphify索引因新增內容已過期，本輪未重建，手機同步未驗證。

'''
project=project.replace('## 工作紀錄\n\n','## 工作紀錄\n\n'+log,1)
status='11/12點雲比較已實際執行，取代同日較早「尚未實作」的探索狀態：RQAOA能解但沒有勝過強古典對照；幾何模型改善才是本輪主要收益。完整範圍、數據與限制見下方「點雲法向量比較實測」，尚未決定最終題目。較早研究與會議紀錄保留為歷史。\n\n'
project=project.replace('## 目前現況\n\n','## 目前現況\n\n'+status,1)
section=f'''## 點雲法向量比較實測（2026-11-12）

**實測判斷：沒有證據支持點雲方向校正特別需要RQAOA。** 已知局部法線軸、面片給定、合成小圖下，RQAOA未勝強古典對照。這修正本頁稍早的應用推薦，不抹除原推論與來源。

108個控制圖的圖最佳解命中（不是表面方向全對次數）：

| 方法 | 圖最佳解命中 | CPU求解時間中位數 |
|---|---:|---:|
| 生成樹 | 54/108 | 0.314 ms |
| 譜方法 | 69/108 | 0.077 ms |
| 最強邊古典縮圖 | 95/108 | 0.711 ms |
| SDP，256取整 | 106/108 | 13.791 ms |
| 局部搜尋，32起點 | 107/108 | 0.967 ms |
| RQAOA，p=1精確相關 | 105/108 | 160.252 ms |

控制圖相對符號由模擬真值產生後加0%、10%、25%錯誤，是隔離求解器的oracle測試，不能當作從幾何推得。各方法固定設定，時間不是等量子資源比較；另按每例RQAOA CPU耗時給古典局部搜尋，108/108達圖最佳。所有RQAOA均在本機CPU解析計算，未使用量子API／真機。

36組幾何案例上，簡單局部內積圖讓圓環與薄板各12例有50%方向錯誤，窮舉也不能救；探索性改用Dipole交互作用圖後，RQAOA、生成樹、譜、SDP及32起點局部搜尋全部36例0%錯誤。改善由幾何資訊帶來，不是RQAOA獨有。Dipole面片傳播NumPy版與WNNC dense NumPy公開稀疏l5版同樣36例全對；WNNC l0平均錯誤13.38%，兩設定均保留，不按例挑最好。

有限shots在24個n=12父圖、每種shots各3抽樣種子：256 shots命中65/72，2048 shots命中69/72，參數仍用精確期望值訓練，沒有閘噪音或端到端含shots優化。3個精確相關RQAOA失敗例已定位不可逆損失步驟，擴大gamma範圍與加密網格仍未救回，沒有宣稱參數全域最佳。

**限制：**36個點雲只有144–240點、局部法線軸與分片已給定，薄板只取上下表面；無真實掃描、PCA／神經法線估計、曲面重建與真機驗收。專門幾何方法是依官方來源公式重實作的條件式版本，不是官方完整程式重現。後續Dipole圖與l5是在同批開發資料上的探索性補測，不能作未見資料驗證。

{link} · [主結果CSV](<./research/normal_orientation/results/main.csv>) · [幾何模型補測CSV](<./research/normal_orientation/results/geometry_graph.csv>) · [驗證收據](<./research/normal_orientation/results/audit.json>)。

AI建議：保留本輪負結果與可重現程式；若繼續，要先證明真實資料有古典強對照仍難解的條件。Johnny只授權開始比較，未決定放棄、採用或承諾研究期限。

'''
anchor='## Recursive QAOA 應用方向研究（2026-11-12）\n\n'
assert project.count(anchor)==1
notice='> 11/12後續實測更新：以下是測試前的探索提案；本輪比較尚未看到RQAOA優勢，主要收益來自幾何模型。請先看上方「點雲法向量比較實測」再判斷是否投入。\n\n'
project=project.replace(anchor,section+anchor+notice,1)
overview=old[O].decode('utf-8-sig').replace('\r\n','\n')
overview,count=re.subn(r'^\| \[\[全國智慧運算與量子資訊創新應用競賽\]\] \|.*$', '| [[全國智慧運算與量子資訊創新應用競賽]] | 11/12已跑點雲面片比較：36組合成點雲、180圖設定含探索補測；RQAOA未勝強古典，主要收益來自幾何模型，未證量子優勢。最終題目／真實掃描驗證未定 |',overview,flags=re.M)
assert count==1
proposal=old[A].decode('utf-8-sig').replace('\r\n','\n')
prefix='# Recursive QAOA 應用方向研究\n\n'
assert proposal.startswith(prefix)
proposal=proposal.replace(prefix,prefix+'> 2026-11-12後續實測：這份文件保留測試前的探索推論。已執行的點雲比較未見RQAOA優勢，改善主要來自幾何模型；最新證據見 [實測報告](<./output/點雲法向量_RQAOA與古典方法比較_2026-11-12.md>)。\n\n',1)
new={P:project.encode(),O:overview.encode(),A:proposal.encode()}
shared=Path.home()/'OneDrive/Desktop/CODEX/AI共享記憶';lock=shared/'.同步鎖'
assert shared.is_dir()
lock.mkdir();owner=lock/'codex-normal-benchmark-owner.json';temps=[]
backup=Path.home()/'Documents/Obsidian Backups'/(datetime.now(TZ).strftime('%Y%m%d-%H%M%S')+'-normal-benchmark')
try:
    owner.write_text(json.dumps(dict(writer='Codex',purpose='Executed normal orientation comparison',started=datetime.now(TZ).isoformat(),pid=os.getpid())),encoding='utf-8')
    for p,b in old.items():
        if p.read_bytes()!=b:raise RuntimeError('Concurrent change: '+str(p))
    assert REPORT.read_bytes()==report and (RESULTS/'summary.json').read_bytes()==summary and (RESULTS/'audit.json').read_bytes()==audit
    for p in old:
        dest=backup/('Vault'/p.relative_to(V) if p.is_relative_to(V) else Path('Workspace')/p.relative_to(ROOT))
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    for p,b in new.items():
        tmp=p.with_name('.'+p.name+'.normal-benchmark.tmp')
        with tmp.open('xb') as f:f.write(b)
        temps.append(tmp)
    for p,tmp in zip(new,temps):
        if p.read_bytes()!=old[p]:raise RuntimeError('Concurrent change before replace: '+str(p))
        os.replace(tmp,p);assert p.read_bytes()==new[p]
    saved=P.read_text(encoding='utf-8')
    restored=saved.replace(next_new,next_old,1).replace(log,'',1).replace(status,'',1).replace(section,'',1).replace(notice,'',1)
    assert restored==old[P].decode('utf-8-sig').replace('\r\n','\n')
    assert REPORT.is_file() and (RESULTS/'main.csv').is_file() and saved.count(title)==1
    result=dict(saved_at=datetime.now(TZ).isoformat(),project=str(P),overview=str(O),report=str(REPORT),backup=str(backup),verified_readback=True,prior_project_preserved_except_ai_next_step=True,graph_rebuilt=False,remote_sync_verified=False,hashes={str(p):sha(p.read_bytes()) for p in new})
    (RESULTS/'vault_save_receipt.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result,ensure_ascii=False,indent=2))
finally:
    for p in temps:
        if p.exists():p.unlink()
    owner.unlink(missing_ok=True);lock.rmdir()
