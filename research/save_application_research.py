"""Save the reviewed research into the existing Obsidian project under its shared lock."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
TAIPEI = timezone(timedelta(hours=8))
import hashlib, json, os, re, shutil

ROOT = Path(__file__).resolve().parents[1]
VAULT = Path(r'<Obsidian 筆記庫>')
PROJECT = VAULT / '10 專案/全國智慧運算與量子資訊創新應用競賽.md'
OVERVIEW = VAULT / '10 專案/專案總覽.md'
REPORT = ROOT / 'output/Recursive_QAOA_應用方向研究_2026-11-12.md'
SHARED = Path.home() / 'OneDrive/Desktop/CODEX/AI共享記憶'
LOCK = SHARED / '.同步鎖'
BASELINE = json.loads((ROOT / 'research/obsidian_research_baseline.json').read_text(encoding='utf-8'))

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

old = {p: p.read_bytes() for p in (PROJECT, OVERVIEW)}
if sha(old[PROJECT]) != BASELINE['sha256']:
    raise RuntimeError('Project changed since research started; re-read before saving.')
report_bytes = REPORT.read_bytes()
report = report_bytes.decode('utf-8')
assert report.startswith('# Recursive QAOA 應用方向研究\n')
heading = '## Recursive QAOA 應用方向研究（2026-11-12）'
project = old[PROJECT].decode('utf-8-sig').replace('\r\n', '\n')
assert heading not in project, 'Research already saved; merge differences instead.'
next_step = '- 11/12新增應用研究（AI建議，尚未定題）：優先評估點雲面片法向量定向、其次感測器極性校正；先做同圖古典強對照與端到端指標，再決定是否實作。完整模型、來源與檢索紀錄見下方「Recursive QAOA 應用方向研究」。本機11/7團隊報告已讀，本輪未重跑。\n'
project = project.replace('## AI 整理的下一步\n\n', '## AI 整理的下一步\n\n' + next_step, 1)
log = '''### 2026-11-12｜Recursive QAOA 應用方向研究（Codex AI交接，11/11開始）

依Johnny本次要求，從RQAOA成對相關、相對符號消去及遞迴縮圖特性查找新應用。公開檢索未找到「點雲面片法向量定向」與「設備通道極性校正」的直接RQAOA研究，暫列優先候選；另列局部二分標籤對齊、基因定相區塊拼接，後者已有QHap近鄰研究。不能保證全球無人研究，這些均為AI提案，Johnny尚未採用或定題。已排除噴漆、無線資源配置、RL、shots控制、回溯／beam等已有先例的泛用新意。

核讀Vault專案原文（與Graphify sync-state SHA256相符）、本機11/7團隊回報與既有文獻查核稿；部分查看加權相關／RQAOA程式入口。11/7報告自述208案例中RQAOA最佳解197例、古典局部搜尋206例、SDP/GW 207例且成本未嚴格匹配，本輪未重跑、未驗收新結果。研究交付包含四方向的變數模型、最小實驗、強對照、失敗條件、實際檢索式及13項可回查來源；詳見本頁下方「Recursive QAOA 應用方向研究」。本機同內容交付：`./output/Recursive_QAOA_應用方向研究_2026-11-12.md`。

下一步建議先做12–20面片的有真值幾何測試，分開求解器與幾何前處理效果；尚未新增應用實驗、真機結果或正式研究結論。保存前依共同鎖、版本檢查與備份規則處理，寫後讀回；一般查詢未重建Graphify，手機／另一端雲端同步未驗證。

'''
project = project.replace('## 工作紀錄\n\n', '## 工作紀錄\n\n' + log, 1)
status = '11/12新增RQAOA應用研究：點雲法向量定向及感測器極性校正列優先候選，未找到直接研究不等於保證首創，尚未定題或執行應用實驗。本次另讀本機11/7報告，其原圖實測與古典比較結果只屬文件回報，未重跑；以下歷史會議／雲端來源狀態保留。\n\n'
project = project.replace('## 目前現況\n\n', '## 目前現況\n\n' + status, 1)
section = heading + '\n\n' + report.split('\n', 1)[1].strip() + '\n\n'
assert project.count('## 資料與來源') == 1
project = project.replace('## 資料與來源\n', section + '## 資料與來源\n', 1)

overview = old[OVERVIEW].decode('utf-8-sig').replace('\r\n', '\n')
overview, count = re.subn(r'^\| \[\[全國智慧運算與量子資訊創新應用競賽\]\] \|.*$', '| [[全國智慧運算與量子資訊創新應用競賽]] | 11/12新增RQAOA應用研究：優先候選為點雲法向量定向／感測器極性校正，尚未定題或實作；11/7本機報告已讀、未重跑，未證量子優勢；原RQAOA／Warm-start分線保留 |', overview, flags=re.M)
assert count == 1
new = {PROJECT: project.encode('utf-8'), OVERVIEW: overview.encode('utf-8')}
stamp = datetime.now(TAIPEI).strftime('%Y%m%d-%H%M%S')
backup = Path.home() / 'Documents/Obsidian Backups' / (stamp + '-rqaoa-applications')
if not SHARED.is_dir():
    raise RuntimeError('Shared folder unavailable; no vault writes.')
LOCK.mkdir()  # Exclusive creation; never remove somebody else's lock.
owner = LOCK / 'codex-rqaoa-owner.json'
temps = []
try:
    owner.write_text(json.dumps({'writer': 'Codex', 'started': datetime.now(TAIPEI).isoformat(), 'purpose': 'RQAOA application research', 'pid': os.getpid()}, ensure_ascii=False), encoding='utf-8')
    for p, raw in old.items():
        if sha(p.read_bytes()) != sha(raw):
            raise RuntimeError(f'Source changed before writing: {p}')
    if REPORT.read_bytes() != report_bytes:
        raise RuntimeError('Research artifact changed before writing.')
    for p in old:
        dest = backup / 'Vault' / p.relative_to(VAULT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dest)
    for p, raw in new.items():
        tmp = p.with_name('.' + p.name + '.codex-rqaoa.tmp')
        if tmp.exists():
            raise RuntimeError('Staging file already exists; do not overwrite.')
        with tmp.open('xb') as f:
            f.write(raw)
        temps.append(tmp)
    for p, tmp in zip(new, temps):
        if sha(p.read_bytes()) != sha(old[p]):
            raise RuntimeError(f'Source changed before replace: {p}')
        os.replace(tmp, p)
        if p.read_bytes() != new[p]:
            raise RuntimeError(f'Write verification failed: {p}')
    saved = PROJECT.read_text(encoding='utf-8')
    assert section in saved and log in saved and next_step in saved
    assert saved.count(heading) == 1
    assert '[S13]' in saved and '## 本次需求原話' in saved
    # All pre-existing lines remain, in order, with only insertions in the project.
    stripped = saved.replace(next_step, '', 1).replace(log, '', 1).replace(status, '', 1).replace(section, '', 1)
    assert stripped == old[PROJECT].decode('utf-8-sig').replace('\r\n', '\n')
    receipt = {'saved_at': datetime.now(TAIPEI).isoformat(), 'project': str(PROJECT), 'overview': str(OVERVIEW), 'report': str(REPORT), 'backup': str(backup), 'source_fresh_at_read': True, 'verified_readback': True, 'original_project_preserved': True, 'graph_rebuilt': False, 'remote_sync_verified': False, 'hashes': {str(p): sha(p.read_bytes()) for p in new}}
    (ROOT / 'research/application_research_save_receipt.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
finally:
    for tmp in temps:
        if tmp.exists():
            tmp.unlink()
    owner.unlink(missing_ok=True)
    LOCK.rmdir()
