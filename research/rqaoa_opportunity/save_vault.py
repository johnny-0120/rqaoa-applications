"""Version-checked Obsidian research log; no graph rebuild or extra scheduling."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,os,re,shutil

ROOT=Path(__file__).resolve().parents[2];V=Path(r'<Obsidian 筆記庫>');P=V/'10 專案/全國智慧運算與量子資訊創新應用競賽.md';O=V/'10 專案/專案總覽.md';OUT=Path(__file__).resolve().parent/'results';REPORT=ROOT/'output/RQAOA_其他優勢方向搜尋與實測_2026-11-12.md';TZ=timezone(timedelta(hours=8))

def main():
    old={p:p.read_bytes() for p in [P,O]};prior=json.loads((ROOT/'research/normal_orientation/results/vault_save_receipt.json').read_text(encoding='utf-8'))
    for p in old:
        if hashlib.sha256(old[p]).hexdigest()!=prior['hashes'][str(p)]:raise RuntimeError('Vault changed since last reviewed save; re-read and merge: '+str(p))
    a=json.loads((OUT/'analysis.json').read_text());assert len(a['dense128'])==4
    dependencies={p:p.read_bytes() for p in [REPORT,OUT/'analysis.json',OUT/'fast_confirmation.csv',OUT/'dense128.csv']}
    project=old[P].decode('utf-8-sig').replace('\r\n','\n');overview=old[O].decode('utf-8-sig').replace('\r\n','\n')
    title='### 2026-11-12｜其他RQAOA優勢搜尋與實測：噴漆有局部收益，未勝強古典（Codex AI交接）'
    assert title not in project
    nxt='- 11/12其他RQAOA方向已實測：四類32例篩選、兩類40例新種子確認、噴漆40例全新確認與128變數4例探測。噴漆參數轉移對作者RSG少約12–14%換色，但同時間禁忌／SDP仍更好。AI建議只保留噴漆作有成本邊界的候選，不主張強古典或量子優勢；最終題目未定。'
    lines=project.splitlines();idx=lines.index('## AI 整理的下一步')+2;assert lines[idx].startswith('- 11/12已依Johnny授權');oldnext=lines[idx];lines[idx]=nxt;project='\n'.join(lines)+'\n' if project.endswith('\n') else '\n'.join(lines)
    reportlink='[其他優勢搜尋與實測報告](<./output/RQAOA_其他優勢方向搜尋與實測_2026-11-12.md>)'
    log=f'''{title}

Johnny要求找其他較有優勢的方向並直接測試。本輪查原始論文與作者程式，執行噴漆BPSP、稠密正權重圖、SK高斯自旋玻璃、三正則正負耦合圖的32例小圖篩選；前兩類各20例新種子確認。參數轉移只從初篩資料學得中位數規則，凍結後測40個全新噴漆例（32／64變數各20），另補128變數稠密圖4例。共116個不同問題實例，pilot與同例不同方法／隨機重複不另算獨立資料；52個小例有窮舉參考，其餘僅最佳已知可行值。

最終新例：32變數RQAOA轉移平均換色10.90、作者RSG12.60、同CPU禁忌與SDP＋局部10.65；64變數分別21.55、24.40、20.30。RQ對RSG勝／平／負分別15／5／0、17／3／0，平均少約13.5%、11.7%，但RQ花更多CPU。對同CPU禁忌為0／16／4、0／9／11，原少量案例相對同時間局部搜尋的好處未穩定保留。參數轉移降低重調參成本的收益與「勝強古典」分開；128四例也未勝SDP／禁忌，非作者原研究完整重現或普遍否證。

程式／協定／圖矩陣／解／縮圖軌跡在`./research/rqaoa_opportunity/`；{reportlink}。作者RSG、RG固定版本已核讀並執行，6車出版例greedy6、RG5、RSG4核對通過；初版插入啟發式誤名RG已更名，不改數值，真正RG／RSG另列。{a['audit']['rescored_outputs']}筆解已以原圖重計分，最大差{a['audit']['max_rescore_error']:.3g}。所有RQ為p=1精確相關CPU模擬，無真機、無完整含shots訓練；沒有新演算法或新應用首創主張。

AI判斷：四類中噴漆／固定順序二元切換最值得保留，但只支持局部品質收益與成本研究，尚不滿足「勝過強古典」條件。Johnny尚未決定最終題目或採用本建議。按共同鎖、備份、版本核對及讀回保存；不改舊點雲結果，不重建全庫Graphify。

'''
    project=project.replace('## 工作紀錄\n\n','## 工作紀錄\n\n'+log,1)
    status='11/12新一輪其他方向比較已完成：噴漆參數轉移對作者RSG有約12–14%換色品質收益，但未勝同CPU禁忌／SDP；SK、三正則及稠密圖亦未形成穩定強古典優勢。116個問題實例的範圍與限制見下方新工作紀錄／報告；沒有量子優勢、原創或最終定題結論。先前點雲負結果與較早探索原樣保留。\n\n'
    project=project.replace('## 目前現況\n\n','## 目前現況\n\n'+status,1)
    section=f'''## 其他RQAOA方向實測（2026-11-12）

**最有機會的是二元噴漆換色，但目前優勢僅相對較簡單策略；強古典優勢未成立。** 這是文獻已有的用途，不是新題目首創。

全新40例、固定轉移參數：

| 方法 | 32變數20例，平均換色 | 64變數20例，平均換色 |
|---|---:|---:|
| 作者RSG專用貪婪 | 12.60 | 24.40 |
| 古典局部搜尋32起點 | 11.60 | 25.20 |
| RQAOA參數轉移 | 10.90 | 21.55 |
| 同CPU局部搜尋 | 10.87 | 21.48 |
| 同CPU禁忌搜尋 | 10.65 | 20.30 |
| SDP＋局部搜尋 | 10.65 | 20.30 |

每例同時間古典方法3種子先在例內平均，不當120個獨立例。RQ對RSG平均少1.70／2.85次換色，但RSG快很多；對同CPU禁忌沒有一例平均更好。最終資料未重調轉移規則。參數轉移可省反覆調參，但p=1相關本就可有效古典計算，不代表量子硬體加速。

其餘結構：32例小圖四類初篩後，噴漆與稠密正權重圖進入40例確認；稠密128變數另以接近文獻的分布與近零第一峰探測4例，未勝SDP或禁忌。這不是作者完整資料／配置重現，不據此否定所有稠密圖或更高p。SK／三正則只做小圖，暫不優先。

{reportlink} · [最終40例CSV](<./research/rqaoa_opportunity/results/fast_confirmation.csv>) · [聚合與驗證](<./research/rqaoa_opportunity/results/analysis.json>)。

若競賽必須勝過強古典，這輪還沒有可據以定題的方向；若接受方法適用性與負結果研究，噴漆有較清楚的應用和可展示收益。這是AI建議，不自動替Johnny採用、放棄或排期限。

'''
    anchor='## 點雲法向量比較實測（2026-11-12）\n';assert project.count(anchor)==1;project=project.replace(anchor,section+anchor,1)
    overview,count=re.subn(r'^\| \[\[全國智慧運算與量子資訊創新應用競賽\]\] \|.*$', '| [[全國智慧運算與量子資訊創新應用競賽]] | 11/12其他RQAOA方向116例實測：噴漆轉移比作者RSG少約12–14%換色，但未勝同CPU禁忌／SDP；點雲負結果保留。無量子優勢／原創／最終定題結論 |',overview,flags=re.M);assert count==1
    new={P:project.encode(),O:overview.encode()};shared=Path.home()/'OneDrive/Desktop/CODEX/AI共享記憶';lock=shared/'.同步鎖';assert shared.is_dir();lock.mkdir();owner=lock/'codex-rqaoa-opportunity-owner.json';temps=[];backup=Path.home()/'Documents/Obsidian Backups'/(datetime.now(TZ).strftime('%Y%m%d-%H%M%S')+'-rqaoa-opportunity')
    try:
        owner.write_text(json.dumps(dict(writer='Codex',purpose='RQAOA opportunity benchmark handoff',started=datetime.now(TZ).isoformat(),pid=os.getpid())),encoding='utf-8')
        for p,b in {**old,**dependencies}.items():
            if p.read_bytes()!=b:raise RuntimeError('Concurrent change: '+str(p))
        for p in old:
            dest=backup/'Vault'/p.relative_to(V);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
        for p,b in new.items():
            tmp=p.with_name('.'+p.name+'.rqaoa-opportunity.tmp')
            with tmp.open('xb') as f:f.write(b)
            temps.append(tmp)
        for p,tmp in zip(new,temps):
            if p.read_bytes()!=old[p]:raise RuntimeError('Concurrent change before replace: '+str(p))
            os.replace(tmp,p);assert p.read_bytes()==new[p]
        saved=P.read_text(encoding='utf-8');restored=saved.replace(nxt,oldnext,1).replace(log,'',1).replace(status,'',1).replace(section,'',1)
        assert restored==old[P].decode('utf-8-sig').replace('\r\n','\n')
        rec=dict(saved_at=datetime.now(TZ).isoformat(),project=str(P),overview=str(O),report=str(REPORT),backup=str(backup),readback_verified=True,historical_project_preserved_except_AI_next_step=True,graph_rebuilt=False,remote_sync_verified=False,hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in new})
        (OUT/'vault_save_receipt.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(rec,ensure_ascii=False,indent=2))
    finally:
        for p in temps:
            if p.exists():p.unlink()
        owner.unlink(missing_ok=True);lock.rmdir()

if __name__=='__main__':main()
