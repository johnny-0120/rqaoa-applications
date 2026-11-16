"""Inspectably summarize the executed benchmark; render static scientific figures."""
import csv,json,statistics,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from benchmark import ROOT,OUT,cloud,error

DEST=ROOT/'output/normal_orientation';DEST.mkdir(parents=True,exist_ok=True)
FONT=Path('C:/Windows/Fonts/msjh.ttc')
if FONT.exists():
    font_manager.fontManager.addfont(str(FONT));plt.rcParams['font.family']=font_manager.FontProperties(fname=str(FONT)).get_name()
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10,'axes.unicode_minus':False,'savefig.facecolor':'white'})
BLUE='#416D91';ORANGE='#C27837';GRAY='#A6ADB4'

def read(name):
    return list(csv.DictReader((OUT/name).open(encoding='utf-8-sig')))

def mean(rows,key):return statistics.mean(float(r[key]) for r in rows)
def median(rows,key):return statistics.median(float(r[key]) for r in rows)

def summary(rows,keys):
    groups={}
    for r in rows:groups.setdefault(tuple(r[k] for k in keys),[]).append(r)
    out=[]
    for key,rr in groups.items():
        out.append(dict(zip(keys,key),cases=len(rr),mean_patch_error=mean(rr,'patch_error'),median_seconds=median(rr,'seconds'),optimal_count=sum(int(r['is_optimal']) for r in rr) if 'is_optimal' in rr[0] else None))
    return out

def run():
    rows=read('main.csv');geo=read('geometry_graph.csv');bud=read('cpu_matched.csv');wn=read('wnnc_sparse_width.csv');shots=read('shots.csv')
    control=[r for r in rows if r['panel']=='control'];geometry=[r for r in rows if r['panel']=='geometry']
    summ=summary(rows,['panel','method']);gsumm=summary(geo,['method']);bsumm=summary(bud,['panel']);ssumm=summary(shots,['shots'])
    result=dict(main=summ,geometry_graph=gsumm,cpu_matched=bsumm,shots=ssumm,geometry_by_shape=summary(geometry,['shape','method']),wnnc_l5_by_shape=summary(wn,['shape']))
    (OUT/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    # Figure 1: controlling graph objective and CPU times, not real orientation.
    names=['mst','spectral','edge_contract','sdp256','local32','rqaoa_p1']
    labels=['生成樹','譜方法','最強邊縮圖','SDP 256取整','局部搜尋32起點','RQAOA p=1']
    rr=[next(x for x in summ if x['panel']=='control' and x['method']==m) for m in names]
    fig,ax=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    colors=[ORANGE if m=='rqaoa_p1' else BLUE for m in names]
    bars=ax[0].barh(labels,[x['optimal_count'] for x in rr],color=colors)
    ax[0].set_xlim(0,118);ax[0].set_xlabel('找到圖模型最佳解的次數／108');ax[0].set_title('控制圖：最佳解命中次數')
    for b,r in zip(bars,rr):ax[0].text(b.get_width()+1,b.get_y()+b.get_height()/2,str(r['optimal_count']),va='center')
    ax[0].invert_yaxis()
    for i,r in enumerate(rr):
        ax[1].scatter(r['median_seconds']*1000,i,color=colors[i],s=60,marker='D' if i==5 else 'o')
        ax[1].annotate(f"{r['median_seconds']*1000:.3g} ms",(r['median_seconds']*1000,i),xytext=(8,5),textcoords='offset points')
    ax[1].set_yticks(range(len(names)),labels);ax[1].invert_yaxis();ax[1].set_xscale('log');ax[1].set_xticks([.1,1,10,100,1000],['0.1','1','10','100','1000']);ax[1].set_xlim(.03,1500);ax[1].set_xlabel('CPU 求解時間中位數，毫秒（對數軸）');ax[1].set_title('同一批控制圖：各方法固定設定')
    fig.suptitle('合成面片控制圖：求解品質與 CPU 成本\n不是等量子資源比較；所有 RQAOA 均為理想 p=1 古典解析模擬',fontsize=12)
    fig.savefig(DEST/'01_solver_comparison.png',dpi=170);plt.close(fig)
    # Figure 2: model specification controls actual orientation.
    fig,ax=plt.subplots(figsize=(8.5,4.7),layout='constrained');x=np.arange(3)
    shapes=['sphere','torus','thin_slab'];labels2=['球面','圓環','薄板上下表面']
    naive=[mean([r for r in geometry if r['shape']==s and r['method']=='rqaoa_p1'],'patch_error')*100 for s in shapes]
    better=[mean([r for r in geo if r['shape']==s and r['method']=='rqaoa_p1'],'patch_error')*100 for s in shapes]
    a=ax.bar(x-.18,naive,.36,label='RQAOA：簡單鄰接圖',color=ORANGE);b=ax.bar(x+.18,better,.36,label='RQAOA：Dipole交互作用圖',color=BLUE)
    for bars,vals in [(a,naive),(b,better)]:
        for bar,val in zip(bars,vals):ax.text(bar.get_x()+bar.get_width()/2,val+1,f'{val:.0f}%',ha='center')
    ax.set_xticks(x,labels2);ax.set_ylim(0,62);ax.set_ylabel('面片方向錯誤率（去除整體翻轉歧義）');ax.set_title('同一批36組點雲：改模型，比換求解器更關鍵');ax.legend(loc='upper left');ax.set_xlabel('Dipole圖上：生成樹、譜方法、SDP、局部搜尋32起點也全部0%錯誤',fontsize=9)
    fig.savefig(DEST/'02_geometry_model.png',dpi=170);plt.close(fig)
    # Figure 3: fixed representative geometry, not selected by best/worst score.
    data=cloud('torus',16,0);main_d=json.loads((OUT/'details.json').read_text(encoding='utf-8'));gd=json.loads((OUT/'geometry_graph_details.json').read_text(encoding='utf-8'))
    d=next(d for d in main_d if d['shape']=='torus' and d['n']==16 and d['seed']==0 and d['panel']=='geometry');g=next(g for g in gd if g['case']==d['case'])
    configs=[('輸入：面片隨機翻轉',np.ones(16)),('RQAOA＋簡單鄰接圖',np.array(d['outputs']['rqaoa_p1']['z'])),('RQAOA＋Dipole圖',np.array(g['outputs']['rqaoa_p1']['z'])),('古典譜方法＋Dipole圖',np.array(g['outputs']['spectral']['z']))]
    fig=plt.figure(figsize=(12,8.6),layout='constrained')
    for k,(title,z) in enumerate(configs,1):
        ax=fig.add_subplot(2,2,k,projection='3d');truth=data['truth']
        if np.mean(z!=truth)>.5:z=-z
        c=data['centers'];v=np.mean(data['normals'],axis=1)*z[:,None];col=np.where(z==truth,BLUE,ORANGE)
        uu,vv=np.meshgrid(np.linspace(0,2*np.pi,36),np.linspace(0,2*np.pi,12));ax.plot_wireframe((1+.35*np.cos(vv))*np.cos(uu),(1+.35*np.cos(vv))*np.sin(uu),.35*np.sin(vv),color=GRAY,alpha=.18,linewidth=.4);ax.scatter(*data['points'].reshape(-1,3).T,s=3,color=GRAY,alpha=.5)
        ax.quiver(c[:,0],c[:,1],c[:,2],v[:,0],v[:,1],v[:,2],length=.35,normalize=True,colors=col,linewidth=1.8)
        ax.set_title(f'{title}\n方向錯誤 {error(z,truth)*100:.0f}%');ax.set_box_aspect((1,1,.45));ax.view_init(24,42);ax.set_axis_off()
    fig.suptitle('固定示例：圓環、16面片、seed=0\n藍箭頭＝方向正確；橘箭頭＝翻反。灰線為生成曲面參考；全局正反僅為評估／顯示而對齊真值。',fontsize=12)
    fig.savefig(DEST/'03_orientation_example.png',dpi=170);plt.close(fig)
    def table(data,methods=None):
        out=['| 方法 | 圖最佳解命中 | 平均方向錯誤率 | CPU時間中位數 |','|---|---:|---:|---:|']
        for r in data:
            if methods and r['method'] not in methods:continue
            out.append(f"| {r['method']} | {r['optimal_count']}/{r['cases']} | {r['mean_patch_error']*100:.2f}% | {r['median_seconds']*1000:.3f} ms |")
        return '\n'.join(out)
    rq=next(r for r in summ if r.get('panel')=='control' and r['method']=='rqaoa_p1');ls=next(r for r in summ if r.get('panel')=='control' and r['method']=='local32')
    report=f'''# RQAOA 與古典方法：點雲面片方向校正實測

日期：2026-11-12，Asia/Taipei。作者：Codex AI。Johnny已授權開始測試，尚未表示採用本題作最終參賽題目。

## 結論

**本輪實測不支持「這個應用特別需要RQAOA，或RQAOA優於古典方法」。** RQAOA可以解這種正反向模型，但本輪最明顯的改善來自幾何模型。納入Dipole交互作用後，RQAOA、生成樹、譜方法、SDP與32起點局部搜尋都在36組合成點雲上得到0%方向錯誤，沒有出現RQAOA獨有的收益。

108個已知相對符號加噪的控制圖中，RQAOA找到圖最佳解105次，32起點局部搜尋107次、SDP 106次；以每例RQAOA的CPU耗時給古典局部搜尋，後者108次全找到最佳。這是小型、特定資料與指定實作的比較，不能外推所有點雲或一般複雜度。

因此修正先前的推薦：本題可作混合求解與幾何建模的學習／負結果研究；目前不足以作「量子方法突破」的競賽主張。若繼續，應先證明真實掃描中存在強古典方法仍難解的具體問題，而非先增加量子模組。

## 實際做了什麼

- 36組可重建合成點雲：球面、圓環、薄板上下表面 × 12/16/20面片 × 4個seeds。每片12點，144–240點／組。
- 已供應局部一致的法線軸，只隨機翻轉每片；點位置高斯擾動標準差0.003、法線向量擾動標準差0.03後正規化。**沒有從原始掃描做PCA／神經網路法線估計，局部面片分割也是給定的。**
- 主比較144個圖設定：108個控制圖、36個簡單幾何圖。後續另加36個Dipole交互作用圖，共180個圖設定，仍來自同36組點雲，不能當180組獨立資料。
- 所有圖求解器共用相同K矩陣、原始計分及窮舉最佳值；圖目標為最大化 Σi<j Kij zi zj。方向真值錯誤率允許全體翻轉，逐案例等權平均，不是跨所有點池化。
- 另有144次有限量測遞迴、144個相近CPU時間的古典搜尋、WNNC兩組公開寬度設定、參數搜尋加密與失敗步驟稽核。這些都是補測，不是新增獨立點雲。

## 1. 已知相對關係加噪：比較求解器

控制圖使用對稱4近鄰聯集，必要時以最近跨分量邊連通。邊權由距離衰減與法線軸內積絕對值決定；**相對正負號由模擬真值產生，再以0%、10%、25%機率翻轉**，是刻意隔離求解器的oracle控制，不是假裝從幾何推得。

同一幾何在三個噪音條件重複使用，同seed亦耦合噪音，沒有把108圖當獨立樣本作顯著性檢定。

{table([r for r in summ if r['panel']=='control'])}

![控制圖求解品質與CPU成本]({DEST.as_posix()}/01_solver_comparison.png)

方向錯誤率與圖最佳解命中是不同問題：當輸入的相對證據有錯誤，圖的最佳解也可能偏離真值。因此某方法方向錯誤略低，不代表它比較會解圖或有較強最佳化能力。不能用105/108直接解讀為105組表面方向全部正確。

RQAOA時間中位數 {rq['median_seconds']*1000:.3f} ms；32起點局部搜尋 {ls['median_seconds']*1000:.3f} ms，兩者中位數比約 {rq['median_seconds']/ls['median_seconds']:.0f} 倍。這是本機小圖CPU實作成本，**不是量子硬體速度比較，也不是相同總量子資源比較**。計時不含共用建圖與窮舉驗收，SDP含求解與256次取整。沒有使用量子API或真機。

每例另以主比較RQAOA的CPU耗時作局部搜尋時限，允許最後一次重啟完成並記實際耗時：控制圖108/108、簡單幾何圖36/36到達圖最佳值。搜尋重啟數、時限與實際耗時可查 `cpu_matched.csv`。這只對齊此CPU程式的耗時，不聲稱對齊真機成本。

## 2. 直接用局部法線相似度建圖：模型失配

簡單幾何圖不使用真值決定符號，只根據面片平均法線內積的正負，以及中心距離建圖。RQAOA與窮舉均在球面12/12組0%方向錯誤、圓環12/12與薄板12/12組50%錯誤；所有36圖仍都達到該圖最佳解。

這說明正反向變數建模自然，卻不代表局部內積就是正確的幾何目標。薄板兩側法線本應相反，靠近也不代表朝同邊；本輪薄板只取上下表面，沒有側壁，是特別容易揭露這個假設的受控案例。

## 3. 加入Dipole幾何交互作用：所有主要方法都能解對

因主比較揭露模型失配，本輪在相同開發資料上增加探索性模型：Kij取第j片形成的dipole場與第i片法線的總交互作用，對稱化、刪除對角線，再按最大絕對權重縮放。這是我們用來比較的圖目標，**不是宣稱Dipole原論文採這個全局Ising最佳化目標，也不宣稱新模型首創**。建圖未使用方向真值。

這是同一批資料上的探索性補測，尚非未見資料驗證。

{table(gsumm)}

![幾何模型影響]({DEST.as_posix()}/02_geometry_model.png)

![固定圓環方向示例]({DEST.as_posix()}/03_orientation_example.png)

譜方法與生成樹等簡單方法已全對，這批例子無法展現需要較複雜求解器的困難性。不能用RQAOA的36/36作量子優勢證據。

## 4. 專門幾何方法的比較範圍

官方程式皆已查閱並保存來源：WNNC commit `175b116580c9c49e8f8d421fb17eee10fcddc198`；Dipole commit `72ba66b6a54858ebdb34de44afaefaeabfeee81b`。本機研究venv缺torch／torch_geometric／Open3D；WNNC原setup同時建CPU/CUDA擴充，沒有直接執行官方完整環境。本輪以NumPy執行相關公式，保留MIT授權，不能稱為官方程式重現。

| 幾何方法 | 實際執行範圍 | 36組平均面片方向錯誤率 |
|---|---|---:|
| Dipole patch NumPy | 依官方field與面片傳播流程重實作，diffuse=False；給定局部面片／法線，省略原學習步驟 | {mean([r for r in geometry if r['method']=='dipole_patch_numpy'],'patch_error')*100:.2f}% |
| WNNC dense NumPy l0 | 官方A／AT／G逐點kernel與40次更新的float64密集計算；公開乾淨點寬度0.002–0.016 | {mean([r for r in geometry if r['method']=='wnnc_dense_numpy_l0'],'patch_error')*100:.2f}% |
| WNNC dense NumPy l5 | 同公式，另用官方針對稀疏點的公開寬度0.05–0.2，40次更新 | {mean(wn,'patch_error')*100:.2f}% |

WNNC從點位置估計新法線，再以每片與輸入法線軸的總內積投影成正反符號，供條件式面片比較；沒有驗收其連續法線精度或重建曲面品質。使用密集精確逐對kernel取代treecode近似與加速，未與官方binary逐輸出對照。l0/l5兩組全部回報，不按個別例子挑最好參數；l5是看見l0不適合稀疏輸入後，依官方文件補上的探索性設定。

Dipole也未執行原神經網路分片與方向估計，因此這些結果只支持本條件式階段，不等同原論文端到端效能。各類方法使用的幾何資訊不同，不能把它們和「只看K矩陣」的求解器混成單一公平排行榜。

## 5. RQAOA有限量測、參數與失敗稽核

RQAOA固定p=1，每次選全部點對中絕對相關最大者，字典序處理1e-9內平手；一次消去一點，到4點窮舉。每輪按最大耦合重縮放，gamma先掃0到2π的129點，最多細化4峰，beta在每個gamma解析求最佳。實數權重沒有假定2π為全域週期，這只是有界啟發式搜尋，沒有參數全域最佳證書。

有限量測選定所有n=12、控制噪音10%與25%的24個父圖，各跑3個抽樣種子；共同Z基底bitstrings估計所有相關。**參數訓練仍用理想精確期望值，只在消去決策抽樣，不含閘噪音／讀出錯誤／排隊／端到端含shots訓練。**

| 每階段決策shots | 找到圖最佳解 | 平均方向錯誤率 |
|---|---:|---:|
| 256 | {ssumm[0]['optimal_count']}/{ssumm[0]['cases']} | {ssumm[0]['mean_patch_error']*100:.2f}% |
| 2048 | {ssumm[1]['optimal_count']}/{ssumm[1]['cases']} | {ssumm[1]['mean_patch_error']*100:.2f}% |

72次不是72個獨立圖，是24圖各3次抽樣。詳細總shots逐次可查 `shots.csv`。

所有n=16、seed=0、控制噪音10%／25%的6例，改257網格與gamma上界4π，最終分數均未變；另對主比較3個失敗例改513網格與4π，仍未救回原圖最佳。這些補測不能證明已找到全域最佳參數，只排除這兩種簡單擴搜就翻轉結論的情況。

3個失敗例均逐步窮舉受限問題，驗證選定關係何時不可逆排除最佳解。例case11：原圖最佳3.007819，RQAOA最終2.881404；各步最佳值與損失見 `failure_audit.json`。這些分數是各圖自身單位，不用跨圖相除宣稱近似比。

## 驗證與限制

- 相關公式對Qiskit statevector最大誤差約2e-15；期望值與縮圖含常數代入已核對。WNNC的A／AT伴隨性及G的有限差分符號核對通過。
- 1,656筆儲存輸出以原始邊清單重新計分，最大差約2.14e-14。幾何建圖移除真值資料後不變；隨機翻轉輸入法線時，Dipole圖及兩種幾何方法符合符號轉換，沒有在解算中使用真值。
- SDP為SCS數值解；最小特徵值約-5.53e-7，relaxation低於枚舉最佳的最大數值偏差約8.71e-7，不把它當嚴格最佳性證書。所有最佳值判定用枚舉原圖，不用SDP值。
- 合成點雲、少量seeds、小面片數、已知局部法線軸、缺失側壁及稀疏取樣都限制外推。沒有真實掃描資料、曲面重建、真機或p>1測試。
- 初始小型試跑與後續模型修訂均留檔；本輪屬探索性實驗，沒有宣稱正式預先註冊或獨立測試集。

## 建議下一步

先保留本輪負結果與可重現基準。若仍選這個方向，下一個門檻是公開真實掃描的端到端比較：使用實際法線估計與面片化，補視點／可見性，跑官方幾何方法，測方向與重建品質；需先找到古典強對照仍困難的案例。現在不建議只因RQAOA能解而把它放進產品或主打量子優勢。

## 檔案與重跑

工作目錄：`./`。沿用既有`.venv`及`requirements-lock.txt`，沒有新增套件或改動舊圖研究結果。

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONIOENCODING='utf-8'
.\\.venv\\Scripts\\python.exe research/normal_orientation/benchmark.py validate
.\\.venv\\Scripts\\python.exe research/normal_orientation/benchmark.py run
.\\.venv\\Scripts\\python.exe research/normal_orientation/followup.py
.\\.venv\\Scripts\\python.exe research/normal_orientation/geometry_graph.py
.\\.venv\\Scripts\\python.exe research/normal_orientation/audit.py
.\\.venv\\Scripts\\python.exe research/normal_orientation/report.py
```

結果在 `research/normal_orientation/results/`；`main.csv`、`details.json`是主比較；`geometry_graph.csv`、`geometry_graph_details.json`是探索性幾何模型補測；`summary.json`是聚合結果；`validation.json`與`audit.json`是驗證收據。重跑會更新本輪同名輸出，舊`results/`及`research/q9_results/`不受影響。

## 研究來源

- [RQAOA原始論文與附錄C](https://arxiv.org/html/1910.08980v1)：p=1相關、消去及古典可模擬性。此CPU解析實作不等於量子加速。
- [Dipole官方程式](https://github.com/galmetzer/dipole-normal-prop/tree/72ba66b6a54858ebdb34de44afaefaeabfeee81b)：field_utils.py、面片傳播；授權在vendor_inspection保存。
- [WNNC官方程式](https://github.com/jsnln/WNNC/tree/175b116580c9c49e8f8d421fb17eee10fcddc198)：main_wnnc.py、CPU kernels、l0/l5公開設定；授權在vendor_inspection保存。
- [WNNC論文](https://arxiv.org/abs/2405.16634)與[Dipole論文](https://arxiv.org/abs/2105.01604)：作方法定位，沒有把本輪重實作當論文結果驗收。

來源與狀態：本次使用者要求「請你開始來做測試，直接開始實施這幾個的比較。」本輪已實際執行上述本機測試；未定最終題目，也沒有將AI建議寫成Johnny決策。
'''
    path=ROOT/'output/點雲法向量_RQAOA與古典方法比較_2026-11-12.md';path.write_text(report,encoding='utf-8')
    receipt=dict(report=str(path),figures=[str(p) for p in DEST.glob('*.png')],summary=str(OUT/'summary.json'),hash=hashlib.sha256(path.read_bytes()).hexdigest())
    (OUT/'report_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(receipt,ensure_ascii=False))

if __name__=='__main__':run()
