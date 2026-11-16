"""Create a source-backed Markdown research handoff and static scientific plots."""
from pathlib import Path
import csv,json,statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from sweep import ROOT,OUT

FIG=ROOT/'output/rqaoa_opportunity';FIG.mkdir(parents=True,exist_ok=True)
font=Path('C:/Windows/Fonts/msjh.ttc')
if font.exists():font_manager.fontManager.addfont(str(font));plt.rcParams['font.family']=font_manager.FontProperties(fname=str(font)).get_name()
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'font.size':10,'savefig.facecolor':'white'})
BLUE='#426F91';ORANGE='#C27B3A';GRAY='#C0C6CC'
def read(name):return list(csv.DictReader((OUT/name).open(encoding='utf-8-sig')))
def mean(a):return float(statistics.mean(a))

def run():
    a=json.loads((OUT/'analysis.json').read_text(encoding='utf-8'));summ=a['fresh_paint_summary'];ds=json.loads((OUT/'discovery_details.json').read_text());cs=json.loads((OUT/'confirmation_details.json').read_text());main=read('confirmation.csv');transfer=read('transfer.csv');dense=read('dense128.csv')
    methods=['paint_RSG_author','edge_contract','local32','rqaoa_transfer','local_cpu_matched','tabu_cpu_matched','sdp1024_local']
    labels=['作者RSG','最強邊縮圖','局部32起點','RQAOA參數轉移','同時間局部搜尋','同時間禁忌搜尋','SDP＋局部搜尋']
    fig,axes=plt.subplots(1,2,figsize=(13.5,5.6),layout='constrained')
    for ax,n in zip(axes,[32,64]):
        rr=[next(r for r in summ if r['n']==n and r['method']==m) for m in methods]
        vals=[r['mean_swaps'] for r in rr];bars=ax.barh(labels,vals,color=[ORANGE if m=='rqaoa_transfer' else BLUE for m in methods]);ax.invert_yaxis();ax.set_xlim(0,max(vals)*1.18);ax.set_xlabel('平均換色次數（越少越好）');ax.set_title(f'{n}變數：20個全新案例')
        for bar,v in zip(bars,vals):ax.text(v+.12,bar.get_y()+bar.get_height()/2,f'{v:.2f}',va='center')
    fig.suptitle('噴漆換色：RQAOA優於簡單策略，但未勝強古典搜尋\n參數只用舊篩選資料學得；局部／禁忌的隨機重複先在每圖內平均。',fontsize=12)
    fig.savefig(FIG/'01_fresh_paint_quality.png',dpi=165);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.6),layout='constrained')
    for ax,n in zip(axes,[32,64]):
        ms=['paint_RSG_author','rqaoa_transfer','tabu_cpu_matched','sdp1024_local'];rr=[next(r for r in summ if r['n']==n and r['method']==m) for m in ms]
        for r,m in zip(rr,ms):
            ax.scatter(r['median_seconds']*1000,r['mean_swaps'],s=65,color=ORANGE if m=='rqaoa_transfer' else BLUE,marker='D' if m=='rqaoa_transfer' else 'o')
            label={'paint_RSG_author':'RSG','rqaoa_transfer':'RQ轉移','tabu_cpu_matched':'禁忌','sdp1024_local':'SDP＋局部'}[m]
            ax.annotate(label,(r['median_seconds']*1000,r['mean_swaps']),xytext=(5,8),textcoords='offset points')
        ax.set_xscale('log');ax.set_xticks([.1,1,10,100,1000],['0.1','1','10','100','1000']);ax.set_xlim(.1,1500);ax.set_xlabel('CPU時間中位數，毫秒（對數軸）');ax.set_ylabel('平均換色次數（越少越好）');ax.set_title(f'{n}變數，同一批20例');ax.grid(axis='y',alpha=.15)
    fig.suptitle('噴漆換色：品質與CPU成本\n所有RQAOA為p=1解析CPU模擬；轉移訓練成本另列，不是量子硬體加速。',fontsize=12)
    fig.savefig(FIG/'02_quality_cost.png',dpi=165);plt.close(fig)
    # Small-case family table avoids averaging incompatible application scores.
    family_names={'paintshop':'噴漆換色','dense_positive':'稠密正權重切割','sk_gaussian':'SK高斯自旋玻璃','signed_regular3':'三正則正負耦合圖'}
    screen=['| 候選 | 案例數 | RQAOA找到確定最佳解 | 對SDP1024：勝／平／負 |','|---|---:|---:|---:|']
    for f in ['paintshop','dense_positive','sk_gaussian','signed_regular3']:
        rr=[d for d in ds if d['family']==f];diff=[d['outputs']['rqaoa_p1']['score']-d['outputs']['sdp1024']['score'] for d in rr]
        screen.append(f"| {family_names[f]} | {len(rr)} | {sum(abs(d['outputs']['rqaoa_p1']['score']-d['exact_optimum'])<1e-7 for d in rr)} | {sum(v>1e-7 for v in diff)}／{sum(abs(v)<=1e-7 for v in diff)}／{sum(v< -1e-7 for v in diff)} |")
    fresh=['| 方法 | 32變數平均換色 | 64變數平均換色 |','|---|---:|---:|']
    for m,label in zip(methods,labels):
        x=[next(r['mean_swaps'] for r in summ if r['n']==n and r['method']==m) for n in [32,64]];fresh.append(f'| {label} | {x[0]:.2f} | {x[1]:.2f} |')
    newcomparisons=a['fresh_paint_comparisons'];rsg32=next(r for r in newcomparisons if r['n']==32 and r['comparator']=='paint_RSG_author');rsg64=next(r for r in newcomparisons if r['n']==64 and r['comparator']=='paint_RSG_author')
    full64=[r for r in main if r['family']=='paintshop' and int(r['n'])==64 and r['method']=='rqaoa_p1'];fast64=[r for r in transfer if r['family']=='paintshop' and int(r['n'])==64 and r['method']=='rqaoa_transfer']
    fulltime=statistics.median(float(r['seconds']) for r in full64);fasttime=statistics.median(float(r['seconds']) for r in fast64)
    training=sum(d['outputs']['rqaoa_p1']['info']['evaluations'] for d in ds if d['family']=='paintshop')
    largetable=['| 密度 | seed | RQAOA切割值 | SDP1024 | SDP＋局部 | 同CPU禁忌 |','|---|---:|---:|---:|---:|---:|']
    for p,seed in [(.3,400),(.3,401),(.7,400),(.7,401)]:
        rr=[r for r in dense if float(r['density'])==p and int(r['seed'])==seed]
        if not rr:continue
        vals=[next(float(r['maxcut']) for r in rr if r['method']==m) for m in ['rqaoa_first_peak','sdp1024','sdp1024_local','tabu_cpu_matched']]
        largetable.append(f'| {p} | {seed} | '+' | '.join(f'{v:.0f}' for v in vals)+' |')
    report='\n'.join([
    '# RQAOA 還有哪些較有機會的方向？搜尋與實測',
    '',
    '日期：2026-11-12，Asia/Taipei。作者：Codex AI。本輪依Johnny要求直接找候選、實作與比較；不是已採用的最終題目，也不宣稱創新或量子優勢。',
    '',
    '## 判斷',
    '',
    '**四類中，二元噴漆換色最值得保留；但目前可確認的是「優於簡單專用貪婪／少量局部搜尋」，不是「優於強古典方法」。** 參數轉移能大幅省掉RQAOA反覆調參成本；在最終40個全新案例上，它仍優於作者RSG，但同時間禁忌搜尋與SDP更好。',
    '',
    '這個結果比點雲更有結構上的理由：噴漆本來就有成對必須異色的約束，可用一個變數代表同車型的兩次出現，沒有另外加入大量硬限制。相鄰車的換色目標會形成稀疏正負耦合；RQAOA能逐步固定關係。但這是已有人研究的應用，不能當作未被做過的新用途。',
    '',
    '## 文獻如何影響選擇',
    '',
    '- [噴漆RQAOA研究](https://arxiv.org/abs/2507.10908)將參數轉移和遞迴縮圖結合，適合測試「是否能省掉調參又保留品質」。本輪是自行實作的簡單轉移規則，沒有冒充作者方法。',
    '- [較新的噴漆比較](https://arxiv.org/abs/2509.15294)指出RQAOA的規模收益可能減弱；因此納入作者原始RSG程式與較大的新案例，避免只看基本greedy。',
    '- [近最優p=1調參研究](https://arxiv.org/abs/2501.16419)報告128/256變數加權圖的表現。除小圖篩選外，本輪也做128變數、近零第一峰、edge-only、cutoff8的四例探測；不是相同資料／軟體的完整重現。',
    '- [原始RQAOA](https://arxiv.org/abs/1910.08980)及[beam RQAOA](https://arxiv.org/abs/2611.08646)顯示挫折性關係值得研究，也指出不可逆消去風險。本輪不新增beam或RL來掩蓋基本方法失敗。',
    '',
    '## 1. 四類初篩：32個16／20變數案例',
    '',
    '\n'.join(screen),
    '',
    '每類16／20變數各4種子。所有小圖用窮舉最佳值核對；SK高斯圖與三正則正負圖是結構測試，不是假裝已驗證某個商業應用。稠密圖按文獻分布生成：ER p=0.7、Gaussian均值50、變異數30後取整；不含外場。',
    '',
    '標準QAOA以精確分布計算「256次量測最好答案的期望」，另算256×(n−4)次的最好答案期望，避免拿QAOA平均值和RQAOA單一結果相比。初篩依「RQAOA相對後者的平均正規化差」取前兩類：噴漆、稠密正權重圖。第二名的平均差仍為負，入選只是多測，不代表已顯示優勢。只對齊決策抽樣數，不等於完整訓練與量子資源公平。',
    '',
    '## 2. 新種子確認與參數轉移',
    '',
    '兩類各用20個新案例：20變數10例、32變數6例、64變數4例。RQAOA使用p=1精確相關、每輪重新找參數、cutoff4；完整CPU代價與所有殘餘圖軌跡保存。大於20變數只稱最佳已知可行值，不假定全域最佳。',
    '',
    f'這批噴漆64變數的每例重調參時間中位數 {fulltime:.3f} 秒；轉移版 {fasttime:.3f} 秒，約省 {fulltime/fasttime:.1f} 倍單次推論時間。兩者平均換色均20次。這是CPU演算法內部比較，不是量子硬體加速；轉移的前置訓練成本不在單次時間內。',
    '',
    f'轉移參數只用初篩噴漆8例、112個殘餘圖，完整RQ調參共 {training} 次期望函數評估；凍結每類gamma×RMS加權度、beta的中位數。新例每輪以alpha／RMS度取得gamma，不再搜尋。沒有用確認資料調參、挑最好的參數組或每例選版本。',
    '',
    '最初64變數4例的轉移版曾以平均換色20次優於同時間局部搜尋20.75次；不能只留下這個好看的結果。故新增下面40例，檢查能否保留。',
    '',
    '## 3. 最終40個全新噴漆案例：局部收益成立，強古典優勢未成立',
    '',
    '32變數seed100–119、64變數seed200–219，各20例。參數規則完全凍結；RQ耗時取3次中位數，據此給局部搜尋／禁忌搜尋每例相同CPU預算，各3個隨機種子先在例內平均。沒有把120次隨機重複當成120個獨立案例。',
    '',
    '\n'.join(fresh),
    '',
    f'RQAOA對作者RSG：32變數勝／平／負 = {rsg32["wins"]}／{rsg32["ties"]}／{rsg32["losses"]}，平均少1.70次（約13.5%）；64變數 = {rsg64["wins"]}／{rsg64["ties"]}／{rsg64["losses"]}，平均少2.85次（約11.7%）。這是品質改善，RQAOA花的CPU時間更多，不是對RSG的速度優勢。',
    '',
    '但在新例中，RQAOA對同時間局部搜尋的平均差是32變數多0.033次、64變數多0.067次；先前少量案例的同時間優勢沒有穩定保留。對同時間禁忌搜尋，32變數0勝／16平／4負，64變數0勝／9平／11負。這使「能勝強古典」的說法不成立。',
    '',
    f'![新噴漆案例品質比較]({FIG.as_posix()}/01_fresh_paint_quality.png)',
    '',
    f'![品質與CPU成本]({FIG.as_posix()}/02_quality_cost.png)',
    '',
    '作者RSG／RG程式來自[固定版本](https://github.com/vijeycreative/XQAOA-BPSP/tree/c380c935f48d328cba85d259f244cd4866100404)，已完整查閱相關演算法並直接執行，不是自行放一個弱版本。出版6車例核對：greedy=6次、RG=5次、RSG=4次；自寫canonical RG另與作者程式對照通過。初版另有「按首次出現插入」啟發式被誤名為RG，正式分析前已改名，保留數值及修正紀錄；最終比較使用真正的作者RSG。',
    '',
    '禁忌搜尋為自寫古典方法：任期5–9步、允許創新最佳的aspiration、10n步未改善則重啟；增益公式已逐步重計分核對。未宣稱它是工業最佳求解器，但已足以排除本輪「勝過強古典」主張。',
    '',
    '## 4. 稠密圖128變數探測',
    '',
    '\n'.join(largetable),
    '',
    '切割值越大越好，所有值回原圖計算。RQAOA靠近零尋找第一個峰、每gamma解析最佳beta，edge-only消去，到8變數窮舉。權重分布／規模接近文獻，但隨機種子、搜尋細節、實作和資料不是完整重現；沒有GUROBI最佳性證明。4例只能作可行性探測，不能推翻原論文或證明此方向普遍無效。',
    '',
    '## 哪個值得繼續',
    '',
    '1. **若一定要在RQAOA中選應用：先保留噴漆／固定順序的二元切換最小化。** 有真實目標、自然的成對限制，且在新例上相對專用貪婪仍有品質收益。參數轉移能降低重複調參成本，適合做可重現的資源－品質研究。這是已有研究的應用，需另找真正的領域缺口才能作創新。',
    '2. **稠密加權圖只列方法研究候選。** 本輪未得到穩定收益，且新例／128探測不支持直接拿文獻摘要當本機優勢。',
    '3. **SK與三正則挫折圖暫不優先。** 本輪小圖未勝強古典；尚未評估更高p、其他遞迴規則或真機，不能外推所有RQAOA變體。',
    '',
    '若競賽的必要條件是「實驗上勝過強古典」，本輪還沒有可據以選定的方向。若接受研究問題與負結果，噴漆可展示清楚的局部優點及成本邊界。不要改稱量子優勢，也不要只報RQ對greedy而隱藏tabu／SDP。',
    '',
    '## 驗證與證據邊界',
    '',
    '- 所有RQAOA均為p=1解析CPU模擬，未呼叫量子API、未跑真機、未計入硬體排隊或完整訓練shots。p=1相關可有效古典計算，故不是量子計算優勢證據。',
    '- 四類小圖的自寫QAOA分布與Qiskit statevector核對，誤差低於1e-15；噴漆映射對全部6變數配置逐項核對換色成本。',
    f'- 已獨立重算 {a["audit"]["rescored_outputs"]} 筆儲存解，最大原圖計分差 {a["audit"]["max_rescore_error"]:.3g}。SDP為SCS數值relaxation，不當成嚴格最佳性證書。',
    '- 每次搜尋與RSG轉換皆驗證成對異色可行性。大小案例、獨立種子與重複計時分開；沒有把不同規模的原始Ising能量相除作通用近似比。',
    '- 有界gamma掃描與第一峰策略都沒有每一殘餘圖的全域參數最優證書。最終挑選經過初篩，因此新40例才是主要確認證據；仍限於合成二元噴漆分布，沒有真實工廠容量／多色／可重排訂單限制。',
    '- 本輪只把噴漆模型中固定順序、每型兩次且兩次異色的問題稱為BPSP，不把它等同完整工業排程。',
    '',
    '## 可重跑資料與程式',
    '',
    f'主目錄：`{(ROOT/"research/rqaoa_opportunity").as_posix()}`。所有結果在`results/`，包含`discovery.csv`、`confirmation.csv`、`transfer.csv`、`fast_confirmation.csv`、`paint_canonical_RG.csv`、`dense128.csv`及完整K矩陣／解／縮圖軌跡。協定在各`*_protocol.json`，模型在`transfer_model.json`，聚合與配對bootstrap在`analysis.json`。',
    '',
    '重跑順序（沿用既有venv，設OPENBLAS_NUM_THREADS=1、OMP_NUM_THREADS=1、PYTHONIOENCODING=utf-8）：`sweep.py validate` → `sweep.py discovery` → `sweep.py confirmation` → `paint_baselines.py` → `transfer.py` → `fast_confirmation.py` → `dense_probe.py` → `analyze.py` → `report.py`。本輪同名輸出會更新，舊點雲與第9題資料不變。',
    '',
    '本次需求原話：',
    '',
    '> 幫我找找看有沒有其他 RQAOA 有優勢的地方，並且直接執行然後測試給我看看，哪一個比較有優勢？或是有機會、有優勢的地方',
    '',
    '以上是AI執行研究與判斷，未將結果或後續建議寫成Johnny已選定的題目。',
    ])
    path=ROOT/'output/RQAOA_其他優勢方向搜尋與實測_2026-11-12.md';path.write_text(report,encoding='utf-8')
    print(str(path));print('Figures',len(list(FIG.glob('*.png'))))

if __name__=='__main__':run()
